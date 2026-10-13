use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::Polynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{end_timer, start_timer};
use w3f_pcs::pcs::PCS;

use crate::domain::DomainSet;
use crate::domains::{Domains, Evals};
use crate::piop::affine_addition::{AffineAdditionEvaluations, PartialSumsAndBitmaskCommitments};
use crate::piop::{
    RegisterCommitments, RegisterEvaluations, RegisterPolynomials, VerifierProtocol,
};
use crate::utils::LagrangeEvaluations;
use crate::{utils, AccumulatorSeed, ApkError, Bitmask};

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct BitmaskPackingCommitments<G: AffineRepr> {
    pub c_comm: G,
    pub acc_comm: G,
}

impl<G: AffineRepr> BitmaskPackingCommitments<G> {
    pub fn new(c_comm: G, acc_comm: G) -> Self {
        BitmaskPackingCommitments { c_comm, acc_comm }
    }
}

impl<G: AffineRepr> RegisterCommitments<G> for BitmaskPackingCommitments<G> {
    fn as_vec(&self) -> Vec<G> {
        vec![self.c_comm, self.acc_comm]
    }
}

#[derive(Clone)]
pub struct BitmaskPackingPolynomials<F: Field> {
    pub c_poly: DensePolynomial<F>,
    pub acc_poly: DensePolynomial<F>,
}

impl<F: Field> BitmaskPackingPolynomials<F> {
    //TODO: &self
    pub fn to_vec(self) -> Vec<DensePolynomial<F>> {
        vec![self.c_poly, self.acc_poly]
    }
}

impl<G: AffineRepr> RegisterPolynomials<G> for BitmaskPackingPolynomials<G::ScalarField> {
    type C = BitmaskPackingCommitments<G>;

    fn commit<E, F: Fn(&DensePolynomial<G::ScalarField>) -> Result<G, E>>(
        &self,
        f: F,
    ) -> Result<Self::C, E> {
        Ok(BitmaskPackingCommitments::<G>::new(
            f(&self.c_poly)?,
            f(&self.acc_poly)?,
        ))
    }
}

//TODO: remove pubs
#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct SuccinctAccountableRegisterEvaluations<F: FftField> {
    pub c: F,
    pub acc: F,
    pub basic_evaluations: AffineAdditionEvaluations<F>,
}

impl<F: FftField> RegisterEvaluations<F> for SuccinctAccountableRegisterEvaluations<F> {
    fn as_vec(&self) -> Vec<F> {
        let mut res = self.basic_evaluations.as_vec();
        res.extend(vec![self.c, self.acc]);
        res
    }
}

impl<F: PrimeField> SuccinctAccountableRegisterEvaluations<F> {
    pub fn evaluate_constraint_polynomials<IC, OC>(
        &self,
        apk: &IC::Affine,
        evals_at_zeta: &LagrangeEvaluations<F>,
        r: F,
        bitmask: &Bitmask,
        domain_size: u64,
    ) -> Result<Vec<F>, ApkError>
    where
        IC: AccumulatorSeed,
        OC: CurveGroup<ScalarField = F>,
        OC::ScalarField: From<IC::BaseField>,
    {
        let bits_in_bitmask_chunk = 256;
        let bits_in_big_int_limb = 64;
        assert_eq!(bits_in_bitmask_chunk % bits_in_big_int_limb, 0);
        let limbs_in_chunk = bits_in_bitmask_chunk / bits_in_big_int_limb;
        assert_eq!(domain_size % bits_in_bitmask_chunk, 0);
        let chunks_in_bitmask = domain_size / bits_in_bitmask_chunk; // TODO: bitmask should be right-padded with 0s to domain_size

        #[allow(clippy::expect_used, reason = "invariant argued in the message")]
        let bits_in_bitmask_chunk_inv = F::from(256u16)
            .inverse()
            .expect("invariant: 256 is nonzero in any field of odd characteristic");

        let powers_of_r = utils::powers(r, (chunks_in_bitmask - 1) as usize);
        #[allow(clippy::expect_used, reason = "invariant argued in the message")]
        let r_pow_m = r * powers_of_r
            .last()
            .expect("invariant: `powers` always returns at least r^0");
        let mut bitmask_chunks = bitmask.to_chunks_as_field_elements::<F>(limbs_in_chunk as usize);
        //TODO: pad in Bitmask
        bitmask_chunks.resize_with(chunks_in_bitmask as usize, || F::zero());
        assert_eq!(powers_of_r.len(), bitmask_chunks.len());
        let aggregated_bitmask = bitmask_chunks
            .into_iter()
            .zip(powers_of_r)
            .map(|(bj, rj)| bj * rj)
            .sum::<F>();

        let t_a_zeta_omega1 = start_timer!(|| "A(zw) as fraction");
        let zeta_omega_pow_m = evals_at_zeta.zeta_omega.pow([chunks_in_bitmask]); // m = chunks_in_bitmask
        let zeta_omega_pow_n = zeta_omega_pow_m.pow([bits_in_bitmask_chunk]); // n = domain_size
        let a_zeta_omega1 = bits_in_bitmask_chunk_inv * (zeta_omega_pow_n - F::one())
            / (zeta_omega_pow_m - F::one());
        end_timer!(t_a_zeta_omega1);

        let t_a_zeta_omega2 = start_timer!(|| "A(zw) as polynomial");
        let zeta_omega_pow_m = evals_at_zeta.zeta_omega.pow([chunks_in_bitmask]); // m = chunks_in_bitmask
        let a_zeta_omega2 = bits_in_bitmask_chunk_inv
            * utils::powers(zeta_omega_pow_m, (bits_in_bitmask_chunk - 1) as usize)
                .iter()
                .sum::<F>();
        end_timer!(t_a_zeta_omega2);

        assert_eq!(a_zeta_omega1, a_zeta_omega2);
        let two = F::from(2u8);
        let a = two + (r / two.pow([255u64]) - two) * a_zeta_omega1;

        let b = self.basic_evaluations.bitmask;
        let acc = self.acc;
        let c = self.c;

        let a6 = evaluate_inner_product_constraint_linearized(
            aggregated_bitmask,
            &evals_at_zeta,
            b,
            c,
            acc,
        );

        let a7 = evaluate_multipacking_mask_constraint_linearized(a, r_pow_m, &evals_at_zeta, c);

        let mut res = self
            .basic_evaluations
            .evaluate_constraint_polynomials::<IC, OC>(apk, evals_at_zeta)?;
        res.extend(vec![a6, a7]);
        Ok(res)
    }
}

impl<IC, OC, S> VerifierProtocol<IC, OC, S>
    for SuccinctAccountableRegisterEvaluations<OC::ScalarField>
where
    IC: AccumulatorSeed,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField>,
    S: PCS<OC::ScalarField>,
{
    type C1 = PartialSumsAndBitmaskCommitments<OC::Affine>;
    type C2 = BitmaskPackingCommitments<OC::Affine>;

    const POLYS_OPENED_AT_ZETA: usize = 8;

    fn restore_commitment_to_linearization_polynomial(
        &self,
        phi: OC::ScalarField,
        zeta_minus_omega_inv: OC::ScalarField,
        commitments: &PartialSumsAndBitmaskCommitments<OC::Affine>,
        extra_commitments: &BitmaskPackingCommitments<OC::Affine>,
    ) -> OC {
        let powers_of_phi = utils::powers(phi, 6);
        let mut r_comm = <AffineAdditionEvaluations<OC::ScalarField> as VerifierProtocol<
            IC,
            OC,
            S,
        >>::restore_commitment_to_linearization_polynomial(
            &self.basic_evaluations,
            phi,
            zeta_minus_omega_inv,
            &commitments.partial_sums,
            &(),
        );
        r_comm += extra_commitments.acc_comm * powers_of_phi[5];
        r_comm += extra_commitments.c_comm * powers_of_phi[6];
        r_comm
    }
}

// Constraint evaluators shared by prover and verifier. Free functions rather than inherent
// methods: the verifier evaluates constraints without ever constructing a domain, so it has
// no `D` to name.
pub(crate) fn evaluate_inner_product_constraint<F: PrimeField>(
    bitmask_chunks_aggregated: F,
    evals_at_zeta: &LagrangeEvaluations<F>,
    b_zeta: F,
    c_zeta: F,
    acc_zeta: F,
    acc_zeta_omega: F,
) -> F {
    acc_zeta_omega - acc_zeta - b_zeta * c_zeta + bitmask_chunks_aggregated * evals_at_zeta.l_last
}

pub(crate) fn evaluate_inner_product_constraint_linearized<F: PrimeField>(
    bitmask_chunks_aggregated: F,
    evals_at_zeta: &LagrangeEvaluations<F>,
    b_zeta: F,
    c_zeta: F,
    acc_zeta: F,
) -> F {
    evaluate_inner_product_constraint(
        bitmask_chunks_aggregated,
        evals_at_zeta,
        b_zeta,
        c_zeta,
        acc_zeta,
        F::zero(),
    )
}

pub(crate) fn evaluate_multipacking_mask_constraint<F: PrimeField>(
    a: F,
    r_pow_m: F,
    evals_at_zeta: &LagrangeEvaluations<F>,
    c_zeta: F,
    c_zeta_omega: F,
) -> F {
    c_zeta_omega - c_zeta * a - (F::one() - r_pow_m) * evals_at_zeta.l_last
}

pub(crate) fn evaluate_multipacking_mask_constraint_linearized<F: PrimeField>(
    a: F,
    r_pow_m: F,
    evals_at_zeta: &LagrangeEvaluations<F>,
    c_zeta: F,
) -> F {
    evaluate_multipacking_mask_constraint(a, r_pow_m, evals_at_zeta, c_zeta, F::zero())
}

pub(crate) struct BitmaskPackingRegisters<F: PrimeField, D: DomainSet<F>> {
    domains: Domains<F, D>,

    bitmask: Evals<F>,
    c: Evals<F>,
    c_shifted: Evals<F>,
    acc: Evals<F>,
    acc_shifted: Evals<F>,

    bitmask_chunks_aggregated: F,
    polynomials: BitmaskPackingPolynomials<F>,
    r: F,
}

impl<F: PrimeField, D: DomainSet<F>> BitmaskPackingRegisters<F, D> {
    // TODO: remove bitmask arg
    pub fn new(
        domains: Domains<F, D>,
        bitmask: &Bitmask,
        bitmask_chunks_aggregation_challenge: F, // denoted 'r' in the write-ups
    ) -> Self {
        let n = domains.size;
        let bits_in_bitmask_chunk = 256; //256 is the highest power of 2 that fits field bit capacity //TODO: const
        assert_eq!(n % bits_in_bitmask_chunk, 0); // n is a power of 2

        let mut bitmask = bitmask.to_bits_as_field_elements();
        bitmask.resize(domains.size, F::zero());

        let r = bitmask_chunks_aggregation_challenge;
        let c = Self::build_multipacking_mask_register(n, bits_in_bitmask_chunk, r);
        let acc = Self::build_partial_inner_products_register(n, &bitmask, &c);
        let bitmask_chunks_aggregated =
            bitmask.iter().zip(c.iter()).map(|(&b, c)| b * c).sum::<F>();

        let mut c_shifted = c.clone();
        c_shifted.rotate_left(1);
        let mut acc_shifted = acc.clone();
        acc_shifted.rotate_left(1);

        Self::new_unchecked(
            domains,
            bitmask,
            c,
            c_shifted,
            acc,
            acc_shifted,
            bitmask_chunks_aggregated,
            r,
        )
    }

    fn new_unchecked(
        domains: Domains<F, D>,

        bitmask: Vec<F>,
        c: Vec<F>,
        c_shifted: Vec<F>,
        acc: Vec<F>,
        acc_shifted: Vec<F>,
        bitmask_chunks_aggregated: F,
        r: F,
    ) -> Self {
        let c_polynomial = domains.interpolate(c);
        let acc_polynomial = domains.interpolate(acc);
        Self {
            domains: domains.clone(),

            bitmask: domains.amplify(bitmask),
            c: domains.amplify_polynomial(&c_polynomial),
            c_shifted: domains.amplify(c_shifted),
            acc: domains.amplify_polynomial(&acc_polynomial),
            acc_shifted: domains.amplify(acc_shifted),
            bitmask_chunks_aggregated,
            polynomials: BitmaskPackingPolynomials {
                c_poly: c_polynomial,
                acc_poly: acc_polynomial,
            },
            r,
        }
    }

    //TODO: comment
    fn build_multipacking_mask_register(
        domain_size: usize,
        chunk_size: usize,
        randomizer: F,
    ) -> Vec<F> {
        let powers_of_2 = utils::powers(F::from(2u8), chunk_size - 1);
        let powers_of_r = utils::powers(randomizer, domain_size / chunk_size - 1);
        // tensor product (powers_of_r X powers_of_2)
        powers_of_r
            .iter()
            .flat_map(|rj| powers_of_2.iter().map(move |_2k| *rj * _2k))
            .collect::<Vec<F>>()
    }

    /// Returns length n vec (0, a[0]b[0],...,a[n-2]b[n-2]), where n is domain size
    fn build_partial_inner_products_register(domain_size: usize, a: &Vec<F>, b: &Vec<F>) -> Vec<F> {
        // we ignore the last elements but still...
        assert_eq!(a.len(), domain_size);
        assert_eq!(b.len(), domain_size);
        let mut acc = Vec::with_capacity(domain_size);
        acc.push(F::zero());
        a.iter()
            .zip(b.iter())
            .map(|(a, b)| *a * b)
            .take(domain_size - 1)
            .for_each(|x| {
                #[allow(clippy::expect_used, reason = "invariant argued in the message")]
                let prev = *acc.last().expect("invariant: acc starts with one element");
                acc.push(x + prev);
            });
        acc
    }

    pub fn compute_inner_product_constraint_polynomial(&self) -> DensePolynomial<F> {
        let bc_ln_x4 = self
            .domains
            .l_last_scaled_by(self.bitmask_chunks_aggregated);
        let constraint =
            &(&(&self.acc_shifted - &self.acc) - &(&self.bitmask * &self.c)) + &bc_ln_x4;
        self.domains.interpolate_4x(&constraint)
    }

    pub fn compute_multipacking_mask_constraint_polynomial(&self) -> DensePolynomial<F> {
        let n = self.domains.size;
        let chunks = n / 256; //TODO: consts
        let mut a = vec![F::from(2u8); n];
        a.iter_mut()
            .step_by(256)
            .for_each(|a| *a = self.r / F::from(2u8).pow([255u64]));
        a.rotate_left(1);
        let a_x4 = self.domains.amplify(a);

        let x_todo = F::one() - self.r.pow([chunks as u64]); //TODO: name
        let ln_x4 = self.domains.l_last_scaled_by(x_todo);

        let a7 = &(&self.c_shifted - &(&self.c * &a_x4)) - &ln_x4;
        self.domains.interpolate_4x(&a7)
    }
}

impl<F: PrimeField, D: DomainSet<F>> BitmaskPackingRegisters<F, D> {
    pub fn evaluate_register_polynomials(&self, point: F) -> (F, F) {
        //TODO: struct
        (
            self.polynomials.c_poly.evaluate(&point),
            self.polynomials.acc_poly.evaluate(&point),
        )
    }

    pub fn compute_constraints_linearized(&self) -> Vec<DensePolynomial<F>> {
        vec![
            self.polynomials.acc_poly.clone(),
            self.polynomials.c_poly.clone(),
        ]
    }

    pub fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<F>> {
        let a6_poly = self.compute_inner_product_constraint_polynomial();
        let a7_poly = self.compute_multipacking_mask_constraint_polynomial();
        vec![a6_poly, a7_poly]
    }

    pub fn get_register_polynomials(&self) -> BitmaskPackingPolynomials<F> {
        self.polynomials.clone()
    }
}

#[cfg(test)]
mod tests {
    type TestDomain = crate::Radix2DomainSet<ark_bw6_761::Fr>;
    type TestDomains = crate::domains::Domains<ark_bw6_761::Fr, TestDomain>;
    use ark_bw6_761::Fr;
    use ark_ff::One;
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};

    use crate::test_helpers::_random_bits;

    use super::*;

    #[test]
    fn test_multipacking_mask_register() {
        let r = Fr::rand(&mut test_rng());
        let two = Fr::from(2u8);
        let multipacking_mask =
            BitmaskPackingRegisters::<Fr, TestDomain>::build_multipacking_mask_register(4, 2, r);
        assert_eq!(multipacking_mask, vec![Fr::one(), two, r, r * two]);
    }

    #[test]
    fn test_partial_inner_products_register() {
        let from_u8_vec = |v: [u8; 4]| v.iter().map(|&x| Fr::from(x)).collect::<Vec<Fr>>();
        let a = from_u8_vec([1, 2, 3, 4]);
        let b = from_u8_vec([5, 6, 7, 8]);
        let partial_inner_product =
            BitmaskPackingRegisters::<Fr, TestDomain>::build_partial_inner_products_register(
                4, &a, &b,
            );
        assert_eq!(
            partial_inner_product,
            from_u8_vec([0, 1 * 5, 1 * 5 + 2 * 6, 1 * 5 + 2 * 6 + 3 * 7])
        );
    }

    #[test]
    fn test_multipacking_mask_constraint() {
        let rng = &mut test_rng();
        let n = 256;
        let m = n - 1;
        let domains = TestDomains::new(n);

        let bitmask = Bitmask::from_bits(&_random_bits(m, 0.5, rng));

        let r = Fr::rand(rng);
        let acc_registers = BitmaskPackingRegisters::new(domains.clone(), &bitmask, r);
        let constraint_poly = acc_registers.compute_multipacking_mask_constraint_polynomial();
        assert_eq!(constraint_poly.degree(), 2 * n - 2);
        assert!(domains.is_zero(&constraint_poly));
    }
}
