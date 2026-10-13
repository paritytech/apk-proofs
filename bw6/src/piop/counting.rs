use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{FftField, PrimeField};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use w3f_pcs::pcs::PCS;

use crate::domain::DomainSet;
use crate::domains::Domains;
use crate::piop::affine_addition::{
    AffineAdditionEvaluations, AffineAdditionRegisters, PartialSumsAndBitmaskCommitments,
    PartialSumsAndBitmaskPolynomials,
};
use crate::piop::bit_counting::{BitCountingEvaluation, BitCountingRegisters};
use crate::piop::{
    ProverProtocol, RegisterCommitments, RegisterEvaluations, RegisterPolynomials, VerifierProtocol,
};
use crate::utils::LagrangeEvaluations;
use crate::{utils, AccumulatorSeed, ApkError, Bitmask, CountingPublicInput, Keyset};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct CountingCommitments<G: AffineRepr> {
    affine_addition_commitments: PartialSumsAndBitmaskCommitments<G>,
    partial_counts_commitment: G,
}

impl<G: AffineRepr> RegisterCommitments<G> for CountingCommitments<G> {
    fn as_vec(&self) -> Vec<G> {
        let mut commitments = self.affine_addition_commitments.as_vec();
        commitments.push(self.partial_counts_commitment);
        commitments
    }
}

pub struct CountingPolynomials<F: FftField> {
    affine_addition_polynomials: PartialSumsAndBitmaskPolynomials<F>,
    partial_counts_polynomial: DensePolynomial<F>,
}

impl<G: AffineRepr> RegisterPolynomials<G> for CountingPolynomials<G::ScalarField> {
    type C = CountingCommitments<G>;

    fn commit<E, FN: Clone + Fn(&DensePolynomial<G::ScalarField>) -> Result<G, E>>(
        &self,
        f: FN,
    ) -> Result<Self::C, E> {
        Ok(CountingCommitments {
            affine_addition_commitments: self.affine_addition_polynomials.commit(f.clone())?,
            partial_counts_commitment: f(&self.partial_counts_polynomial)?,
        })
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct CountingEvaluations<F: FftField> {
    affine_addition_evaluations: AffineAdditionEvaluations<F>,
    partial_counts_evaluation: BitCountingEvaluation<F>,
}

impl<F: FftField> RegisterEvaluations<F> for CountingEvaluations<F> {
    fn as_vec(&self) -> Vec<F> {
        let mut evals = self.affine_addition_evaluations.as_vec();
        evals.push(self.partial_counts_evaluation.0);
        evals
    }
}

pub struct CountingScheme<F: PrimeField, D: DomainSet<F>> {
    affine_addition_registers: AffineAdditionRegisters<F, D>,
    bit_counting_registers: BitCountingRegisters<F, D>,
    register_evaluations: Option<CountingEvaluations<F>>,
}

impl<IC, OC, S, D> ProverProtocol<IC, OC, S, D> for CountingScheme<OC::ScalarField, D>
where
    IC: AccumulatorSeed,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField> + FftField,
    S: PCS<OC::ScalarField>,
    D: DomainSet<OC::ScalarField>,
{
    type P1 = CountingPolynomials<OC::ScalarField>;
    type P2 = ();
    type E = CountingEvaluations<OC::ScalarField>;
    type PI = CountingPublicInput<IC>;

    fn init(
        domains: Domains<OC::ScalarField, D>,
        bitmask: Bitmask,
        keyset: Keyset<IC, OC, D>,
    ) -> Self {
        CountingScheme {
            affine_addition_registers: AffineAdditionRegisters::new(
                domains.clone(),
                keyset,
                &bitmask.to_bits(),
            ),
            bit_counting_registers: BitCountingRegisters::new(domains, &bitmask),
            register_evaluations: None,
        }
    }

    fn get_register_polynomials_to_commit1(&self) -> Self::P1 {
        CountingPolynomials {
            affine_addition_polynomials: self
                .affine_addition_registers
                .get_partial_sums_and_bitmask_polynomials(),
            partial_counts_polynomial: self.bit_counting_registers.get_partial_counts_polynomial(),
        }
    }

    fn get_register_polynomials_to_commit2(
        &mut self,
        _verifier_challenge: OC::ScalarField,
    ) -> Self::P2 {
        ()
    }

    fn get_register_polynomials_to_open(self) -> Vec<DensePolynomial<OC::ScalarField>> {
        [
            self.affine_addition_registers
                .get_register_polynomials()
                .to_vec(),
            vec![self.bit_counting_registers.get_partial_counts_polynomial()],
        ]
        .concat()
    }

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<OC::ScalarField>> {
        [
            self.affine_addition_registers
                .compute_constraint_polynomials::<IC, OC>(),
            self.bit_counting_registers.constraints(),
        ]
        .concat()
    }

    fn evaluate_register_polynomials(&mut self, point: OC::ScalarField) -> Self::E {
        let affine_addition_evaluations = self
            .affine_addition_registers
            .evaluate_register_polynomials(point);
        let partial_counts_evaluation = self
            .bit_counting_registers
            .evaluate_partial_counts_register(point);
        let evals = CountingEvaluations {
            affine_addition_evaluations,
            partial_counts_evaluation,
        };
        self.register_evaluations = Some(evals.clone());
        evals
    }

    fn compute_linearization_polynomial(
        &self,
        phi: OC::ScalarField,
        zeta: OC::ScalarField,
    ) -> DensePolynomial<OC::ScalarField> {
        #[allow(clippy::expect_used, reason = "see ROUNDS_IN_ORDER")]
        let evals = self
            .register_evaluations
            .as_ref()
            .expect(crate::piop::ROUNDS_IN_ORDER);
        let parts = [
            self.affine_addition_registers
                .compute_constraints_linearized(&evals.affine_addition_evaluations, zeta),
            self.bit_counting_registers.constraints_lin(),
        ]
        .concat();
        utils::randomize(phi, &parts)
    }
}

impl<IC, OC, S> VerifierProtocol<IC, OC, S> for CountingEvaluations<OC::ScalarField>
where
    IC: AccumulatorSeed,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField> + FftField,
    S: PCS<OC::ScalarField>,
{
    type C1 = CountingCommitments<OC::Affine>;
    type C2 = ();

    const POLYS_OPENED_AT_ZETA: usize = 7;

    fn restore_commitment_to_linearization_polynomial(
        &self,
        phi: OC::ScalarField,
        zeta_minus_omega_inv: OC::ScalarField,
        commitments: &Self::C1,
        _extra_commitments: &Self::C2,
    ) -> OC {
        let powers_of_phi = utils::powers(phi, 6);
        let partial_sums_commitments = &commitments.affine_addition_commitments.partial_sums;
        let mut r_comm = <AffineAdditionEvaluations<OC::ScalarField> as VerifierProtocol<
            IC,
            OC,
            S,
        >>::restore_commitment_to_linearization_polynomial(
            &self.affine_addition_evaluations,
            phi,
            zeta_minus_omega_inv,
            partial_sums_commitments,
            &(),
        );
        r_comm += commitments.partial_counts_commitment * powers_of_phi[5];
        r_comm
    }
}

impl<F: FftField> CountingEvaluations<F> {
    pub fn evaluate_constraint_polynomials<IC, OC>(
        &self,
        apk: IC::Affine,
        count: OC::ScalarField,
        evals_at_zeta: &LagrangeEvaluations<OC::ScalarField>,
    ) -> Result<Vec<OC::ScalarField>, ApkError>
    where
        IC: AccumulatorSeed,
        OC: CurveGroup<ScalarField = F>,
        F: From<IC::BaseField>,
    {
        let b_at_zeta = self.affine_addition_evaluations.bitmask;
        Ok([
            self.affine_addition_evaluations
                .evaluate_constraint_polynomials::<IC, OC>(&apk, evals_at_zeta)?,
            self.partial_counts_evaluation.evaluate_constraints_at_zeta(
                count,
                b_at_zeta,
                evals_at_zeta.l_last,
            ),
        ]
        .concat())
    }
}

#[cfg(test)]
mod tests {
    type TestDomain = crate::Radix2DomainSet<ark_bw6_761::Fr>;
    type TestDomains = crate::domains::Domains<ark_bw6_761::Fr, TestDomain>;
    use super::*;
    use crate::instances::bls12_377_bw6_761::kzg::Pcs;
    use crate::test_helpers::{_random_bits, random_pks};
    use ark_bls12_377::G1Projective;
    use ark_bw6_761::{Fr, G1Projective as OuterCurve};
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};
    use w3f_pcs::pcs::PcsParams;

    #[test]
    fn test_polynomial_ordering() {
        let rng = &mut test_rng();
        let n = 16;
        let m = n - 1;

        let pcs_params = Pcs::setup(m, rng);
        let mut keyset =
            Keyset::<G1Projective, OuterCurve, TestDomain>::new(random_pks(m, rng)).unwrap();
        keyset.amplify();

        let mut scheme: CountingScheme<Fr, TestDomain> =
            ProverProtocol::<G1Projective, OuterCurve, Pcs, TestDomain>::init(
                TestDomains::new(n),
                Bitmask::from_bits(&_random_bits(m, 0.5, rng)),
                keyset,
            );

        let zeta = Fr::rand(rng);

        let actual_commitments = <CountingScheme<Fr, TestDomain> as ProverProtocol<
            G1Projective,
            OuterCurve,
            Pcs,
            TestDomain,
        >>::get_register_polynomials_to_commit1(&scheme)
        .commit(|p| Pcs::commit(&pcs_params.ck(), &p).map(|c| c.0))
        .unwrap()
        .as_vec();
        let actual_evaluations = <CountingScheme<Fr, TestDomain> as ProverProtocol<
            G1Projective,
            OuterCurve,
            Pcs,
            TestDomain,
        >>::evaluate_register_polynomials(&mut scheme, zeta)
        .as_vec();
        let polynomials = <CountingScheme<Fr, TestDomain> as ProverProtocol<
            G1Projective,
            OuterCurve,
            Pcs,
            TestDomain,
        >>::get_register_polynomials_to_open(scheme);

        let expected_evaluations = polynomials
            .iter()
            .map(|p| p.evaluate(&zeta))
            .collect::<Vec<_>>();
        assert_eq!(actual_evaluations, expected_evaluations);

        let expected_commitments = polynomials
            .iter()
            .skip(2) // keyset commitment is publicly known
            .map(|p| Pcs::commit(&pcs_params.ck(), &p).unwrap().0)
            .collect::<Vec<_>>();
        assert_eq!(actual_commitments, expected_commitments);
    }
}
