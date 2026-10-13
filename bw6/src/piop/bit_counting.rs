use ark_ff::{FftField, Field, PrimeField, Zero};
use ark_poly::univariate::DensePolynomial;
use ark_poly::Polynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::iter::once;

use crate::domain::DomainSet;
use crate::domains::Domains;
use crate::Bitmask;

// This "gadget" is used in the 'counting' scheme to constraint the number of set bits in the bitmask.

// It adds a single register tracking partial bitmask sums:
// r[0] = 0, r[i] = b[0] + ... + b[i-1], i = 1,...,n-1,
// and a single constraint operating on this register and the bitmask register b.
// Let S = b[0]+...+b[n-1] - the value we want to constraint and
// let SL[i] = 0, i = 0,...,n-2, and SL[n-1] = S.
// The constraint then is
// r[i+1] = r[i] + b[i] - SL[i], i = 0,...,n-1
// For i = 0,...,n-2 it ensures that the register is well-formed: r[i+1] = r[i] + b[i],
// for i = n-1, r[n] = r[n-1] + b[n-1] - S, but r[n] = r[0] and r[n-1] + b[n-1] = S.
// It follows that r[0] = 0 and S is properly constrained.

// After interpolation SL this defined becomes S * L_{n-1}(Z) that is known by the verifier, and easy to compute.

// What makes this constraint different from any other in the codebase is that it is linear in the registers involved,
// so by construction deg(C) <= n-1, and C(w^i) = 0 for i = 0,...,n-1, that results in C = 0 identically.
// To check this constraint holds, the verifier doesn't need to check C(Z) = Q(Z)(Z^n - 1), it could just check C(zeta) = 0,
// But we use the former check not to handle this case differently.

pub(crate) struct BitCountingRegisters<F: PrimeField, D: DomainSet<F>> {
    domains: Domains<F, D>,
    bitmask: Vec<F>,
    partial_counts: DensePolynomial<F>,
}

impl<F: PrimeField, D: DomainSet<F>> BitCountingRegisters<F, D> {
    pub fn new(domains: Domains<F, D>, bitmask: &Bitmask) -> Self {
        let mut bitmask = bitmask.to_bits_as_field_elements();
        bitmask.resize_with(domains.size, || F::zero());
        let partial_counts = Self::build_partial_counts_register(&bitmask);
        Self::new_unchecked(domains, bitmask, partial_counts)
    }

    fn new_unchecked(domains: Domains<F, D>, bitmask: Vec<F>, partial_counts: Vec<F>) -> Self {
        let partial_counts = domains.interpolate(partial_counts);
        Self {
            domains,
            bitmask,
            partial_counts,
        }
    }

    /// Returns length n vec (0, b[0], b[0] + b[1], ..., b[0] + b[1] + b[n-2])
    fn build_partial_counts_register(bitmask: &[F]) -> Vec<F> {
        let partial_counts = bitmask.iter().scan(F::zero(), |state, bit| {
            *state += bit;
            Some(*state)
        });
        once(F::zero())
            .chain(partial_counts.take(bitmask.len() - 1))
            .collect()
    }

    pub fn constraints(&self) -> Vec<DensePolynomial<F>> {
        vec![
            BitCount::constraint_poly::<F>(),
            BitmaskEndsWithZero::constraint_poly(self),
        ]
    }

    pub fn constraints_lin(&self) -> Vec<DensePolynomial<F>> {
        vec![
            BitCount::linearization(&self),
            BitmaskEndsWithZero::linearization::<F>(),
        ]
    }

    pub fn get_partial_counts_polynomial(&self) -> DensePolynomial<F> {
        self.partial_counts.clone()
    }

    pub fn evaluate_partial_counts_register(&self, zeta: F) -> BitCountingEvaluation<F> {
        let eval = self.partial_counts.evaluate(&zeta);
        BitCountingEvaluation(eval)
    }

    #[cfg(test)]
    fn get_bitmask_polynomial(&self) -> DensePolynomial<F> {
        self.domains.interpolate(self.bitmask.clone())
    }
}

struct BitmaskEndsWithZero {}

// Constraints the last bitmask element to 0
// by showing that the polynomial b(X) * L_{n-1}(X) is zero over the domain.
impl BitmaskEndsWithZero {
    // C = b * L_{n-1}
    fn constraint_poly<F: PrimeField, D: DomainSet<F>>(
        registers: &BitCountingRegisters<F, D>,
    ) -> DensePolynomial<F> {
        let n = registers.domains.size;
        let mut ln = vec![F::zero(); n];
        ln[n - 1] = F::one();
        let ln_x2 = registers.domains.amplify_x2(ln);
        let b_x2 = registers.domains.amplify_x2(registers.bitmask.clone());
        let c = &b_x2 * &ln_x2;
        registers.domains.interpolate_2x(&c)
    }

    fn linearization<F: FftField>() -> DensePolynomial<F> {
        DensePolynomial::zero()
    }

    // C(z) = b(z) * L_{n-1}(z)
    fn _evaluate_full<F: FftField>(bitmask_at_zeta: F, l_last_at_zeta: F) -> F {
        bitmask_at_zeta * l_last_at_zeta
    }

    fn evaluate_main<F: FftField>(bitmask_at_zeta: F, l_last_at_zeta: F) -> F {
        Self::_evaluate_full::<F>(bitmask_at_zeta, l_last_at_zeta)
    }
}

struct BitCount {}

// partial_counts(wZ) - partial_counts(Z) - bitmask(Z) + count*L_{n-1}(Z)
impl BitCount {
    fn constraint_poly<F: Field>() -> DensePolynomial<F> {
        DensePolynomial::zero()
    }

    // Though the constraint is zero, the verifier still needs the opening of the register in zeta * omega.
    fn linearization<F: PrimeField, D: DomainSet<F>>(
        registers: &BitCountingRegisters<F, D>,
    ) -> DensePolynomial<F> {
        registers.get_partial_counts_polynomial()
    }

    fn _evaluate_full<F: FftField>(
        eval: &BitCountingEvaluation<F>,
        count: F,
        partial_counts_at_zeta_omega: F,
        bitmask_at_zeta: F,
        l_last_at_zeta: F,
    ) -> F {
        partial_counts_at_zeta_omega - eval.0 - bitmask_at_zeta + count * l_last_at_zeta
    }

    fn evaluate_main<F: FftField>(
        eval: &BitCountingEvaluation<F>,
        count: F,
        bitmask_at_zeta: F,
        l_last_at_zeta: F,
    ) -> F {
        Self::_evaluate_full::<F>(eval, count, F::zero(), bitmask_at_zeta, l_last_at_zeta)
    }
}

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct BitCountingEvaluation<F: FftField>(pub F);

impl<F: FftField> BitCountingEvaluation<F> {
    pub fn evaluate_constraints_at_zeta(
        &self,
        count: F,
        bitmask_at_zeta: F,
        l_last_at_zeta: F,
    ) -> Vec<F> {
        vec![
            BitCount::evaluate_main(&self, count, bitmask_at_zeta, l_last_at_zeta),
            BitmaskEndsWithZero::evaluate_main::<F>(bitmask_at_zeta, l_last_at_zeta),
        ]
    }
}

#[cfg(test)]
mod tests {
    type TestDomain = crate::Radix2DomainSet<ark_bw6_761::Fr>;
    type TestDomains = crate::domains::Domains<ark_bw6_761::Fr, TestDomain>;
    use crate::domain::FftDomain;
    use crate::test_helpers::{_random_bitmask, _random_bits};
    use crate::utils;
    use crate::utils::lagrange_evaluations;
    use ark_bw6_761::Fr;
    use ark_poly::Polynomial;
    use ark_std::{test_rng, UniformRand};

    use super::*;

    #[test]
    fn test_partial_counts_register() {
        let rng = &mut test_rng();
        let n = 16;
        let bitmask = _random_bitmask::<_, ark_bw6_761::G1Projective>(n, rng);
        let partial_counts =
            BitCountingRegisters::<Fr, TestDomain>::build_partial_counts_register(&bitmask);

        assert_eq!(partial_counts.len(), 16);
        assert_eq!(partial_counts[0], Fr::zero());
        for i in 1..n - 1 {
            assert_eq!(partial_counts[i], bitmask[..i].iter().sum());
        }
    }

    #[test]
    fn test_bit_counting_constraint() {
        let rng = &mut test_rng();
        let n = 16;
        let domains = TestDomains::new(n);

        let bitmask = Bitmask::from_bits(&_random_bits(n, 2.0 / 3.0, rng));
        let count = Fr::from(bitmask.count_ones() as u32);
        let registers = BitCountingRegisters::<Fr, TestDomain>::new(domains.clone(), &bitmask);

        let z = Fr::rand(rng);
        let w = domains.omega;

        let acc_z = registers.evaluate_partial_counts_register(z);
        let acc_zw = registers.get_partial_counts_polynomial().evaluate(&(z * w));
        let bitmask_z = registers.get_bitmask_polynomial().evaluate(&z);
        let domain_z = utils::lagrange_evaluations(z, domains.domain());

        let x_full = BitCount::constraint_poly::<Fr>();
        let x_lin = BitCount::linearization(&registers);
        let x_eval_full =
            BitCount::_evaluate_full(&acc_z, count, acc_zw, bitmask_z, domain_z.l_last);
        let x_eval_main = BitCount::evaluate_main(&acc_z, count, bitmask_z, domain_z.l_last);
        let x_lin_zw = x_lin.evaluate(&(z * w));

        assert_eq!(x_eval_full, x_full.evaluate(&z));
        assert_eq!(x_eval_full, x_eval_main + x_lin_zw);
        assert!(domains
            .domain()
            .divide_by_vanishing_poly(&x_full)
            .1
            .is_zero()); // actually x_full is 0 over the field
    }

    #[test]
    fn test_bitmask_ends_with_zero_constraint() {
        let rng = &mut test_rng();
        let n = 16;
        let domains = TestDomains::new(n);
        let domain = domains.domain().clone();

        let bits = _random_bits(n, 2.0 / 3.0, rng);

        let mut good_bitmask = bits.clone();
        good_bitmask[n - 1] = false;
        let good_bitmask = Bitmask::from_bits(&good_bitmask);
        let registers = BitCountingRegisters::<Fr, TestDomain>::new(domains.clone(), &good_bitmask);
        let constraint = BitmaskEndsWithZero::constraint_poly(&registers);
        assert_eq!(constraint.degree(), 2 * (n - 1));
        assert!(domain.divide_by_vanishing_poly(&constraint).1.is_zero());

        let zeta = Fr::rand(rng);
        let prover_eval = constraint.evaluate(&zeta);
        let domain_evals = lagrange_evaluations(zeta, &domain);
        let bitmask_at_zeta = registers.get_bitmask_polynomial().evaluate(&zeta);
        let verifier_eval =
            BitmaskEndsWithZero::_evaluate_full::<Fr>(bitmask_at_zeta, domain_evals.l_last);
        assert_eq!(prover_eval, verifier_eval);

        let mut bad_bitmask = bits.clone();
        bad_bitmask[n - 1] = true;
        let bad_bitmask = Bitmask::from_bits(&bad_bitmask);
        let registers = BitCountingRegisters::<Fr, TestDomain>::new(domains.clone(), &bad_bitmask);
        let constraint = BitmaskEndsWithZero::constraint_poly(&registers);
        assert!(!domain.divide_by_vanishing_poly(&constraint).1.is_zero());
    }
}
