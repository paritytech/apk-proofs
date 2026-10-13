//! The O(n^2) reference domain. Compiled for tests and under `test-utils` only: nothing on
//! the prover's or verifier's path uses it.

use ark_ff::PrimeField;

use super::subgroup::subgroup_generator;
use super::types::{evaluate_at_powers, FftDomain};

/// An evaluation domain that transforms by naive O(n^2) evaluation.
///
/// This works for any `n` dividing `q - 1`, with no smoothness requirement at all, which makes it
/// the reference implementation: the Cooley-Tukey and Rader domains are differential-tested
/// against it, and benchmarked against it in `benches/`.
///
/// It is not the domain to prove with. At n = 1551 a single transform is on the order of a
/// million field multiplications.
#[derive(Clone, Debug)]
pub struct NaiveDomain<F: PrimeField> {
    size: usize,
    w: F,
    w_inv: F,
    size_inv: F,
}

impl<F: PrimeField> NaiveDomain<F> {
    /// Returns `None` if `F*` has no subgroup of order `size`, or if `size` is not invertible.
    pub fn new(size: usize) -> Option<Self> {
        let w = subgroup_generator::<F>(size)?;
        let w_inv = w.inverse()?;
        let size_inv = F::from(size as u64).inverse()?;
        debug_assert!(
            w.pow([size as u64]).is_one(),
            "generator has the wrong order"
        );
        Some(NaiveDomain {
            size,
            w,
            w_inv,
            size_inv,
        })
    }
}

impl<F: PrimeField> FftDomain<F> for NaiveDomain<F> {
    fn size(&self) -> usize {
        self.size
    }

    fn generator(&self) -> F {
        self.w
    }

    fn generator_inv(&self) -> F {
        self.w_inv
    }

    fn size_inv(&self) -> F {
        self.size_inv
    }

    /// `coeffs` shorter than the domain is fine and means the high coefficients are zero;
    /// longer would alias modulo `X^size - 1`, so it is rejected.
    fn fft(&self, coeffs: &[F]) -> Vec<F> {
        assert!(
            coeffs.len() <= self.size,
            "{} coefficients do not fit a domain of size {}",
            coeffs.len(),
            self.size
        );
        evaluate_at_powers(coeffs, self.w, self.size)
    }

    fn interpolate(&self, evals: &[F]) -> Vec<F> {
        assert_eq!(
            evals.len(),
            self.size,
            "interpolation needs exactly `size` evaluations"
        );
        let mut coeffs = evaluate_at_powers(evals, self.w_inv, self.size);
        coeffs.iter_mut().for_each(|c| *c *= self.size_inv);
        coeffs
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Radix2Domain;
    use ark_std::{test_rng, UniformRand, Zero};

    // BW6-761's scalar field: two-adicity 46, so radix-2 sizes exist and can be cross-checked.
    type Fr761 = ark_bw6_761::Fr;
    // BW6-767's scalar field (= BLS12-381's base field): two-adicity 1.
    type Fr767 = ark_bw6_767::Fr;

    #[test]
    fn roundtrips_on_a_non_power_of_two_domain() {
        let rng = &mut test_rng();
        // 11 * 23, a size no radix-2 domain can represent.
        let domain = NaiveDomain::<Fr767>::new(253).expect("253 divides q-1");
        let evals: Vec<Fr767> = (0..domain.size()).map(|_| Fr767::rand(rng)).collect();

        let coeffs = domain.interpolate(&evals);
        assert_eq!(domain.fft(&coeffs), evals);
    }

    /// The whole point of this domain: it agrees with arkworks wherever both are defined.
    /// If the generators disagreed, every downstream evaluation would silently differ.
    #[test]
    fn agrees_with_radix2_where_both_exist() {
        let rng = &mut test_rng();
        let n = 256;
        let naive = NaiveDomain::<Fr761>::new(n).unwrap();
        let radix2 = Radix2Domain::<Fr761>::new(n);

        assert_eq!(naive.generator(), radix2.generator());
        assert_eq!(naive.generator_inv(), radix2.generator_inv());
        assert_eq!(naive.size_inv(), radix2.size_inv());

        let coeffs: Vec<Fr761> = (0..n).map(|_| Fr761::rand(rng)).collect();
        assert_eq!(naive.fft(&coeffs), radix2.fft(&coeffs));

        let evals: Vec<Fr761> = (0..n).map(|_| Fr761::rand(rng)).collect();
        assert_eq!(naive.interpolate(&evals), radix2.interpolate(&evals));

        let z = Fr761::rand(rng);
        assert_eq!(
            naive.evaluate_vanishing_polynomial(z),
            radix2.evaluate_vanishing_polynomial(z)
        );
        assert_eq!(
            naive.evaluate_all_lagrange_coefficients(z),
            radix2.evaluate_all_lagrange_coefficients(z)
        );
    }

    /// A polynomial of degree < n amplified to a larger domain must evaluate to the same values.
    /// This is the pattern `Domains::_amplify` relies on, and it is exactly what frees APK-381
    /// from needing nested domains: the pair below is deliberately NOT nested, since
    /// 1551 / 253 is not an integer, and amplification still round-trips.
    #[test]
    fn amplifying_to_a_larger_domain_preserves_the_polynomial() {
        let rng = &mut test_rng();
        let small = NaiveDomain::<Fr767>::new(253).unwrap(); // 11 * 23
        let large = NaiveDomain::<Fr767>::new(1551).unwrap(); // 3 * 11 * 47, and 1551 >= 4*253 - 3
        assert_ne!(
            large.size() % small.size(),
            0,
            "the pair must not be nested"
        );

        let evals: Vec<Fr767> = (0..small.size()).map(|_| Fr767::rand(rng)).collect();
        let coeffs = small.interpolate(&evals);
        let amplified = large.fft(&coeffs);

        // Interpolating back over the large domain recovers the same polynomial, padded.
        let recovered = large.interpolate(&amplified);
        assert_eq!(&recovered[..coeffs.len()], &coeffs[..]);
        assert!(recovered[coeffs.len()..].iter().all(|c| c.is_zero()));
    }

    #[test]
    fn rejects_sizes_that_are_not_subgroup_orders() {
        // q - 1 = 2 * 3^2 * 11 * 23 * 47 * 10177 * 859267 * 52437899 * M for BW6-767's Fr,
        // so 4 does not divide it: two-adicity is 1.
        assert!(NaiveDomain::<Fr767>::new(4).is_none());
        assert!(NaiveDomain::<Fr767>::new(256).is_none());
        assert!(NaiveDomain::<Fr767>::new(5).is_none());
        // ...but these do.
        assert!(NaiveDomain::<Fr767>::new(2).is_some());
        assert!(NaiveDomain::<Fr767>::new(9306).is_some());
    }

    /// Two-adicity 1 is the constraint the whole APK-381 design turns on. Pin it.
    #[test]
    fn bw6_767_scalar_field_has_two_adicity_one() {
        use ark_ff::FftField;
        assert_eq!(Fr767::TWO_ADICITY, 1);
        assert_eq!(Fr761::TWO_ADICITY, 46);
    }
}
