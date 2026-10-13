use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_std::Zero;

/// A multiplicative evaluation domain: a subgroup `H = <w>` of `F*` of order `n`, together with
/// forward and inverse transforms between coefficient and evaluation form.
///
/// This exists because the two APK configurations cannot share an FFT strategy. BW6-761's scalar
/// field has two-adicity 46, so radix-2 works for any power-of-two size. BW6-767's has two-adicity
/// **1**: no power-of-two domain beyond size 2 exists there, and no domain whose size is divisible
/// by 4 exists at all. Sizes there are divisors of `2 * 3^2 * 11 * 23 * 47 * 10177`, transformed
/// with mixed-radix Cooley-Tukey (and Rader for the 10177 factor).
///
/// Implementations must return evaluations in natural order, i.e. index `i` holds the value at
/// `w^i`. Several callers depend on that ordering.
pub trait FftDomain<F: PrimeField>: Clone + Sized {
    /// Domain size `n`.
    fn size(&self) -> usize;

    /// The generator `w`, a primitive `n`-th root of unity.
    fn generator(&self) -> F;

    /// `w^{-1}`.
    fn generator_inv(&self) -> F;

    /// `n^{-1}`, used to normalise the inverse transform.
    fn size_inv(&self) -> F;

    /// Coefficients -> evaluations at `(1, w, ..., w^{n-1})`.
    fn fft(&self, coeffs: &[F]) -> Vec<F>;

    /// Evaluations -> coefficients. Exact only when the underlying polynomial has degree `< n`;
    /// otherwise the result is the polynomial reduced modulo `X^n - 1`.
    fn interpolate(&self, evals: &[F]) -> Vec<F>;

    /// The `i`-th domain element, `w^i`.
    fn element(&self, i: usize) -> F {
        self.generator().pow([i as u64])
    }

    /// `Z_H(x) = x^n - 1`.
    ///
    /// Generic exponentiation rather than a chain of squarings: `n` is not a power of two in
    /// general.
    fn evaluate_vanishing_polynomial(&self, x: F) -> F {
        x.pow([self.size() as u64]) - F::one()
    }

    /// Divides by `Z_H(X) = X^n - 1`, returning `(quotient, remainder)`.
    ///
    /// Done in coefficient form rather than pointwise over a coset, which keeps this correct
    /// whether or not `H` is contained in the larger evaluation domain.
    ///
    /// Written out rather than delegated to `DenseOrSparsePolynomial::divide_with_q_and_r`,
    /// which since ark-poly 0.6 switches to Hensel division once the divisor's degree reaches
    /// 256 and multiplies by FFT to do it — so it panics with "field is not smooth enough to
    /// construct domain" on exactly the field this crate exists to support. Dividing by
    /// `X^n - 1` needs none of that: `X^i = X^(i-n) * (X^n - 1) + X^(i-n)`, so sweeping the
    /// coefficients from the top down carries each one into the quotient and into position
    /// `i - n`, in `O(deg p)` with no multiplications at all.
    fn divide_by_vanishing_poly(
        &self,
        poly: &DensePolynomial<F>,
    ) -> (DensePolynomial<F>, DensePolynomial<F>) {
        let n = self.size();
        let mut remainder = poly.coeffs.clone();
        if remainder.len() <= n {
            return (
                DensePolynomial::zero(),
                DensePolynomial::from_coefficients_vec(remainder),
            );
        }

        let mut quotient = vec![F::zero(); remainder.len() - n];
        // Top down, so that a coefficient carried into position `i - n` is itself carried
        // further when `i - n` is still at least `n` (i.e. when deg p >= 2n).
        for i in (n..remainder.len()).rev() {
            let c = remainder[i];
            if c.is_zero() {
                continue;
            }
            quotient[i - n] += c;
            remainder[i - n] += c;
        }
        remainder.truncate(n);

        (
            DensePolynomial::from_coefficients_vec(quotient),
            DensePolynomial::from_coefficients_vec(remainder),
        )
    }

    /// All Lagrange basis polynomials evaluated at `x`:
    /// `L_i(x) = (x^n - 1) / (n * (x - w^i)) * w^i`.
    fn evaluate_all_lagrange_coefficients(&self, x: F) -> Vec<F> {
        let n = self.size();
        let vanishing_scaled = self.evaluate_vanishing_polynomial(x) * self.size_inv();

        let mut denominators = Vec::with_capacity(n);
        let mut power = F::one();
        for _ in 0..n {
            denominators.push(x - power);
            power *= self.generator();
        }
        ark_ff::batch_inversion(&mut denominators);

        let mut result = Vec::with_capacity(n);
        power = F::one();
        for denom_inv in denominators {
            result.push(vanishing_scaled * denom_inv * power);
            power *= self.generator();
        }
        result
    }
}

/// Why a domain could not be built.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DomainError {
    /// The field has no evaluation domain this large. Its multiplicative group has order
    /// `q - 1`, so the available sizes are that number's divisors, and they run out.
    TooLarge { requested: usize },
    /// The field has no domain of *exactly* this size. `nearest` is the smallest one at least
    /// as large.
    NotExact { requested: usize, nearest: usize },
}

impl core::fmt::Display for DomainError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            DomainError::TooLarge { requested } => {
                write!(f, "no evaluation domain of size >= {} exists in this field", requested)
            }
            DomainError::NotExact { requested, nearest } => write!(
                f,
                "no evaluation domain of size exactly {} exists in this field; the smallest at least that large is {}",
                requested, nearest
            ),
        }
    }
}

impl std::error::Error for DomainError {}

/// Marks domains that can supply the sizes the `packed` scheme needs.
///
/// `packed` splits the bitmask into 256-bit chunks and asserts `256 | n`. Only a radix-2 domain
/// can guarantee that. Over BW6-767 it is outright impossible: `q - 1` carries a single factor
/// of 2, so no domain size there is even a multiple of 4, let alone 256.
///
/// This is a marker rather than a runtime check so that asking for a packed proof on a
/// configuration that cannot produce one fails to compile, instead of panicking inside the
/// prover after the caller has already built a keyset.
pub trait SupportsPackedScheme {}

/// The naive O(n^2) DFT: the polynomial with coefficients `coeffs`, evaluated at
/// `1, w, w^2, ..., w^(n-1)` by Horner's rule at each point.
///
/// Shared by `NaiveDomain` and the naive leaves of the Cooley-Tukey plan.
pub(crate) fn evaluate_at_powers<F: PrimeField>(coeffs: &[F], w: F, n: usize) -> Vec<F> {
    debug_assert!(
        coeffs.len() <= n,
        "{} coefficients for {} points",
        coeffs.len(),
        n
    );
    let mut result = Vec::with_capacity(n);
    let mut point = F::one();
    for _ in 0..n {
        result.push(
            coeffs
                .iter()
                .rev()
                .fold(F::zero(), |acc, &c| acc * point + c),
        );
        point *= w;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::{CooleyTukeyDomain, Radix2Domain};
    use ark_poly::univariate::DenseOrSparsePolynomial;
    use ark_std::{test_rng, One, UniformRand};

    /// The hand-written sweep must agree with arkworks' generic division, coefficient for
    /// coefficient, on a field where arkworks can actually run.
    ///
    /// Sizes are chosen either side of ark-poly's `SWITCH_TO_HENSEL_DIV` (256), so both of
    /// arkworks' own branches are the reference at least once.
    #[test]
    fn dividing_by_the_vanishing_polynomial_agrees_with_arkworks() {
        use ark_bw6_761::Fr;

        let rng = &mut test_rng();
        for n in [16usize, 128, 256, 512] {
            let domain = Radix2Domain::<Fr>::new(n);
            for degree in [0, n - 1, n, n + 1, 4 * n - 3, 5 * n] {
                let poly = DensePolynomial::from_coefficients_vec(
                    (0..=degree).map(|_| Fr::rand(rng)).collect(),
                );

                let mut vanishing = vec![Fr::zero(); n + 1];
                vanishing[0] = -Fr::one();
                vanishing[n] = Fr::one();
                let vanishing = DensePolynomial::from_coefficients_vec(vanishing);
                let expected = DenseOrSparsePolynomial::from(&poly)
                    .divide_with_q_and_r(&DenseOrSparsePolynomial::from(&vanishing))
                    .unwrap();

                assert_eq!(
                    domain.divide_by_vanishing_poly(&poly),
                    expected,
                    "n = {}, degree = {}",
                    n,
                    degree
                );
            }
        }
    }

    /// The regression. arkworks' division reaches for an FFT once the divisor's degree hits
    /// 256, which BW6-767 cannot supply, so every domain above 256 used to panic inside the
    /// prover. Nothing caught it because every APK-381 test ran at domain 253.
    #[test]
    fn dividing_by_the_vanishing_polynomial_works_on_a_non_smooth_field() {
        use ark_bw6_767::Fr;

        let rng = &mut test_rng();
        // 517 = 11 * 47, the first table entry past the threshold.
        let domain = CooleyTukeyDomain::<Fr>::new(517).unwrap();

        // Degree 4n - 3 is the highest the prover's constraint polynomial reaches, once
        // multiplied by its selector; the quotient is then degree 3n - 3.
        let quotient = DensePolynomial::from_coefficients_vec(
            (0..3 * 517 - 2).map(|_| Fr::rand(rng)).collect(),
        );
        let remainder =
            DensePolynomial::from_coefficients_vec((0..517).map(|_| Fr::rand(rng)).collect());

        // Reconstruct q * (X^n - 1) + r without multiplying polynomials, which is itself the
        // operation this field cannot do by FFT.
        let mut product = vec![Fr::zero(); quotient.coeffs.len() + 517];
        for (i, c) in quotient.coeffs.iter().enumerate() {
            product[i + 517] += c;
            product[i] -= c;
        }
        for (i, c) in remainder.coeffs.iter().enumerate() {
            product[i] += c;
        }
        let dividend = DensePolynomial::from_coefficients_vec(product);

        assert_eq!(
            domain.divide_by_vanishing_poly(&dividend),
            (quotient, remainder)
        );
    }

    /// A polynomial that vanishes on `H` divides exactly. This is the property the prover
    /// asserts on every proof.
    #[test]
    fn a_polynomial_vanishing_on_the_domain_leaves_no_remainder() {
        use ark_bw6_767::Fr;

        let rng = &mut test_rng();
        let domain = CooleyTukeyDomain::<Fr>::new(517).unwrap();
        let large = CooleyTukeyDomain::<Fr>::new(6 * 517).unwrap();

        // Interpolating evaluations that are zero over H gives a multiple of X^n - 1.
        let mut evals = vec![Fr::zero(); large.size()];
        for (i, e) in evals.iter_mut().enumerate() {
            if i % 6 != 0 {
                *e = Fr::rand(rng);
            }
        }
        let poly = DensePolynomial::from_coefficients_vec(large.interpolate(&evals));

        let (_, remainder) = domain.divide_by_vanishing_poly(&poly);
        assert!(remainder.is_zero());
    }
}
