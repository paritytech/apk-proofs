use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::DenseUVPolynomial;
use ark_std::ops::{Add, AddAssign, Deref, Mul, MulAssign, Sub, SubAssign};
use ark_std::Zero;

use crate::domain::{DomainSet, FftDomain};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Evaluations of a polynomial over some domain, in natural order.
///
/// Replaces arkworks' `Evaluations`, which is tied to `EvaluationDomain` and so cannot be used
/// once the domain type is abstract. The domain is deliberately *not* carried in the type: the
/// length assert on every pointwise operation is what catches a domain mix-up.
#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Evals<F> {
    pub evals: Vec<F>,
}

impl<F> Deref for Evals<F> {
    type Target = Vec<F>;
    fn deref(&self) -> &Self::Target {
        &self.evals
    }
}

impl<F> From<Vec<F>> for Evals<F> {
    fn from(evals: Vec<F>) -> Self {
        Evals { evals }
    }
}

macro_rules! impl_pointwise_op {
    ($trait:ident, $method:ident, $assign_trait:ident, $assign_method:ident, $op:tt) => {
        impl<'a, F: PrimeField> $assign_trait<&'a Evals<F>> for Evals<F> {
            #[inline]
            fn $assign_method(&mut self, rhs: &'a Evals<F>) {
                assert_eq!(self.evals.len(), rhs.evals.len(), "evaluations over different domains");
                ark_std::cfg_iter_mut!(self.evals)
                    .zip(ark_std::cfg_iter!(rhs.evals))
                    .for_each(|(a, b)| *a $op b)
            }
        }

        impl<'a, 'b, F: PrimeField> $trait<&'a Evals<F>> for &'b Evals<F> {
            type Output = Evals<F>;

            #[inline]
            fn $method(self, rhs: &'a Evals<F>) -> Evals<F> {
                let mut result = self.clone();
                result.$assign_method(rhs);
                result
            }
        }
    };
}

impl_pointwise_op!(Add, add, AddAssign, add_assign, +=);
impl_pointwise_op!(Sub, sub, SubAssign, sub_assign, -=);
impl_pointwise_op!(Mul, mul, MulAssign, mul_assign, *=);
/// The three domains the PIOP works over, plus the precomputations that depend on all three.
///
/// The domains themselves come from a [`DomainSet`], which chooses the trace domain `H` of size
/// `n` together with the two larger ones that hold products of register polynomials in
/// evaluation form. The largest product interpolated back has degree `4n - 4`, so the large
/// domain needs at least `4n - 3` points; see [`crate::DomainTriple`].
///
/// The names `domain2x` / `domain4x` are historical. On APK-377 they really are `2n` and `4n`;
/// on APK-381 they are `2n` and `6n`. All the protocol requires is that they be large enough,
/// which is the `DomainSet` contract and what `new` re-asserts.
#[derive(Clone)]
pub struct Domains<F: PrimeField, D: DomainSet<F>> {
    /// The three domains, and the shift derived from their sizes.
    domains: D,

    /// First Lagrange basis polynomial L_0 of degree n evaluated over the large domain; L_0(\omega^0) = 1
    pub l_first_evals_over_4x: Evals<F>,
    /// Last  Lagrange basis polynomial L_{n-1} of degree n evaluated over the large domain; L_{n-1}(\omega^{n-1}}) = 1
    pub l_last_evals_over_4x: Evals<F>,
    /// \omega, a primitive n-th root of unity. Multiplicative generator of the base domain.
    pub omega: F,
    /// \omega^{n-1}
    pub omega_inv: F,
    /// The smaller domain size.
    pub size: usize,
}

impl<F: PrimeField, D: DomainSet<F>> Domains<F, D> {
    /// Builds the domains for a trace of at least `min_domain_size` rows.
    ///
    /// Panics if the field has no domain that large, so it is compiled only for tests and
    /// benchmarks; production code goes through `DomainSet::for_min_size` and [`Self::from_set`].
    #[cfg(any(test, feature = "test-utils"))]
    #[allow(clippy::expect_used, reason = "test-only shorthand")]
    pub fn new(min_domain_size: usize) -> Self {
        Self::from_set(
            D::for_min_size(min_domain_size).expect("no evaluation domain of the requested size"),
        )
    }

    /// Builds the precomputations for an already-chosen triple.
    ///
    /// Used where the triple is already in hand — the keyset owns one, and re-deriving it from a
    /// size would round up a second time.
    pub fn from_set(domains: D) -> Self {
        let n = domains.base().size();
        // The DomainSet contract, restated where it is relied on: everything below assumes the
        // large domain can hold the highest-degree constraint polynomial.
        assert!(
            domains.medium().size() >= 2 * n - 1,
            "domain2x too small: {} < {}",
            domains.medium().size(),
            2 * n - 1
        );
        assert!(
            domains.large().size() >= 4 * n - 3,
            "domain4x too small: {} < {}",
            domains.large().size(),
            4 * n - 3
        );

        let l_first = Self::first_lagrange_basis_polynomial(n);
        let l_last = Self::last_lagrange_basis_polynomial(n);
        let l_first_evals_over_4x = Self::_amplify(l_first, domains.base(), domains.large());
        let l_last_evals_over_4x = Self::_amplify(l_last, domains.base(), domains.large());

        Domains {
            omega: domains.base().generator(),
            omega_inv: domains.base().generator_inv(),
            size: n,
            l_first_evals_over_4x,
            l_last_evals_over_4x,
            domains,
        }
    }

    /// The triple these precomputations were built for.
    pub fn set(&self) -> &D {
        &self.domains
    }

    /// The trace domain `H`, of size `n`.
    pub fn domain(&self) -> &D::Domain {
        self.domains.base()
    }

    /// The domain holding at least `2n - 1` points.
    pub fn domain2x(&self) -> &D::Domain {
        self.domains.medium()
    }

    /// The domain holding at least `4n - 3` points.
    pub fn domain4x(&self) -> &D::Domain {
        self.domains.large()
    }

    /// Interpolates evaluations over the base domain into the unique polynomial of degree `< n`
    /// taking them.
    pub fn interpolate(&self, evals: Vec<F>) -> DensePolynomial<F> {
        // TODO: assert evals.len()
        DensePolynomial::from_coefficients_vec(self.domain().interpolate(&evals))
    }

    /// Interpolates evaluations taken over the 4x domain back into coefficient form.
    /// Exact for polynomials of degree at most `4n - 4`, since `from_set` guarantees
    /// `domain4x.size() >= 4n - 3`; that covers every product of registers the PIOP forms.
    pub fn interpolate_4x(&self, evals: &Evals<F>) -> DensePolynomial<F> {
        DensePolynomial::from_coefficients_vec(self.domain4x().interpolate(&evals.evals))
    }

    /// Interpolates evaluations taken over the 2x domain back into coefficient form.
    pub fn interpolate_2x(&self, evals: &Evals<F>) -> DensePolynomial<F> {
        DensePolynomial::from_coefficients_vec(self.domain2x().interpolate(&evals.evals))
    }

    /// Evaluates a polynomial of degree `< n` over the large domain.
    pub fn amplify_polynomial(&self, poly: &DensePolynomial<F>) -> Evals<F> {
        // TODO: assert poly.degree()
        self.domain4x().fft(&poly.coeffs).into()
    }

    /// Evaluations over the large domain of the polynomial taking `evals` over the base domain.
    pub fn amplify(&self, evals: Vec<F>) -> Evals<F> {
        Self::_amplify(evals, self.domain(), self.domain4x())
    }

    /// Whether the polynomial vanishes at every point of the base domain.
    pub fn is_zero(&self, poly: &DensePolynomial<F>) -> bool {
        self.domain().divide_by_vanishing_poly(poly).1 == DensePolynomial::zero()
    }

    /// Divides by the vanishing polynomial `X^n - 1` of the base domain, returning
    /// `(quotient, remainder)`.
    pub fn compute_quotient(
        &self,
        poly: &DensePolynomial<F>,
    ) -> (DensePolynomial<F>, DensePolynomial<F>) {
        self.domain().divide_by_vanishing_poly(poly)
    }

    /// `c * L_{n-1}`, of degree `n - 1`, evaluated over the large domain.
    pub fn l_last_scaled_by(&self, c: F) -> Evals<F> {
        &self.constant_4x(c) * &self.l_last_evals_over_4x
    }

    /// The constant polynomial `c` evaluated over the large domain.
    pub fn constant_4x(&self, c: F) -> Evals<F> {
        // TODO: ConstantEvaluations to save memory
        vec![c; self.domain4x().size()].into()
    }

    /// Evaluations over `domain_nx` of the degree `< n` polynomial taking `evals` over `domain`.
    /// Products of such polynomials can then be formed pointwise, and are exact as long as
    /// their degree stays below `domain_nx.size()`.
    fn _amplify(evals: Vec<F>, domain: &D::Domain, domain_nx: &D::Domain) -> Evals<F> {
        let coeffs = domain.interpolate(&evals);
        domain_nx.fft(&coeffs).into()
    }

    /// Evaluations of `p(Xw)` over the 4x domain, given `p` and its evaluations there.
    ///
    /// `p(Xw)` represents the left circular shift of the register `p` interpolates, which the
    /// affine-addition constraints need to relate consecutive rows. Whether that is a free
    /// rotation or an extra transform follows from how the two domains' sizes relate, so it is
    /// delegated to the domain set; see [`DomainSet::shift_over_large`].
    pub fn shift_over_4x(&self, poly: &DensePolynomial<F>, evals_over_4x: &Evals<F>) -> Evals<F> {
        self.domains
            .shift_over_large(poly, &evals_over_4x.evals)
            .into()
    }

    fn first_lagrange_basis_polynomial(domain_size: usize) -> Vec<F> {
        Self::li(0, domain_size)
    }

    fn last_lagrange_basis_polynomial(domain_size: usize) -> Vec<F> {
        Self::li(domain_size - 1, domain_size)
    }

    fn li(i: usize, domain_size: usize) -> Vec<F> {
        let mut li = vec![F::zero(); domain_size];
        li[i] = F::one();
        li
    }

    /// For a polynomial p returns a polynomial p' such that p'(H) = p(gH).
    fn coset_polynomial(poly: &DensePolynomial<F>, g: F) -> DensePolynomial<F> {
        let coset_coeffs = poly
            .coeffs
            .iter()
            .scan(F::one(), |pow, &coeff| {
                let coset_coeff = *pow * coeff;
                *pow *= g;
                Some(coset_coeff)
            })
            .collect();
        DensePolynomial::from_coefficients_vec(coset_coeffs)
    }

    /// Evaluations over the 2x domain, given evaluations over the base domain.
    ///
    /// When the domains are nested at index 2 the larger domain is `H ∪ gH`, and since
    /// `p(gH) = p'(H)` for `p'(X) = p(gX)`, the transform can be done as two n-point FFTs whose
    /// results interleave, rather than one 2n-point FFT. Both configurations nest at index 2,
    /// but that is a property of how their sizes are chosen rather than something the trait
    /// promises, so it is checked rather than assumed.
    pub fn amplify_x2(&self, evals: Vec<F>) -> Evals<F> {
        match crate::domain::nesting_index(self.domain(), self.domain2x()) {
            Some(2) => {
                let poly = self.interpolate(evals.clone());
                let coset_poly = Self::coset_polynomial(&poly, self.domain2x().generator());
                let coset_evals = self.domain().fft(&coset_poly.coeffs);
                evals
                    .into_iter()
                    .zip(coset_evals)
                    .flat_map(|(e, ce)| [e, ce])
                    .collect::<Vec<_>>()
                    .into()
            }
            _ => Self::_amplify(evals, self.domain(), self.domain2x()),
        }
    }

    /// Same as [`Self::amplify`].
    pub fn amplify_x4(&self, evals: Vec<F>) -> Evals<F> {
        Self::_amplify(evals, self.domain(), self.domain4x())
    }
}

#[cfg(test)]
mod tests {
    use ark_bw6_761::Fr;
    use ark_poly::{EvaluationDomain, Evaluations, Polynomial, Radix2EvaluationDomain};
    use ark_std::{test_rng, One, UniformRand};

    use crate::domain::Radix2DomainSet;

    use super::*;

    type TestDomains = Domains<Fr, Radix2DomainSet<Fr>>;

    /// Pins, directly against arkworks, the identity the coset fast path in `amplify_x2`
    /// relies on.
    #[test]
    fn test_coset_amplify() {
        let rng = &mut test_rng();
        let n = 64;

        // Let H < G be a subgroup of index 2 (meaning |G| = 2|H|).
        // Then G = H \cup gH, where g is a generator of G.

        // Let |H| = n. Observe that
        // evaluations of a polynomial p(X) = a_0 + ... + a_{n-1}.X^{n-1} over a coset gH
        // are equal to the
        // evaluations of the polynomial p'(X) = a_0 + ... + (a_{n-1}.g^{n-1}).X^{n-1} over the
        // subgroup H: p(gH) = p'(H).

        // Thus p(G) can be computed either with a 2n-FFT,
        // or as p(G) = p(H \cup gH) = p(H) \cup p(gH) = p(H) \cup p'(H) with 2 n-FFTs.
        // In the case when p(H) is already known, the latter approach might be more efficient.

        let domain = Radix2EvaluationDomain::<Fr>::new(n).unwrap(); // H
        let domain2x = Radix2EvaluationDomain::<Fr>::new(2 * n).unwrap(); // G
        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>(); // p(H)
        let poly = Evaluations::from_vec_and_domain(evals.clone(), domain).interpolate(); // p
        let evals2x = poly.evaluate_over_domain_by_ref(domain2x); // p(G)

        let root2x = domain2x.group_gen; // g
        let coset_coeffs = poly
            .coeffs
            .iter()
            .scan(Fr::one(), |pow, &coeff| {
                let coset_coeff = *pow * coeff;
                *pow = *pow * root2x;
                Some(coset_coeff)
            })
            .collect();

        let coset_poly = DensePolynomial::from_coefficients_vec(coset_coeffs); // p'
        let coset_evals = coset_poly.evaluate_over_domain_by_ref(domain); // p'(H)

        let evals2x_2: Vec<_> = evals
            .into_iter()
            .zip(coset_evals.evals)
            .flat_map(|(e, ce)| vec![e, ce])
            .collect(); // p(G)

        assert_eq!(evals2x.evals, evals2x_2);
    }

    #[test]
    fn test_amplify() {
        let rng = &mut test_rng();
        let n = 64;

        let domains = TestDomains::new(n);

        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let poly = domains.interpolate(evals.clone());

        let evals4x_from_poly = domains.amplify_polynomial(&poly);
        let evals4x_from_vec = domains.amplify(evals);

        assert_eq!(evals4x_from_poly, evals4x_from_vec);
        assert_eq!(domains.interpolate_4x(&evals4x_from_poly), poly);
    }

    /// `amplify_x2` must produce the same evaluations as a direct 2n-point FFT, whichever path
    /// it takes.
    #[test]
    fn test_amplify_2x() {
        let rng = &mut test_rng();
        let n = 64;

        let domains = TestDomains::new(n);

        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let poly = domains.interpolate(evals.clone());
        let evals2x: Evals<Fr> = domains.domain2x().fft(&poly.coeffs).into();

        assert_eq!(evals2x, domains.amplify_x2(evals));
    }

    #[test]
    fn test_domains_l_last_scaled_by() {
        let rng = &mut test_rng();
        let n = 64;

        let c = Fr::rand(rng);

        let mut c_ln = vec![Fr::zero(); n];
        c_ln[n - 1] = c;

        let domains = TestDomains::new(n);

        assert_eq!(domains.l_last_scaled_by(c), domains.amplify(c_ln));
    }

    /// Radix-2 domains nest at index 4 (large) and 2 (medium).
    #[test]
    fn radix2_domains_are_nested() {
        let domains = TestDomains::new(64);
        assert_eq!(
            crate::domain::nesting_index(domains.domain(), domains.domain4x()),
            Some(4)
        );
        assert_eq!(
            crate::domain::nesting_index(domains.domain(), domains.domain2x()),
            Some(2)
        );
    }

    /// Ground truth for the shift: the result must be `p(Xw)` evaluated over the 4x domain,
    /// point by point, however it was computed.
    #[test]
    fn shift_over_4x_evaluates_the_shifted_polynomial() {
        let rng = &mut test_rng();
        let n = 64;
        let domains = TestDomains::new(n);

        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let poly = domains.interpolate(evals);
        let evals_4x = domains.amplify_polynomial(&poly);

        let shifted = domains.shift_over_4x(&poly, &evals_4x);

        for i in 0..domains.domain4x().size() {
            let point = domains.domain4x().element(i) * domains.omega;
            assert_eq!(
                shifted.evals[i],
                poly.evaluate(&point),
                "mismatch at index {}",
                i
            );
        }
    }

    /// The rotation fast path and the coefficient-scaling fallback must agree. Both production
    /// triples nest, so every proof takes the rotation branch; this pins the fallback, which
    /// only a non-nesting triple would take, against it.
    #[test]
    fn shift_fallback_agrees_with_the_rotation_fast_path() {
        let rng = &mut test_rng();
        let n = 64;
        let domains = TestDomains::new(n);
        assert!(
            crate::domain::nesting_index(domains.domain(), domains.domain4x()).is_some(),
            "this test is only meaningful when the fast path is taken"
        );

        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let poly = domains.interpolate(evals);
        let evals_4x = domains.amplify_polynomial(&poly);

        let via_rotation = domains.shift_over_4x(&poly, &evals_4x);

        // The fallback, spelled out: p(Xw) = sum_i (c_i w^i) X^i.
        let mut coeffs = poly.coeffs.clone();
        let mut power = Fr::one();
        for c in coeffs.iter_mut() {
            *c *= power;
            power *= domains.omega;
        }
        let via_scaling: Evals<Fr> = domains.domain4x().fft(&coeffs).into();

        assert_eq!(via_rotation, via_scaling);
    }

    /// Likewise for the coset fast path in `amplify_x2`.
    #[test]
    fn amplify_x2_coset_path_agrees_with_the_general_path() {
        let rng = &mut test_rng();
        let n = 64;
        let domains = TestDomains::new(n);

        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
        let via_coset = domains.amplify_x2(evals.clone());

        let poly = domains.interpolate(evals);
        let via_general: Evals<Fr> = domains.domain2x().fft(&poly.coeffs).into();

        assert_eq!(via_coset, via_general);
    }

    /// The expanded domains must be sized off the realized base size, not the requested one.
    #[test]
    fn test_expanded_domains_are_sized_off_realized_base() {
        let domains = TestDomains::new(33); // rounds up to 64
        assert_eq!(domains.size, 64);
        assert!(domains.domain2x().size() >= 2 * 64 - 1);
        assert!(domains.domain4x().size() >= 4 * 64 - 3);
    }
}
