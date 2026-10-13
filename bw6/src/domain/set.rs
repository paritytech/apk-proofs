//! The triple of evaluation domains the PIOP works over, and the two ways of producing it.
//!
//! The protocol needs three domains: the trace domain `H` of size `n`, and two larger ones that
//! hold products of register polynomials in evaluation form.
//!
//! The protocol's demand on the two larger domains is only capacity — `2n - 1` and `4n - 3`
//! points — and that is all this module states. Which sizes a field actually offers, and how it
//! reaches them from `n`, is left to the [`DomainSizes`] table the configuration supplies.
//! Everything downstream — [`crate::Keyset`], [`crate::domains::Domains`], the PIOP, the prover
//! and the verifier — carries a single `D: DomainSet` parameter and never asks which curve it
//! is on.

use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use core::marker::PhantomData;

use super::cooley_tukey::CooleyTukeyDomain;
use super::radix2::Radix2Domain;
use super::subgroup::subgroup_generator;
use super::types::{DomainError, FftDomain, SupportsPackedScheme};

/// The three domain sizes one proof is computed over.
///
/// `base` is the trace domain `H`. The other two hold products of register polynomials in
/// evaluation form, and the protocol's only demand on them is capacity: the constraint
/// polynomials reach degree `2n - 2` and `4n - 4`, so recovering them takes `2n - 1` and
/// `4n - 3` evaluations. Anything at or above those is correct. Which sizes a field actually
/// offers, and which of them are worth choosing, is not the protocol's business.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DomainTriple {
    /// The trace domain `H`, of size `n`.
    pub base: usize,
    /// Holds at least `2n - 1` points.
    pub medium: usize,
    /// Holds at least `4n - 3` points.
    pub large: usize,
}

impl DomainTriple {
    pub const fn new(base: usize, medium: usize, large: usize) -> Self {
        DomainTriple {
            base,
            medium,
            large,
        }
    }

    /// The protocol invariant, in one place. Both bounds are floors, not targets: a domain set
    /// is free to overshoot them, and on a field whose subgroup orders are sparse it must.
    pub const fn meets_protocol_bounds(&self) -> bool {
        self.base >= 1 && self.medium >= 2 * self.base - 1 && self.large >= 4 * self.base - 3
    }
}

/// The three domains a proof is computed over, chosen together.
///
/// Chosen *together* because they are not independent: the protocol requires
/// `medium >= 2n - 1` and `large >= 4n - 3`, and whether a domain of a given size exists at all
/// is a property of the field. Picking `n` first and then asking for `4n` separately is how the
/// sizes silently fail to line up.
pub trait DomainSet<F: PrimeField>: Clone + Sized {
    /// The underlying single-domain implementation. Radix-2 for APK-377, mixed-radix
    /// Cooley-Tukey for APK-381.
    type Domain: FftDomain<F>;

    /// The cheapest usable triple whose base domain holds at least `min_size` points.
    ///
    /// The base size is rounded **up** to something the field supports, so callers must read
    /// back `base().size()` rather than assuming they got what they asked for.
    fn for_min_size(min_size: usize) -> Result<Self, DomainError>;

    /// The trace domain `H`, of size `n`.
    fn base(&self) -> &Self::Domain;

    /// A domain of at least `2n - 1` points.
    fn medium(&self) -> &Self::Domain;

    /// A domain of at least `4n - 3` points.
    fn large(&self) -> &Self::Domain;

    /// The triple whose base domain has *exactly* `size` points.
    ///
    /// The verifier needs the exactness: the domain size arrives inside a keyset commitment,
    /// which for a bridge comes from the chain and is untrusted. Rounding it up there would let
    /// a proof verify against a domain the prover did not use.
    fn for_exact_size(size: usize) -> Result<Self, DomainError> {
        let set = Self::for_min_size(size)?;
        let realised = set.base().size();
        if realised != size {
            return Err(DomainError::NotExact {
                requested: size,
                nearest: realised,
            });
        }
        Ok(set)
    }

    /// Just the base domain of [`for_exact_size`](Self::for_exact_size), without building the
    /// two larger ones.
    ///
    /// The verifier is the only caller and needs nothing else, and the large domain is where the
    /// cost is: at 100,000 validators its transform plan holds on the order of a million
    /// precomputed twiddle factors that the verifier would never touch.
    ///
    /// Must accept exactly the sizes `for_exact_size` accepts, and reject the rest with the same
    /// error — otherwise a proof could verify against a domain no prover would have chosen.
    fn base_for_exact_size(size: usize) -> Result<Self::Domain, DomainError> {
        Ok(Self::for_exact_size(size)?.base().clone())
    }

    /// Evaluations over `large` of `p(Xw)`, where `w` generates `base` — the left circular shift
    /// of the register `p` interpolates, which the affine-addition constraints need in order to
    /// relate consecutive rows.
    ///
    /// Derived from the sizes rather than configured: a rotation when the domains nest, a
    /// coefficient scaling and one extra transform when they do not. No implementation here
    /// overrides it, and one that did would be claiming to know something about its domains
    /// that their sizes do not already say.
    fn shift_over_large(&self, poly: &DensePolynomial<F>, evals_over_large: &[F]) -> Vec<F> {
        shifted_evals(self.base(), self.large(), poly, evals_over_large)
    }
}

/// Evaluations over `large` of `p(Xw)`, where `w` generates `small`.
///
/// When the domains nest — `|large| = k * |small|` with `w_L^k = w` — this is a rotation of the
/// evaluation vector by `k`, and costs nothing: `p(w_L^{i+k}) = p(w_L^i * w)`. Otherwise `p(Xw)`
/// is formed in coefficient form, where it is `sum_i (c_i w^i) X^i`, and transformed. That is
/// always correct but costs one extra FFT over the largest domain the protocol uses.
///
/// Which branch applies follows from the two sizes, so no domain set has to declare it. The
/// nesting relation itself is `debug_assert`ed rather than re-derived in release: both domain
/// families here generate their subgroups as `GENERATOR^((q-1)/n)`, under which divisibility of
/// the sizes implies it, and checking would mean a field exponentiation on a path the prover
/// takes once per register.
pub(crate) fn shifted_evals<F: PrimeField, D: FftDomain<F>>(
    small: &D,
    large: &D,
    poly: &DensePolynomial<F>,
    evals_over_large: &[F],
) -> Vec<F> {
    let (n, m) = (small.size(), large.size());
    if n > 0 && m % n == 0 {
        let k = m / n;
        debug_assert_eq!(
            nesting_index(small, large),
            Some(k),
            "sizes {} and {} divide but the domains do not nest",
            n,
            m
        );
        let mut shifted = evals_over_large.to_vec();
        shifted.rotate_left(k);
        return shifted;
    }

    let mut coeffs = poly.coeffs.clone();
    let mut power = F::one();
    let w = small.generator();
    for c in coeffs.iter_mut() {
        *c *= power;
        power *= w;
    }
    large.fft(&coeffs)
}

/// `k` such that `large.element(k) == small.generator()`, when `small` sits inside `large` with
/// its elements interleaved in the natural way; `None` otherwise.
///
/// Both domain families here generate their subgroups as `GENERATOR^((q-1)/n)`, so a divisibility
/// of sizes does imply nesting — but that is a property of the construction, not a theorem, so
/// it is checked rather than assumed.
pub(crate) fn nesting_index<F: PrimeField, D: FftDomain<F>>(small: &D, large: &D) -> Option<usize> {
    let (n, m) = (small.size(), large.size());
    if n == 0 || m % n != 0 {
        return None;
    }
    let k = m / n;
    (large.element(k) == small.generator()).then_some(k)
}

// -------------------------------------------------------------------------------------------
// APK-377
// -------------------------------------------------------------------------------------------

/// Power-of-two domains `n, 2n, 4n`. Used by APK-377.
///
/// BW6-761's scalar field has two-adicity 46, so every size the prover could afford is available
/// and the classical radix-2 layout applies unchanged.
#[derive(Clone, Copy, Debug)]
pub struct Radix2DomainSet<F: PrimeField> {
    base: Radix2Domain<F>,
    medium: Radix2Domain<F>,
    large: Radix2Domain<F>,
}

impl<F: PrimeField> DomainSet<F> for Radix2DomainSet<F> {
    type Domain = Radix2Domain<F>;

    fn for_min_size(min_size: usize) -> Result<Self, DomainError> {
        let base = Radix2Domain::try_new(min_size)?;
        let n = base.size();
        Ok(Radix2DomainSet {
            base,
            medium: Radix2Domain::try_new(2 * n)?,
            large: Radix2Domain::try_new(4 * n)?,
        })
    }

    fn base(&self) -> &Radix2Domain<F> {
        &self.base
    }

    fn medium(&self) -> &Radix2Domain<F> {
        &self.medium
    }

    fn large(&self) -> &Radix2Domain<F> {
        &self.large
    }

    fn base_for_exact_size(size: usize) -> Result<Radix2Domain<F>, DomainError> {
        let base = Radix2Domain::<F>::try_new(size)?;
        if base.size() != size {
            return Err(DomainError::NotExact {
                requested: size,
                nearest: base.size(),
            });
        }
        Ok(base)
    }
}

/// Radix-2 sizes from 256 up are all multiples of 256, which is what `packed` needs.
impl<F: PrimeField> SupportsPackedScheme for Radix2DomainSet<F> {}

// -------------------------------------------------------------------------------------------
// Mixed-radix: fields without the two-adicity for radix-2
// -------------------------------------------------------------------------------------------

/// A field's precomputed answer to "which domains should a trace of `min_size` rows use?".
///
/// The implementor owns both the storage and the reading: a list of base sizes expanded by some
/// rule, a list of `(base, medium, large)` tuples read straight off, whatever suits the field's
/// multiplicative order. Working it out once and writing it down is the point — deriving it at
/// run time would mean a modular exponentiation per candidate, paid on every keyset, prover and
/// verifier, to recompute a constant.
///
/// All this trait asks of the answer is [`DomainTriple::meets_protocol_bounds`] and that all
/// three sizes be subgroup orders of `F`. Neither is checkable by the type system, so both are
/// asserted where the table is written down and `debug_assert`ed wherever a set is built.
///
/// `F` is a parameter so a table worked out for one field cannot be attached to another.
pub trait DomainSizes<F: PrimeField>: 'static {
    /// The triple for a trace of at least `min_size` rows, or `None` when the field has nothing
    /// that large. `base` is rounded **up**, so callers must read back what they got.
    ///
    /// Must be deterministic: the prover and the verifier resolve the same request
    /// independently and have to land on the same domains.
    fn triple_for(min_size: usize) -> Option<DomainTriple>;
}

/// Mixed-radix domains, with the sizes supplied by `S`. Used by APK-381.
pub struct SmoothDomainSet<F: PrimeField, S: DomainSizes<F>> {
    base: CooleyTukeyDomain<F>,
    medium: CooleyTukeyDomain<F>,
    large: CooleyTukeyDomain<F>,
    _sizes: PhantomData<fn() -> S>,
}

// Hand-written rather than derived: `S` is a marker that is never held by value, so deriving
// would saddle every user of the type with spurious `S: Clone` / `S: Debug` bounds.
impl<F: PrimeField, S: DomainSizes<F>> Clone for SmoothDomainSet<F, S> {
    fn clone(&self) -> Self {
        SmoothDomainSet {
            base: self.base.clone(),
            medium: self.medium.clone(),
            large: self.large.clone(),
            _sizes: PhantomData,
        }
    }
}

impl<F: PrimeField, S: DomainSizes<F>> core::fmt::Debug for SmoothDomainSet<F, S> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("SmoothDomainSet")
            .field("base", &self.base.size())
            .field("medium", &self.medium.size())
            .field("large", &self.large.size())
            .finish()
    }
}

impl<F: PrimeField, S: DomainSizes<F>> SmoothDomainSet<F, S> {
    /// Builds the three transforms for sizes `S` has already chosen.
    ///
    /// The two `debug_assert`s are the `DomainSizes` contract, checked where it is relied on but
    /// kept off the hot path: a table paired with the wrong field, or one entry mistyped, names
    /// itself here in every debug and test build rather than surfacing as a confusing
    /// `TooLarge` at run time.
    fn build(sizes: DomainTriple) -> Result<Self, DomainError> {
        debug_assert!(
            sizes.meets_protocol_bounds(),
            "{:?} is too small to hold the constraint polynomials",
            sizes
        );
        debug_assert!(
            [sizes.base, sizes.medium, sizes.large]
                .iter()
                .all(|&n| subgroup_generator::<F>(n).is_some()),
            "{:?} names a size that is not a subgroup order of this field",
            sizes
        );

        let missing = |n| move || DomainError::TooLarge { requested: n };
        Ok(SmoothDomainSet {
            base: CooleyTukeyDomain::new(sizes.base).ok_or_else(missing(sizes.base))?,
            medium: CooleyTukeyDomain::new(sizes.medium).ok_or_else(missing(sizes.medium))?,
            large: CooleyTukeyDomain::new(sizes.large).ok_or_else(missing(sizes.large))?,
            _sizes: PhantomData,
        })
    }
}

impl<F: PrimeField, S: DomainSizes<F>> DomainSet<F> for SmoothDomainSet<F, S> {
    type Domain = CooleyTukeyDomain<F>;

    /// Whatever `S` answers. Which sizes exist, and how the two larger ones are derived from
    /// the base, is entirely the table's business.
    fn for_min_size(min_size: usize) -> Result<Self, DomainError> {
        Self::build(S::triple_for(min_size).ok_or(DomainError::TooLarge {
            requested: min_size,
        })?)
    }

    fn base_for_exact_size(size: usize) -> Result<CooleyTukeyDomain<F>, DomainError> {
        let chosen = S::triple_for(size).ok_or(DomainError::TooLarge { requested: size })?;
        if chosen.base != size {
            return Err(DomainError::NotExact {
                requested: size,
                nearest: chosen.base,
            });
        }
        CooleyTukeyDomain::new(size).ok_or(DomainError::TooLarge { requested: size })
    }

    fn base(&self) -> &CooleyTukeyDomain<F> {
        &self.base
    }

    fn medium(&self) -> &CooleyTukeyDomain<F> {
        &self.medium
    }

    fn large(&self) -> &CooleyTukeyDomain<F> {
        &self.large
    }
}

// Deliberately **not** `impl SupportsPackedScheme for SmoothDomainSet`: `packed` splits the
// bitmask into 256-bit chunks and so needs `256 | n`, which a table of arbitrary subgroup
// orders cannot promise. On BW6-767, the field this exists for, it is outright impossible —
// `q - 1` carries a single factor of 2, so no domain size there is even a multiple of 4.

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bw6_761::Fr as Fr761;
    use ark_bw6_767::Fr as Fr767;
    use ark_poly::DenseUVPolynomial;
    use ark_std::test_rng;

    // The mechanism is generic, but exercising it needs a concrete size list, so these tests
    // borrow APK-381's. Whether that list is *correct* for BW6-767 is asserted next to the list
    // itself, in `instances::bls12_381_bw6_767`.
    use crate::instances::bls12_381_bw6_767::{Apk381DomainSizes, APK381_DOMAIN_SIZES};

    type Smooth = crate::instances::bls12_381_bw6_767::Domains767;
    type Radix2 = Radix2DomainSet<Fr761>;

    /// What the PIOP requires of a triple, whatever the table chose. The expansion rule is not
    /// asserted here — by what factor a table overshoots the floors is its own business, and
    /// APK-381's answer is pinned in its own module.
    #[test]
    fn every_triple_meets_the_protocol_bounds_and_nests() {
        // Building all 30 means transform plans for domains in the hundreds of millions. The
        // requirements are structural, so the small end proves them.
        for &n in APK381_DOMAIN_SIZES.iter().take(12) {
            let sizes = <Apk381DomainSizes as DomainSizes<Fr767>>::triple_for(n)
                .expect("an entry must resolve to itself");
            assert_eq!(sizes.base, n);
            assert!(sizes.meets_protocol_bounds(), "{:?}", sizes);

            let set = Smooth::build(sizes).expect("entry must be constructible");
            assert_eq!(set.base().size(), n);
            assert!(set.medium().size() >= 2 * n - 1);
            assert!(set.large().size() >= 4 * n - 3);

            // Nested, which is what makes `shifted_evals` take the rotation branch. At what
            // index is the table's choice, not the protocol's.
            assert!(nesting_index(set.base(), set.medium()).is_some());
            assert!(nesting_index(set.base(), set.large()).is_some());
        }
    }

    /// The floors are floors: a triple below either one is rejected, at or above is accepted.
    #[test]
    fn protocol_bounds_are_the_degree_bounds_of_the_constraint_polynomials() {
        let n = 16;
        assert!(DomainTriple::new(n, 2 * n - 1, 4 * n - 3).meets_protocol_bounds());
        assert!(!DomainTriple::new(n, 2 * n - 2, 4 * n - 3).meets_protocol_bounds());
        assert!(!DomainTriple::new(n, 2 * n - 1, 4 * n - 4).meets_protocol_bounds());
        // Overshooting is fine, and on a sparse field unavoidable.
        assert!(DomainTriple::new(n, 2 * n, 6 * n).meets_protocol_bounds());
    }

    /// The selector returns the first entry at least as large as requested...
    #[test]
    fn for_min_size_rounds_up_to_a_table_entry() {
        for (requested, expected) in [(1, 1), (2, 3), (12, 23), (34, 47), (48, 69), (254, 517)] {
            let set = Smooth::for_min_size(requested).unwrap();
            assert_eq!(set.base().size(), expected, "for_min_size({})", requested);
            assert!(APK381_DOMAIN_SIZES.contains(&set.base().size()));
        }
    }

    /// ...and because the dominated sizes are simply not in the table, that plain binary search
    /// never returns one. 5000 validators would otherwise land on 10177, which is prime and has
    /// to go through Rader.
    #[test]
    fn a_request_is_never_served_by_a_dominated_size() {
        assert_eq!(Smooth::for_min_size(5000).unwrap().base().size(), 11891);
        assert_eq!(Smooth::for_min_size(20000).unwrap().base().size(), 35673);
    }

    /// Sizes of interest land where the table says they do.
    #[test]
    fn validator_counts_of_interest_land_on_the_documented_entries() {
        // A keyset of k needs k + 1 points: the extra slot is the accumulator's seed.
        assert_eq!(Smooth::for_min_size(1000 + 1).unwrap().base().size(), 1081);
        assert_eq!(Smooth::for_min_size(1500 + 1).unwrap().base().size(), 1551);
    }

    #[test]
    fn asking_for_more_than_the_field_holds_is_an_error() {
        let too_big = APK381_DOMAIN_SIZES.last().unwrap() + 1;
        assert_eq!(
            Smooth::for_min_size(too_big).unwrap_err(),
            DomainError::TooLarge { requested: too_big }
        );
    }

    /// The verifier's constructor: a size off the table must not be silently rounded.
    #[test]
    fn for_exact_size_rejects_a_size_not_in_the_table() {
        assert_eq!(Smooth::for_exact_size(253).unwrap().base().size(), 253);
        assert_eq!(
            Smooth::for_exact_size(254).unwrap_err(),
            DomainError::NotExact {
                requested: 254,
                nearest: 517
            }
        );
    }

    /// The verifier's cheap path must accept and reject exactly what the full one does. If it
    /// drifted, a proof could verify against a domain no prover would have chosen.
    #[test]
    fn base_for_exact_size_agrees_with_for_exact_size() {
        fn check<F: ark_ff::PrimeField, D: DomainSet<F>>(sizes: &[usize]) {
            for &size in sizes {
                match D::for_exact_size(size) {
                    Ok(set) => assert_eq!(
                        D::base_for_exact_size(size).unwrap().size(),
                        set.base().size(),
                        "size {}",
                        size
                    ),
                    Err(e) => assert_eq!(
                        D::base_for_exact_size(size).map(|d| d.size()).err(),
                        Some(e),
                        "size {}",
                        size
                    ),
                }
            }
        }

        // Table entries, non-entries, the entry the cost rule makes unreachable (10177), and
        // sizes past the end of the field.
        check::<Fr767, Smooth>(&[1, 3, 11, 23, 33, 46, 47, 253, 254, 517, 759, 1081, 10177]);
        check::<Fr761, Radix2>(&[1, 2, 16, 17, 256, 300, 1024]);
    }

    fn random_poly<F: ark_ff::PrimeField>(degree_bound: usize) -> DensePolynomial<F> {
        let rng = &mut test_rng();
        DensePolynomial::from_coefficients_vec((0..degree_bound).map(|_| F::rand(rng)).collect())
    }

    /// The rotation branch of `shifted_evals` must agree with the coefficient-scaling one, and
    /// both with direct evaluation of `p(Xw)`, on both configurations' domains.
    #[test]
    fn rotation_agrees_with_coefficient_scaling() {
        fn check<F: ark_ff::PrimeField, D: DomainSet<F>>(set: &D) {
            let poly = random_poly::<F>(set.base().size());
            let evals = set.large().fft(&poly.coeffs);

            let rotated = set.shift_over_large(&poly, &evals);

            // The coefficient-scaling branch, spelled out: p(Xw) has coefficients c_i * w^i.
            let mut coeffs = poly.coeffs.clone();
            let mut power = F::one();
            for c in coeffs.iter_mut() {
                *c *= power;
                power *= set.base().generator();
            }
            let scaled = set.large().fft(&coeffs);

            assert_eq!(rotated, scaled);

            // And both are really the evaluations of p(Xw) over the large domain.
            let w = set.base().generator();
            for i in 0..set.large().size() {
                let x = set.large().element(i);
                assert_eq!(rotated[i], evaluate(&poly, x * w));
            }
        }

        fn evaluate<F: ark_ff::PrimeField>(p: &DensePolynomial<F>, x: F) -> F {
            p.coeffs.iter().rev().fold(F::zero(), |acc, c| acc * x + c)
        }

        check(&Radix2::for_min_size(16).unwrap());
        check(&Smooth::for_min_size(23).unwrap());
        check(&Smooth::for_min_size(33).unwrap());
    }

    /// A size list paired with a field that cannot realise it is a programming error, and the
    /// `debug_assert` in `build` is what turns it into a loud one. Without that net, dropping
    /// the old runtime filter would have made such a mistake surface as a puzzling `TooLarge`.
    ///
    /// Debug-only, because that is exactly the point: the check costs a modular exponentiation
    /// and must not run in a release prover or verifier.
    #[test]
    #[cfg(debug_assertions)]
    #[should_panic(expected = "not a subgroup order of this field")]
    fn a_size_list_that_the_field_cannot_realise_trips_the_debug_assert() {
        struct WrongForThisField;
        impl DomainSizes<Fr767> for WrongForThisField {
            fn triple_for(_min_size: usize) -> Option<DomainTriple> {
                // Powers of two. BW6-767's scalar field has two-adicity 1, so no such subgroup
                // exists there; this table belongs to a radix-2 field.
                Some(DomainTriple::new(256, 512, 1024))
            }
        }

        let _ = SmoothDomainSet::<Fr767, WrongForThisField>::for_min_size(200);
    }

    /// APK-377 is unchanged: powers of two, and the classical `n, 2n, 4n` layout.
    #[test]
    fn radix2_set_is_the_classical_layout() {
        let set = Radix2::for_min_size(200).unwrap();
        assert_eq!(set.base().size(), 256);
        assert_eq!(set.medium().size(), 512);
        assert_eq!(set.large().size(), 1024);
        assert_eq!(nesting_index(set.base(), set.large()), Some(4));
    }
}
