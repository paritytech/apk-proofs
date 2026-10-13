//! APK-381: BLS12-381 signatures, proofs over BW6-767.
//!
//! This module provides type aliases and constants for APK proofs using:
//! - **Inner curve**: BLS12-381 G1, where the BLS public keys live.
//! - **Outer curve**: BW6-767 G1, where the proof's commitments live. Both the prover and the
//!   verifier do arithmetic here.
//!
//! The two form a 2-chain: BW6-767's scalar field is BLS12-381's base field, so the inner
//! curve's coordinates are native field elements of the proof system. BW6-767 is the BW6 curve
//! over BLS12-381 from El Housni and Guillevic, "Families of SNARK-friendly 2-chains of
//! elliptic curves", <https://eprint.iacr.org/2021/1359>; see also
//! <https://hackmd.io/@gnark/bw6_bls12381>.
//!
//! Unlike BW6-761, its scalar field has two-adicity 1, which is what the mixed-radix domain
//! machinery in [`crate::domain`] exists for.
//!
//! ## Polynomial Commitment Schemes
//!
//! - [`kzg`]: KZG commitments on BW6-767, the only scheme implemented.

use ark_bls12_381::G1Projective as Bls12_381_G1;
use ark_bw6_767::{Fr, G1Affine as BW6_767_G1Affine, G1Projective as BW6_767_G1};

use crate::{AccountablePublicInput, CountingPublicInput, Keyset};

// ============================================================================
// Polynomial Commitment Schemes
// ============================================================================

/// KZG polynomial commitment scheme types for this pairing
pub mod kzg;

// ============================================================================
// Curve Type Aliases
// ============================================================================

/// Inner curve: BLS12-381 G1
///
/// Used for:
/// - BLS signature public keys
/// - BLS signature aggregate public keys
/// - Elements being committed to in the keyset
pub type InnerCurve = Bls12_381_G1;

/// The pairing the inner curve belongs to. Needed to sign with these keys; see [`crate::bls`].
pub type InnerPairing = ark_bls12_381::Bls12_381;

/// Outer curve: BW6-767 G1, in projective form.
///
/// The group the commitments and opening proofs live in, used by both sides of the protocol:
/// - the prover commits to its polynomials and computes the opening proofs here;
/// - the verifier checks every proof point for membership in G1 (see
///   [`crate::PrimeSubgroup`]), rebuilds the commitment to the linearization polynomial and
///   aggregates the commitments opened at `zeta` as linear combinations here, and hands the
///   result to the KZG pairing check on BW6-767.
///
/// The projective form is what the group arithmetic uses; points are stored and serialized in
/// [`OuterAffine`] form.
pub type OuterCurve = BW6_767_G1;

/// Outer curve: BW6-767 G1, in affine form.
///
/// The representation of every outer-curve point that is stored or serialized: keyset and
/// register commitments, the quotient commitment, and KZG opening proofs.
pub type OuterAffine = BW6_767_G1Affine;

/// Outer curve scalar field: BW6-767 Fr.
///
/// Equal to BLS12-381's base field `Fq`, so BLS12-381 point coordinates are elements of it.
/// All polynomial arithmetic in the proof system happens over this field.
pub type OuterScalar = Fr;

// ============================================================================
// PCS-Independent Type Aliases
// ============================================================================

/// The base domain sizes available to APK-381, ascending.
///
/// BW6-767's scalar field has
/// `q - 1 = 2 * 3^2 * 11 * 23 * 47 * 10177 * 859267 * 52437899 * (a 305-bit remainder)`. The
/// two primes after 10177 are out of reach: Rader's algorithm would need a convolution domain of
/// at least `2p - 3` points with prime factors at most 64, and the largest such subgroup order
/// here is `2 * 3^2 * 11 * 23 * 47 = 214038`. So the usable part stops at 10177.
/// Two of those factors are **reserved** rather than spent on the base domain: one 2 and one 3,
/// so that `6n` divides `q - 1` whenever `n` does. That fixes the triple as `n, 2n, 6n` —
/// `2n >= 2n - 1` and `6n >= 4n - 3`, each a multiple of `n` so the domains nest — and leaves
/// the base sizes as the divisors of
///
/// ```text
/// 3 * 11 * 23 * 47 * 10177 = 363,044,121
/// ```
///
/// listed here from the smallest to the largest the field admits. The largest entry's `6n` is
/// `2 * 3^2 * 11 * 23 * 47 * 10177`, the entire usable smooth part of `q - 1`, so the table
/// cannot be extended.
///
/// The price of reserving those factors is padding: the gaps are wide, and the worst case in the
/// range that matters is `n = 3244` rounding up to 11891, a factor of 3.67. Spending the 2 and
/// the 3 on the base domain would close the gaps, but then `2n` and `6n` would not exist and the
/// three domains would have to be chosen independently and would not nest.
///
/// **Two divisors are deliberately absent: 10177 and 30531 = 3 * 10177.** Both are dominated —
/// the next entry up is larger *and* cheaper to transform, because 10177 is prime and has to go
/// through Rader's algorithm rather than a naive small-radix DFT.
///
/// Entries from 111947 up still carry the factor 10177 and still need Rader. There is no smooth
/// alternative that high, so nothing dominates them and they stay.
///
/// The structural claims above are asserted in this module's tests: that each entry is a
/// divisor, that `2n` and `6n` really are subgroup orders of this field, and that the two
/// omissions are the only divisors missing. The cost figures are not.
pub const APK381_DOMAIN_SIZES: &[usize] = &[
    1,         // 1
    3,         // 3
    11,        // 11
    23,        // 23
    33,        // 3 * 11
    47,        // 47
    69,        // 3 * 23
    141,       // 3 * 47
    253,       // 11 * 23
    517,       // 11 * 47
    759,       // 3 * 11 * 23
    1081,      // 23 * 47
    1551,      // 3 * 11 * 47
    3243,      // 3 * 23 * 47
    11891,     // 11 * 23 * 47
    35673,     // 3 * 11 * 23 * 47
    111947,    // 11 * 10177
    234071,    // 23 * 10177
    335841,    // 3 * 11 * 10177
    478319,    // 47 * 10177
    702213,    // 3 * 23 * 10177
    1434957,   // 3 * 47 * 10177
    2574781,   // 11 * 23 * 10177
    5261509,   // 11 * 47 * 10177
    7724343,   // 3 * 11 * 23 * 10177
    11001337,  // 23 * 47 * 10177
    15784527,  // 3 * 11 * 47 * 10177
    33004011,  // 3 * 23 * 47 * 10177
    121014707, // 11 * 23 * 47 * 10177
    363044121, // 3 * 11 * 23 * 47 * 10177
];

/// Binds [`APK381_DOMAIN_SIZES`] to BW6-767's scalar field, and expands a base size to the
/// triple the protocol needs.
pub struct Apk381DomainSizes;

impl crate::DomainSizes<OuterScalar> for Apk381DomainSizes {
    /// Reads the base size off [`APK381_DOMAIN_SIZES`] and expands it to `n, 2n, 6n`.
    ///
    /// `2n` and `6n` rather than the protocol's bare floors of `2n - 1` and `4n - 3` because
    /// these are the sizes this field actually has: `4n` never is, since `4 ∤ q - 1`. One
    /// factor of 2 and one of 3 are held back out of `q - 1` when the table is built, precisely
    /// so that `2n` and `6n` remain subgroup orders for every entry. They clear the floors —
    /// `2n >= 2n - 1` and `6n >= 4n - 3` — and, being whole multiples of `n`, they nest, which
    /// is what lets the shifted register be a rotation of the evaluation vector rather than an
    /// extra transform over the largest domain in the protocol.
    ///
    /// A field with a different factorisation would answer differently, and nothing above this
    /// function knows or cares which numbers come back.
    fn triple_for(min_size: usize) -> Option<crate::DomainTriple> {
        let i = APK381_DOMAIN_SIZES.partition_point(|&n| n < min_size);
        let n = *APK381_DOMAIN_SIZES.get(i)?;
        Some(crate::DomainTriple::new(n, 2 * n, 6 * n))
    }
}

/// Evaluation domains for BW6-767's scalar field.
///
/// Two-adicity is 1 there, so radix-2 does not apply: no power-of-two domain beyond size 2, and
/// no domain size divisible by 4 at all. Sizes come from [`APK381_DOMAIN_SIZES`] and are
/// transformed by mixed-radix Cooley-Tukey, with Rader for the factor 10177.
pub type Domains767 = crate::SmoothDomainSet<OuterScalar, Apk381DomainSizes>;

/// A single mixed-radix domain. Prefer [`Domains767`]: the protocol works over the triple, and
/// the triple is what makes the sizes `n, 2n, 6n` line up.
pub type Domain767 = crate::CooleyTukeyDomain<OuterScalar>;

/// Keyset of BLS12-381 public keys, interpolated over BW6-767's scalar field.
///
/// Independent of the polynomial commitment scheme.
pub type Keyset381 = Keyset<InnerCurve, OuterCurve, Domains767>;

/// Accountable public input for simple and packed proof schemes
///
/// Contains:
/// - Aggregate public key (APK) on the inner curve
/// - Bitmask identifying which keys participated
pub type AccountablePublicInput381 = AccountablePublicInput<InnerCurve>;

/// Counting public input for counting proof scheme
///
/// Contains:
/// - Aggregate public key (APK) on the inner curve  
/// - Count of participating keys (instead of full bitmask)
pub type CountingPublicInput381 = CountingPublicInput<InnerCurve>;

// ============================================================================
// Subgroup Membership
// ============================================================================

/// BW6-767 G1 membership by multiplication by `r`: arkworks'
/// `is_in_correct_subgroup_assuming_on_curve`, after checking the point is on the curve.
///
/// This is the test the curve's authors' implementation uses. BW6-767 was generated by El
/// Housni and Guillevic (<https://hackmd.io/@gnark/bw6_bls12381>, ePrint 2021/1359), and El
/// Housni's implementation lives in gnark-crypto, on the branch `feat/bw6_on_bls12-381` of his
/// fork, pinned here at commit df9b9fe (2021-08-11) and not merged into Consensys/gnark-crypto:
/// - `G1Affine.IsInSubGroup`, on the curve and then the Jacobian test:
///   <https://github.com/yelhousni/gnark-crypto/blob/df9b9feabb6ada6024e5c4648ca54cb968bd69bf/ecc/bw6-767/g1.go#L137-L141>
/// - `G1Jac.IsInSubGroup`, `[r]P` is the identity:
///   <https://github.com/yelhousni/gnark-crypto/blob/df9b9feabb6ada6024e5c4648ca54cb968bd69bf/ecc/bw6-767/g1.go#L369-L376>
///
/// `ark-bw6-767` 0.6 provides no endomorphism (no `GLVConfig`) for G1, so the faster test used
/// for BW6-761 in [`crate::endo`] has no ready-made constants here.
impl crate::PrimeSubgroup for ark_ec::short_weierstrass::Projective<ark_bw6_767::g1::Config> {
    fn is_in_prime_subgroup(p: &BW6_767_G1Affine) -> bool {
        crate::generic_subgroup_check(p)
    }
}

/// arkworks' own check, which for BLS12-381 G1 is already endomorphism-based: `ark-bls12-381`
/// overrides `is_in_correct_subgroup_assuming_on_curve` with Scott's test, checking
/// `φ(P) = -[x^2]P`; Section 6 of <https://eprint.iacr.org/2021/1130>.
impl crate::PrimeSubgroup for ark_ec::short_weierstrass::Projective<ark_bls12_381::g1::Config> {
    fn is_in_prime_subgroup(p: &ark_bls12_381::G1Affine) -> bool {
        crate::generic_subgroup_check(p)
    }
}

/// Seeds the affine-addition accumulator at `(0, 2)`: on the curve `y^2 = x^3 + 4`,
/// of order 3, and so outside G1. See [`crate::AccumulatorSeed`].
// Spelled as the concrete projective type rather than the `InnerCurve` alias: the alias goes
// through `Bls12Config::G1Config`, and coherence cannot tell the two curves' projections apart.
impl crate::AccumulatorSeed for ark_ec::short_weierstrass::Projective<ark_bls12_381::g1::Config> {
    fn accumulator_seed() -> ark_ec::short_weierstrass::Affine<ark_bls12_381::g1::Config> {
        crate::point_in_g1_complement::<ark_bls12_381::g1::Config>()
    }
}

#[cfg(test)]
mod prime_subgroup {
    use crate::PrimeSubgroup;
    use ark_ec::CurveGroup;

    /// Accepts G1, and rejects points of every small order BW6-767's cofactor admits: it is
    /// 2^2 * 3 * 1801 * 10429 * (a 358-bit number with no factor below 2^21).
    #[test]
    fn bw6_767_check_is_exact() {
        crate::test_helpers::check_exact::<ark_bw6_767::g1::Config>(
            |p| super::OuterCurve::is_in_prime_subgroup(&p.into_affine()),
            &[2, 3, 1801, 10429],
        );
    }
}

#[cfg(test)]
mod accumulator_seed {
    use super::*;
    use crate::AccumulatorSeed;
    use ark_ec::AffineRepr;
    use ark_std::Zero;

    /// The seed the prover and verifier start the accumulator from must be a curve point that
    /// is not in G1 — pinned by value, so changing it is a deliberate protocol change.
    #[test]
    fn seed_is_on_the_curve_and_outside_g1() {
        let h = InnerCurve::accumulator_seed();
        assert_eq!(
            h.xy(),
            Some((ark_bls12_381::Fq::from(0u8), ark_bls12_381::Fq::from(2u8)))
        );
        assert!(h.is_on_curve());
        assert!(!h.is_in_correct_subgroup_assuming_on_curve());
        // Order 3: h + h + h is the identity, and h itself is not.
        assert!(!h.is_zero());
        assert!((h.into_group() + h + h).is_zero());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::subgroup_generator;
    use crate::{DomainSizes, DomainTriple};
    use num_bigint::BigUint;

    /// The list is the reserved-factor rule made explicit: every entry divides
    /// `3 * 11 * 23 * 47 * 10177`, and every divisor of it is an entry bar the two documented
    /// omissions. Pinning the omissions here is what keeps them a deliberate choice rather than
    /// a typo nobody would notice.
    #[test]
    fn table_is_the_divisors_of_the_unreserved_part_less_the_dominated_ones() {
        const DOMINATED: [usize; 2] = [10177, 30531];

        let unreserved: usize = 3 * 11 * 23 * 47 * 10177;
        let mut divisors: Vec<usize> = vec![1];
        for f in [3usize, 11, 23, 47, 10177] {
            divisors = divisors.iter().flat_map(|&d| [d, d * f]).collect();
        }
        divisors.sort_unstable();
        divisors.dedup();
        assert_eq!(divisors.len(), 32);

        let expected: Vec<usize> = divisors
            .iter()
            .copied()
            .filter(|n| !DOMINATED.contains(n))
            .collect();
        assert_eq!(APK381_DOMAIN_SIZES, &expected[..]);
        assert!(APK381_DOMAIN_SIZES.iter().all(|&n| unreserved % n == 0));
        assert_eq!(*APK381_DOMAIN_SIZES.last().unwrap(), unreserved);

        // Each omitted size really is dominated: the next entry up is larger, so no request is
        // left unserved, and it is cheaper, which is the reason to pass the smaller one over.
        for d in DOMINATED {
            let next = APK381_DOMAIN_SIZES
                .iter()
                .copied()
                .find(|&n| n > d)
                .expect("a dominated entry must have a successor");
            assert!(
                next < 2 * d,
                "{} would double the padding, not a fair swap",
                next
            );
        }
    }

    /// The reason the 2 and the 3 are held back: `6n` has to exist for every entry, and it has
    /// to be the whole usable smooth part of `q - 1` at the top of the table.
    #[test]
    fn every_entry_admits_its_doubled_and_sextupled_domain() {
        let order: BigUint =
            Into::<BigUint>::into(<OuterScalar as ark_ff::PrimeField>::MODULUS) - 1u8;
        for &n in APK381_DOMAIN_SIZES {
            assert!(
                (&order % BigUint::from(6 * n)) == BigUint::from(0u32),
                "6 * {} does not divide q - 1",
                n
            );
        }
        assert_eq!(
            6 * APK381_DOMAIN_SIZES.last().unwrap(),
            2 * 9 * 11 * 23 * 47 * 10177,
            "the largest triple should exhaust the smooth part of q - 1"
        );
        // ...and every size `triple_for` hands out really is a subgroup order of this field.
        // That is the `DomainSizes` contract, which `SmoothDomainSet::build` only
        // `debug_assert`s so the production path stays a binary search. Asserting it here is
        // what earns that.
        for &n in APK381_DOMAIN_SIZES {
            let t = Apk381DomainSizes::triple_for(n).expect("an entry must resolve to itself");
            for size in [t.base, t.medium, t.large] {
                assert!(
                    subgroup_generator::<OuterScalar>(size).is_some(),
                    "{} is not a subgroup order of BW6-767's scalar field",
                    size
                );
            }
        }
    }

    /// The expansion rule, which is this module's choice and not the protocol's: every entry
    /// resolves to `n, 2n, 6n`, and every such triple clears the protocol's floors.
    #[test]
    fn every_entry_expands_to_n_2n_6n() {
        for &n in APK381_DOMAIN_SIZES {
            let t = Apk381DomainSizes::triple_for(n).expect("an entry must resolve to itself");
            assert_eq!(t, DomainTriple::new(n, 2 * n, 6 * n));
            assert!(t.meets_protocol_bounds(), "{:?}", t);
        }
    }

    /// A request between entries rounds up to the next one; one past the top has no answer.
    #[test]
    fn triple_for_rounds_up_and_runs_out() {
        assert_eq!(Apk381DomainSizes::triple_for(254).unwrap().base, 517);
        assert_eq!(Apk381DomainSizes::triple_for(1).unwrap().base, 1);
        let past_the_end = APK381_DOMAIN_SIZES.last().unwrap() + 1;
        assert_eq!(Apk381DomainSizes::triple_for(past_the_end), None);
    }

    /// No size this configuration uses is a multiple of 4, which is why the `packed` scheme is
    /// unavailable on APK-381: it splits the bitmask into 256-bit chunks and needs `256 | n`.
    #[test]
    fn no_size_is_a_multiple_of_four() {
        for &n in APK381_DOMAIN_SIZES {
            let t = Apk381DomainSizes::triple_for(n).unwrap();
            for size in [t.base, t.medium, t.large] {
                assert_ne!(size % 4, 0, "{} is a multiple of 4", size);
            }
        }
    }
}
