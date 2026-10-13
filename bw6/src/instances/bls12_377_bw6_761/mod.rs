//! APK-377: BLS12-377 signatures, proofs over BW6-761.
//!
//! This module provides type aliases and constants for APK proofs using:
//! - **Inner curve**: BLS12-377 G1, where the BLS public keys live.
//! - **Outer curve**: BW6-761 G1, where the proof's commitments live. Both the prover and the
//!   verifier do arithmetic here.
//!
//! The two form a 2-chain: BW6-761's scalar field is BLS12-377's base field, so the inner
//! curve's coordinates are native field elements of the proof system. Ref "Optimized and secure pairing-friendly
//! elliptic curves suitable for one layer proof composition",
//! <https://eprint.iacr.org/2020/351>.
//!
//! ## Polynomial Commitment Schemes
//!
//! - [`kzg`]: KZG commitments on BW6-761, the only scheme implemented.

use crate::{AccountablePublicInput, CountingPublicInput, Keyset};
use ark_bls12_377::G1Projective as Bls12_377_G1;
use ark_bw6_761::{Fq, Fr, G1Affine as BW6_761_G1Affine, G1Projective as BW6_761_G1};
use ark_ec::bls12::Bls12Config;
use ark_ff::MontFp;

// ============================================================================
// Polynomial Commitment Schemes
// ============================================================================

/// KZG polynomial commitment scheme types for this pairing
pub mod kzg;

// ============================================================================
// Curve Type Aliases
// ============================================================================

/// Inner curve: BLS12-377 G1
///
/// Used for:
/// - BLS signature public keys
/// - BLS signature aggregate public keys
/// - Elements being committed to in the keyset
pub type InnerCurve = Bls12_377_G1;

/// The pairing the inner curve belongs to. Needed to sign with these keys; see [`crate::bls`].
pub type InnerPairing = ark_bls12_377::Bls12_377;

/// Outer curve: BW6-761 G1, in projective form.
///
/// The group the commitments and opening proofs live in, used by both sides of the protocol:
/// - the prover commits to its polynomials and computes the opening proofs here;
/// - the verifier checks every proof point for membership in G1 (see
///   [`crate::PrimeSubgroup`]), rebuilds the commitment to the linearization polynomial and
///   aggregates the commitments opened at `zeta` as linear combinations here, and hands the
///   result to the KZG pairing check on BW6-761.
///
/// The projective form is what the group arithmetic uses; points are stored and serialized in
/// [`OuterAffine`] form.
pub type OuterCurve = BW6_761_G1;

/// Outer curve: BW6-761 G1, in affine form.
///
/// The representation of every outer-curve point that is stored or serialized: keyset and
/// register commitments, the quotient commitment, and KZG opening proofs.
pub type OuterAffine = BW6_761_G1Affine;

/// Outer curve scalar field: BW6-761 Fr.
///
/// Equal to BLS12-377's base field `Fq`, so BLS12-377 point coordinates are elements of it.
/// All polynomial arithmetic in the proof system happens over this field.
pub type OuterScalar = Fr;

// ============================================================================
// PCS-Independent Type Aliases
// ============================================================================

/// The evaluation domains for BW6-761's scalar field.
///
/// Two-adicity is 46 there, so power-of-two domains exist for any size the prover could afford
/// and the classical `n, 2n, 4n` layout applies.
pub type Domains761 = crate::Radix2DomainSet<OuterScalar>;

/// A single radix-2 domain. Prefer [`Domains761`]: the protocol works over the triple.
pub type Domain761 = crate::Radix2Domain<OuterScalar>;

/// Keyset of BLS12-377 public keys, interpolated over BW6-761's scalar field.
///
/// Independent of the polynomial commitment scheme.
pub type Keyset377 = Keyset<InnerCurve, OuterCurve, Domains761>;

/// Accountable public input for simple and packed proof schemes
///
/// Contains:
/// - Aggregate public key (APK) on the inner curve
/// - Bitmask identifying which keys participated
pub type AccountablePublicInput377 = AccountablePublicInput<InnerCurve>;

/// Counting public input for counting proof scheme
///
/// Contains:
/// - Aggregate public key (APK) on the inner curve  
/// - Count of participating keys (instead of full bitmask)
pub type CountingPublicInput377 = CountingPublicInput<InnerCurve>;

// ============================================================================
// Endomorphism Constants
// ============================================================================

/// Eigenvalue `λ` of the endomorphism `φ: (x, y) ↦ (ωx, y)` on BW6-761 G1 (not BLS12-377):
/// `φ(P) = [λ]P` for every `P` in G1, with `ω` = [`OMEGA`]. An element of BW6-761's scalar
/// field, and a primitive cube root of unity there.
///
/// Only the tests use it, to check [`OMEGA`] against it. The same `(ω, λ)` pair is in the
/// Zexe BW6-761 GLV parameters,
/// <https://github.com/celo-org/zexe/blob/master/algebra/src/bw6_761/curves/g1.rs#L37-L71>.
pub const LAMBDA: Fr = MontFp!(
    "80949648264912719408558363140637477264845294720710499478137287262712535938301461879813459410945"
);

/// The BLS12-377 seed `u` (arkworks' `X`), as little-endian limbs.
///
/// BW6-761 is built from BLS12-377, so its group order is a polynomial in the same seed. The
/// G1 membership test in [`crate::endo::subgroup_check`] multiplies by `u` three times.
pub const U: &[u64] = ark_bls12_377::Config::X;

/// A primitive cube root of unity `ω` in BW6-761's base field `Fq`, defining the endomorphism
/// `φ: (x, y) ↦ (ωx, y)` of BW6-761 G1. On G1 it acts as multiplication by [`LAMBDA`].
///
/// Used by the G1 membership test in [`crate::endo::subgroup_check`]. Same value as `OMEGA` in
/// the Zexe BW6-761 GLV parameters,
/// <https://github.com/celo-org/zexe/blob/master/algebra/src/bw6_761/curves/g1.rs#L37-L71>.
pub const OMEGA: Fq = MontFp!(
    "196898582409020929727861073970057715139766638230382572845074161156680037021882725775086501\
     3421937292370006175842381275743914023380727582819905021229583192207421122272650305267822868\
     639090213645505120388400344940985710520836292650"
);

/// BW6-761 G1 membership through the endomorphism, rather than arkworks' default
/// multiplication by `r`. See [`crate::endo::subgroup_check`].
impl crate::PrimeSubgroup for ark_ec::short_weierstrass::Projective<ark_bw6_761::g1::Config> {
    fn is_in_prime_subgroup(p: &BW6_761_G1Affine) -> bool {
        use ark_ec::AffineRepr;
        p.is_on_curve()
            && crate::endo::subgroup_check::<ark_bw6_761::Config>(&p.into_group(), OMEGA, U)
    }
}

/// arkworks' own check: double-and-add by `r`. `ark-bls12-377` 0.6 does not override
/// `is_in_correct_subgroup_assuming_on_curve` for G1, so the default applies.
impl crate::PrimeSubgroup for ark_ec::short_weierstrass::Projective<ark_bls12_377::g1::Config> {
    fn is_in_prime_subgroup(p: &ark_bls12_377::G1Affine) -> bool {
        crate::generic_subgroup_check(p)
    }
}

/// Seeds the affine-addition accumulator at `(0, 1)`: on the curve `y^2 = x^3 + 1`,
/// of order 3, and so outside G1. See [`crate::AccumulatorSeed`].
// Spelled as the concrete projective type rather than the `InnerCurve` alias: the alias goes
// through `Bls12Config::G1Config`, and coherence cannot tell the two curves' projections apart.
impl crate::AccumulatorSeed for ark_ec::short_weierstrass::Projective<ark_bls12_377::g1::Config> {
    fn accumulator_seed() -> ark_ec::short_weierstrass::Affine<ark_bls12_377::g1::Config> {
        crate::point_in_g1_complement::<ark_bls12_377::g1::Config>()
    }
}

#[cfg(test)]
mod prime_subgroup {
    use crate::PrimeSubgroup;
    use ark_ec::CurveGroup;

    /// The endomorphism test accepts G1, and rejects points of every small order BW6-761's
    /// cofactor admits: it is 2^2 * 127 * (a 375-bit number with no factor below 2^21).
    #[test]
    fn bw6_761_check_is_exact() {
        crate::test_helpers::check_exact::<ark_bw6_761::g1::Config>(
            |p| super::OuterCurve::is_in_prime_subgroup(&p.into_affine()),
            &[2, 127],
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
            Some((ark_bls12_377::Fq::from(0u8), ark_bls12_377::Fq::from(1u8)))
        );
        assert!(h.is_on_curve());
        assert!(!h.is_in_correct_subgroup_assuming_on_curve());
        // Order 3: h + h + h is the identity, and h itself is not.
        assert!(!h.is_zero());
        assert!((h.into_group() + h + h).is_zero());
    }
}
