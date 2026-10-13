//! KZG polynomial commitment scheme for BLS12-381 + BW6-767
//!
//! Type aliases for APK proofs using KZG (Kate-Zaverucha-Goldberg) polynomial commitments on
//! BW6-767, as implemented by `w3f-pcs` (<https://github.com/paritytech/fflonk>).
//!
//! - **Trusted setup**: the parameters embed powers of a secret `tau` and must come from a
//!   setup ceremony; see [`crate::setup`].
//! - **Commitments and opening proofs**: one BW6-767 G1 point each, 97 bytes compressed.
//! - **Verification**: the verifier's two opening claims (at `zeta` and `zeta * omega`) are
//!   batched into one pairing-product check by `w3f-pcs`'s `batch_verify`.

use ark_bw6_767::BW6_767;
use w3f_pcs::pcs::commitment::WrappedAffine;
use w3f_pcs::pcs::kzg::KZG;

use super::*;
use crate::{CountingProof, KeysetCommitment, PackedProof, Prover, SimpleProof, Verifier};

// ============================================================================
// KZG PCS Type Aliases
// ============================================================================

/// KZG polynomial commitment scheme on BW6-767
pub type Pcs = KZG<BW6_767>;

/// KZG commitment (a single BW6-767 G1 point)
pub type Commitment = WrappedAffine<OuterCurve>;

// ============================================================================
// Core Types with KZG
// ============================================================================

/// Keyset commitment using KZG on BW6-767
///
/// Commitments to the two polynomials interpolating the x and y coordinates of the public
/// keys, plus the domain size and key count the verifier needs.
pub type KeysetCommitment381 = KeysetCommitment<OuterScalar, Commitment>;

/// Prover for BLS12-381 + BW6-767 with KZG commitments
pub type Prover381 = Prover<InnerCurve, OuterCurve, Pcs, super::Domains767>;

/// Verifier for BLS12-381 + BW6-767 with KZG commitments
pub type Verifier381 = Verifier<InnerCurve, OuterCurve, Pcs, super::Domains767>;

// ============================================================================
// Proof Type Aliases
// ============================================================================

/// The 'basic' accountable scheme's proof: the bitmask is public and the verifier evaluates it
/// itself.
pub type SimpleProof381 = SimpleProof<OuterScalar, OuterAffine, Commitment, OuterAffine>;

/// The 'packed' scheme's proof type, spelled for this configuration.
///
/// **No such proof can be produced on APK-381**: `packed` needs `256 | n`, and no BW6-767
/// domain size is even a multiple of 4. `prove_packed` does not compile for this
/// configuration; see [`crate::SupportsPackedScheme`]. The alias exists only for symmetry with
/// APK-377.
pub type PackedProof381 = PackedProof<OuterScalar, OuterAffine, Commitment, OuterAffine>;

/// The 'counting' scheme's proof: the public input is the number of signers, not the bitmask,
/// so the proof shows that many keys were aggregated without saying which.
pub type CountingProof381 = CountingProof<OuterScalar, OuterAffine, Commitment, OuterAffine>;
