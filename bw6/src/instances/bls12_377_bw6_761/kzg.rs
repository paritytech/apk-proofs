//! KZG polynomial commitment scheme for BLS12-377 + BW6-761
//!
//! Type aliases for APK proofs using KZG (Kate-Zaverucha-Goldberg) polynomial commitments on
//! BW6-761, as implemented by `w3f-pcs` (<https://github.com/paritytech/fflonk>).
//!
//! - **Trusted setup**: the parameters embed powers of a secret `tau` and must come from a
//!   setup ceremony; see [`crate::setup`].
//! - **Commitments and opening proofs**: one BW6-761 G1 point each.
//! - **Verification**: the verifier's two opening claims (at `zeta` and `zeta * omega`) are
//!   batched into one pairing-product check by `w3f-pcs`'s `batch_verify`.

use ark_bw6_761::BW6_761;
use w3f_pcs::pcs::commitment::WrappedAffine;
use w3f_pcs::pcs::kzg::KZG;

use super::*;
use crate::{CountingProof, KeysetCommitment, PackedProof, Prover, SimpleProof, Verifier};

// ============================================================================
// KZG PCS Type Aliases
// ============================================================================

/// KZG polynomial commitment scheme on BW6-761
pub type Pcs = KZG<BW6_761>;

/// KZG commitment (a single BW6-761 G1 point)
pub type Commitment = WrappedAffine<OuterCurve>;

// ============================================================================
// Core Types with KZG
// ============================================================================

/// Keyset commitment using KZG on BW6-761
///
/// Commitments to the two polynomials interpolating the x and y coordinates of the public
/// keys, plus the domain size and key count the verifier needs.
pub type KeysetCommitment377 = KeysetCommitment<OuterScalar, Commitment>;

/// Prover for BLS12-377 + BW6-761 with KZG commitments
pub type Prover377 = Prover<InnerCurve, OuterCurve, Pcs, super::Domains761>;

/// Verifier for BLS12-377 + BW6-761 with KZG commitments
pub type Verifier377 = Verifier<InnerCurve, OuterCurve, Pcs, super::Domains761>;

// ============================================================================
// Proof Type Aliases
// ============================================================================

/// The 'basic' accountable scheme's proof: the bitmask is public and the verifier evaluates it
/// itself.
pub type SimpleProof377 = SimpleProof<OuterScalar, OuterAffine, Commitment, OuterAffine>;

/// The 'packed' accountable scheme's proof: the prover commits to the bitmask, and the verifier
/// checks it against the public bitmask packed into 256-bit chunks: one field operation per
/// chunk instead of per bit. See the "packed" scheme in <https://eprint.iacr.org/2022/1205>.
pub type PackedProof377 = PackedProof<OuterScalar, OuterAffine, Commitment, OuterAffine>;

/// The 'counting' scheme's proof: the public input is the number of signers, not the bitmask,
/// so the proof shows that many keys were aggregated without saying which.
pub type CountingProof377 = CountingProof<OuterScalar, OuterAffine, Commitment, OuterAffine>;
