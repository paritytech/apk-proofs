//! Errors returned by keyset construction, proving and verification.
//!
//! Everything a caller can get wrong, and everything that arrives from outside the process — a
//! public key read off the chain, a bitmask, a keyset commitment — surfaces here rather than as a
//! panic, so a long-running service can log the failure and move on to the next block.

use crate::domain::DomainError;

/// Why a public key was rejected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PublicKeyFault {
    /// The point at infinity. It has no affine coordinates, so it cannot be interpolated.
    Identity,
    /// Not in the prime-order subgroup G1, or not on the curve at all. The accumulator's
    /// soundness argument, and the guarantee that it never reaches the identity or the
    /// doubling case, rely on every key being in G1; see [`crate::AccumulatorSeed`].
    NotInSubgroup,
}

impl core::fmt::Display for PublicKeyFault {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            PublicKeyFault::Identity => write!(f, "is the point at infinity"),
            PublicKeyFault::NotInSubgroup => write!(f, "is not in the prime-order subgroup"),
        }
    }
}

/// Why a keyset could not be built or committed to, or a proof could not be produced or checked.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ApkError {
    /// No evaluation domain fits: the keyset is larger than the field supports, or a keyset
    /// commitment names a size the field cannot realise.
    Domain(DomainError),
    /// The public key at `index` cannot be used.
    InvalidPublicKey { index: usize, fault: PublicKeyFault },
    /// The commitment-scheme parameters do not reach the degree this keyset needs. Regenerate
    /// them for at least as many keys as the keyset holds.
    SrsTooSmall {
        domain_size: usize,
        required_degree: usize,
        available_degree: usize,
    },
    /// The bitmask does not have one bit per key.
    BitmaskLengthMismatch { bitmask: usize, keyset: usize },
    /// The bitmask selects no key. The aggregate would be the identity, which has no affine
    /// representation.
    NoSigners,
    /// The commitment scheme failed to commit or open. Upstream reports no detail; the payload
    /// names the step.
    Pcs(&'static str),
    /// The prover's witness does not satisfy the constraints, so the quotient polynomial does
    /// not exist. With a validated keyset and bitmask this indicates a bug.
    ConstraintsNotSatisfied,
    /// A keyset commitment whose sizes cannot go together: it must hold at least one key, and
    /// fewer keys than the domain has rows, the last row being reserved.
    InvalidKeysetCommitment { keyset_size: u64, domain_size: u64 },
    /// A keyset commitment whose points are not in the outer curve's G1, or not on it.
    KeysetCommitmentNotInG1,
    /// A counting public input claiming no signers, or more than the keyset holds.
    CountOutOfRange { count: usize, keyset_size: usize },
    /// A public input the verifier cannot evaluate; the payload says which part.
    InvalidPublicInput(&'static str),
    /// The configuration's accumulator seed has no affine form. It is a constant of the inner
    /// curve, so this is a misconfigured curve, not bad input.
    InvalidAccumulatorSeed,
}

impl core::fmt::Display for ApkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ApkError::Domain(e) => write!(f, "{}", e),
            ApkError::InvalidPublicKey { index, fault } => {
                write!(f, "public key at index {} {}", index, fault)
            }
            ApkError::SrsTooSmall {
                domain_size,
                required_degree,
                available_degree,
            } => write!(
                f,
                "SRS too small: a keyset over a domain of {} needs commitments up to degree {}, \
                 but these parameters stop at {}",
                domain_size, required_degree, available_degree
            ),
            ApkError::BitmaskLengthMismatch { bitmask, keyset } => write!(
                f,
                "bitmask has {} bits but the keyset has {} keys",
                bitmask, keyset
            ),
            ApkError::NoSigners => write!(f, "bitmask selects no signer"),
            ApkError::Pcs(step) => write!(f, "commitment scheme failed: {}", step),
            ApkError::ConstraintsNotSatisfied => {
                write!(f, "witness does not satisfy the constraints")
            }
            ApkError::InvalidKeysetCommitment {
                keyset_size,
                domain_size,
            } => write!(
                f,
                "keyset commitment claims {} keys over a domain of {}; needs at least one key \
                 and fewer keys than the domain size",
                keyset_size, domain_size
            ),
            ApkError::KeysetCommitmentNotInG1 => {
                write!(f, "keyset commitment is not in the outer curve's G1")
            }
            ApkError::CountOutOfRange { count, keyset_size } => write!(
                f,
                "public input claims {} signers; must be between 1 and the keyset size {}",
                count, keyset_size
            ),
            ApkError::InvalidPublicInput(what) => write!(f, "invalid public input: {}", what),
            ApkError::InvalidAccumulatorSeed => {
                write!(f, "the accumulator seed is the point at infinity")
            }
        }
    }
}

impl std::error::Error for ApkError {}

impl From<DomainError> for ApkError {
    fn from(e: DomainError) -> Self {
        ApkError::Domain(e)
    }
}
