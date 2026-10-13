use merlin::{Transcript, TranscriptRng};
use rand::{CryptoRng, Error, RngCore};

/// An "rng" that only yields zeros, for [`fiat_shamir_rng`], which must be deterministic.
/// Merlin's `TranscriptRngBuilder::finalize` only calls `fill_bytes`, so the rest is unreachable.
struct DummyRng;

impl RngCore for DummyRng {
    fn next_u32(&mut self) -> u32 {
        unimplemented!()
    }

    fn next_u64(&mut self) -> u64 {
        unimplemented!()
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        dest.iter_mut().for_each(|byte| *byte = 0u8);
    }

    fn try_fill_bytes(&mut self, _dest: &mut [u8]) -> Result<(), Error> {
        unimplemented!()
    }
}

impl CryptoRng for DummyRng {}

/// A deterministic rng bound to the transcript, for the randomness of batched opening
/// verification.
///
/// Merlin's `TranscriptRng` is designed for provers, who mix in secret witness bytes and
/// external randomness (<https://merlin.cool/transcript/rng.html>). The verifier has no secret
/// and needs none: the randomness only has to be unpredictable to the prover when it fixes the
/// proof, which binding it to the transcript already ensures. So it rekeys with a fixed
/// placeholder and finalizes with [`DummyRng`]'s zeros, and the stream is a function of the
/// transcript alone.
pub fn fiat_shamir_rng(transcript: &mut Transcript) -> TranscriptRng {
    transcript
        .build_rng()
        .rekey_with_witness_bytes(b"verifier_secret", &[42]) //TODO: Does verifier know secrets?
        .finalize(&mut DummyRng)
}
