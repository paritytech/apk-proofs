use ark_ec::CurveGroup;
use ark_ff::PrimeField;
use ark_serialize::CanonicalSerialize;
use merlin::Transcript;
use w3f_pcs::pcs::{Commitment, RawVerifierKey};

use crate::domain::FftDomain;
use crate::piop::{RegisterCommitments, RegisterEvaluations};
use crate::{KeysetCommitment, PublicInput};

pub(crate) trait ApkTranscript<F: PrimeField> {
    fn set_protocol_params<D: FftDomain<F>, VK: RawVerifierKey>(
        &mut self,
        domain: &D,
        verifier_key: &VK,
    ) {
        // Bind the domain by size and generator rather than by serialising a concrete domain
        // type. Sizes are not powers of two on BW6-767, and prover and verifier disagreeing
        // about which domain a proof is over must be a verification failure, not silent.
        self._append_serializable(b"domain_size", &(domain.size() as u64));
        self._append_serializable(b"domain_generator", &domain.generator());
        self._append_serializable(b"vk", verifier_key);
    }

    fn set_keyset_commitment<C>(&mut self, keyset_commitment: &KeysetCommitment<F, C>)
    where
        C: Commitment<F>,
    {
        self._append_serializable(b"keyset_commitment", keyset_commitment);
    }

    fn append_public_input<IC: CurveGroup>(&mut self, public_input: &impl PublicInput<IC>) {
        self._append_serializable(b"public_input", public_input);
    }

    fn append_register_commitments<G: ark_ec::AffineRepr>(
        &mut self,
        register_commitments: &impl RegisterCommitments<G>,
    ) {
        self._append_serializable(b"register_commitments", register_commitments);
    }

    fn get_bitmask_aggregation_challenge(&mut self) -> F {
        self._get_128_bit_challenge(b"bitmask_aggregation")
    }

    fn append_2nd_round_register_commitments<G: ark_ec::AffineRepr>(
        &mut self,
        register_commitments: &impl RegisterCommitments<G>,
    ) {
        self._append_serializable(b"2nd_round_register_commitments", register_commitments);
    }

    fn get_constraints_aggregation_challenge(&mut self) -> F {
        self._get_128_bit_challenge(b"constraints_aggregation")
    }

    fn append_quotient_commitment<C: Commitment<F>>(&mut self, commitment: &C) {
        self._append_serializable(b"quotient", commitment);
    }

    fn get_evaluation_point(&mut self) -> F {
        self._get_128_bit_challenge(b"evaluation_point")
    }

    fn append_evaluations(
        &mut self,
        evals: &impl RegisterEvaluations<F>,
        q_at_zeta: &F,
        r_at_zeta_omega: &F,
    ) {
        self._append_serializable(b"register_evaluations", evals);
        self._append_serializable(b"quotient_evaluation", q_at_zeta);
        self._append_serializable(b"shifted_linearization_evaluation", r_at_zeta_omega);
    }

    fn get_opening_aggregation_challenges(&mut self, n: usize) -> Vec<F> {
        self._get_128_bit_challenges(b"opening_aggregation", n)
    }

    fn append_opening_proofs<P: CanonicalSerialize>(&mut self, at_zeta: &P, at_zeta_omega: &P) {
        self._append_serializable(b"opening_proof_at_zeta", at_zeta);
        self._append_serializable(b"opening_proof_at_zeta_omega", at_zeta_omega);
    }

    fn _get_128_bit_challenge(&mut self, label: &'static [u8]) -> F;

    fn _get_128_bit_challenges(&mut self, label: &'static [u8], n: usize) -> Vec<F>;

    fn _append_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize);
}

impl<F: PrimeField> ApkTranscript<F> for Transcript {
    fn _get_128_bit_challenge(&mut self, label: &'static [u8]) -> F {
        let mut buf = [0u8; 16];
        self.challenge_bytes(label, &mut buf);
        #[allow(clippy::expect_used, reason = "invariant argued in the message")]
        let challenge = F::from_random_bytes(&buf)
            .expect("invariant: 128 random bits always fit below a modulus this large");
        challenge
    }

    fn _get_128_bit_challenges(&mut self, label: &'static [u8], n: usize) -> Vec<F> {
        (0..n)
            .map(|_| <Self as ApkTranscript<F>>::_get_128_bit_challenge(self, label))
            .collect() //TODO: unlikely secure
    }

    fn _append_serializable(&mut self, label: &'static [u8], message: &impl CanonicalSerialize) {
        let mut buf = vec![0; message.compressed_size()];
        #[allow(clippy::expect_used, reason = "invariant argued in the message")]
        message
            .serialize_compressed(&mut buf)
            .expect("invariant: the buffer is exactly `compressed_size` bytes");
        self.append_message(label, &buf);
    }
}
