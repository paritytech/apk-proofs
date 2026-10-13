use ark_ec::CurveGroup;
use ark_ff::FftField;
use ark_std::{end_timer, start_timer};
use merlin::{Transcript as MerlinTranscript, TranscriptRng};
use w3f_pcs::aggregation::single::aggregate_claims_multiexp;
use w3f_pcs::pcs::{PcsParams, RawVerifierKey, PCS};

use crate::domain::{DomainSet, FftDomain};
use crate::fsrng::fiat_shamir_rng;
use crate::piop::affine_addition::AffineAdditionEvaluations;
use crate::piop::bitmask_packing::SuccinctAccountableRegisterEvaluations;
use crate::piop::counting::CountingEvaluations;
use crate::piop::{RegisterCommitments, RegisterEvaluations, VerifierProtocol};
use crate::transcript::ApkTranscript;
use crate::utils::LagrangeEvaluations;
use crate::{
    utils, AccountablePublicInput, AccumulatorSeed, ApkError, CommitmentExt, CountingProof,
    CountingPublicInput, KeysetCommitment, OpeningProofPoints, PackedProof, PrimeSubgroup, Proof,
    PublicInput, SimpleProof,
};

type Transcript = MerlinTranscript;

/// The verifier challenges, recomputed from the transcript.
pub struct Challenges<F: FftField> {
    /// Bitmask-chunk aggregation challenge (used by 'packed' only).
    pub r: F,
    /// Constraint aggregation challenge.
    pub phi: F,
    /// Evaluation point.
    pub zeta: F,
    /// Opening aggregation challenges, one per polynomial opened at `zeta`.
    pub nus: Vec<F>,
}

/// Checks proofs against one keyset commitment. Construction validates the commitment and
/// binds it into the transcript once, for every proof that follows.
pub struct Verifier<IC, OC, S, D>
where
    IC: AccumulatorSeed,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField> + FftField,
    S: PCS<OC::ScalarField>,
    D: DomainSet<OC::ScalarField>,
{
    domain: D::Domain,
    verifier_key: <S::Params as PcsParams>::RVK,
    pks_comm: KeysetCommitment<OC::ScalarField, S::C>,
    preprocessed_transcript: Transcript,
    _marker: std::marker::PhantomData<(IC, S)>,
}

impl<IC, OC, S, D> Verifier<IC, OC, S, D>
where
    IC: AccumulatorSeed,
    OC: CurveGroup + PrimeSubgroup,
    OC::ScalarField: From<IC::BaseField> + FftField,
    S: PCS<OC::ScalarField>,
    S::C: CommitmentExt<OC::ScalarField, Affine = OC::Affine> + Clone,
    S::Proof: OpeningProofPoints<OC::Affine>,
    D: DomainSet<OC::ScalarField>,
{
    /// The sizes come out of the keyset commitment, which for a bridge arrives from the chain.
    /// The commitment is signed by the validators, so the sizes are trusted to be what they
    /// signed — but not to be usable: the domain size may name one this field cannot realise,
    /// or the key count may not fit it. Both are rejected rather than panicking.
    ///
    /// Everything checked here is checked once per keyset rather than once per proof.
    pub fn try_new(
        verifier_key: <S::Params as PcsParams>::RVK,
        pks_comm: KeysetCommitment<OC::ScalarField, S::C>,
        mut empty_transcript: Transcript,
    ) -> Result<Self, ApkError> {
        // At least one key, and every key on a row below the last, which is reserved. The
        // cheap check goes first: building the domain is not free.
        if pks_comm.keyset_size == 0 || pks_comm.keyset_size >= pks_comm.domain_size {
            return Err(ApkError::InvalidKeysetCommitment {
                keyset_size: pks_comm.keyset_size,
                domain_size: pks_comm.domain_size,
            });
        }
        let domain_size = pks_comm.domain_size as usize;
        let domain = D::base_for_exact_size(domain_size)?;
        // The signers vouch for what the commitment says, but a malformed encoding is not
        // something they signed off on; checking here keeps it off the pairing for good.
        if ![&pks_comm.pks_comm.0, &pks_comm.pks_comm.1]
            .iter()
            .all(|c| OC::is_in_prime_subgroup(&c.to_affine()))
        {
            return Err(ApkError::KeysetCommitmentNotInG1);
        }

        <Transcript as ApkTranscript<OC::ScalarField>>::set_protocol_params(
            &mut empty_transcript,
            &domain,
            &verifier_key,
        );
        <Transcript as ApkTranscript<OC::ScalarField>>::set_keyset_commitment(
            &mut empty_transcript,
            &pks_comm,
        );

        Ok(Self {
            domain,
            verifier_key,
            pks_comm,
            preprocessed_transcript: empty_transcript,
            _marker: std::marker::PhantomData,
        })
    }

    /// `Ok(false)` for a proof that does not verify, including one with a curve point outside
    /// G1; `Err` for a public input that is malformed before any proof is looked at.
    pub fn verify_simple(
        &self,
        public_input: &AccountablePublicInput<IC>,
        proof: &SimpleProof<OC::ScalarField, OC::Affine, S::C, S::Proof>,
    ) -> Result<bool, ApkError> {
        self.check_bitmask_length(&public_input.bitmask)?;
        Self::check_apk(&public_input.apk)?;
        if !Self::proof_points_in_g1(proof) {
            return Ok(false);
        }
        let (challenges, mut fsrng) = self.restore_challenges(
            public_input,
            proof,
            <AffineAdditionEvaluations<OC::ScalarField> as VerifierProtocol<IC, OC, S>>::POLYS_OPENED_AT_ZETA
        );
        let evals_at_zeta = utils::lagrange_evaluations(challenges.zeta, &self.domain);

        let t_linear_accountability = start_timer!(|| "linear accountability check");
        let b_at_zeta =
            utils::barycentric_eval_binary_at(challenges.zeta, &public_input.bitmask, &self.domain);
        end_timer!(t_linear_accountability);

        let evaluations_with_bitmask = AffineAdditionEvaluations {
            keyset: proof.register_evaluations.keyset,
            bitmask: b_at_zeta,
            partial_sums: proof.register_evaluations.partial_sums,
        };

        let openings_valid = self.validate_evaluations(
            proof,
            &evaluations_with_bitmask,
            &challenges,
            &mut fsrng,
            &evals_at_zeta,
        );

        let apk = public_input.apk;
        let constraint_polynomial_evals = evaluations_with_bitmask
            .evaluate_constraint_polynomials::<IC, OC>(&apk, &evals_at_zeta)?;
        let w = utils::horner_field(&constraint_polynomial_evals, challenges.phi);
        Ok(openings_valid
            && (proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial))
    }

    /// See [`Self::verify_simple`].
    pub fn verify_packed(
        &self,
        public_input: &AccountablePublicInput<IC>,
        proof: &PackedProof<OC::ScalarField, OC::Affine, S::C, S::Proof>,
    ) -> Result<bool, ApkError> {
        self.check_bitmask_length(&public_input.bitmask)?;
        Self::check_apk(&public_input.apk)?;
        if !Self::proof_points_in_g1(proof) {
            return Ok(false);
        }
        let (challenges, mut fsrng) = self.restore_challenges(
            public_input,
            proof,
            <SuccinctAccountableRegisterEvaluations<OC::ScalarField> as VerifierProtocol<
                IC,
                OC,
                S,
            >>::POLYS_OPENED_AT_ZETA,
        );
        let evals_at_zeta = utils::lagrange_evaluations(challenges.zeta, &self.domain);

        let openings_valid = self.validate_evaluations(
            proof,
            &proof.register_evaluations,
            &challenges,
            &mut fsrng,
            &evals_at_zeta,
        );

        let apk = public_input.apk;
        let constraint_polynomial_evals = proof
            .register_evaluations
            .evaluate_constraint_polynomials::<IC, OC>(
                &apk,
                &evals_at_zeta,
                challenges.r,
                &public_input.bitmask,
                self.domain.size() as u64,
            )?;
        let w = utils::horner_field(&constraint_polynomial_evals, challenges.phi);
        Ok(openings_valid
            && (proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial))
    }

    /// See [`Self::verify_simple`].
    pub fn verify_counting(
        &self,
        public_input: &CountingPublicInput<IC>,
        proof: &CountingProof<OC::ScalarField, OC::Affine, S::C, S::Proof>,
    ) -> Result<bool, ApkError> {
        let keyset_size = self.keyset_size();
        if public_input.count == 0 || public_input.count > keyset_size {
            return Err(ApkError::CountOutOfRange {
                count: public_input.count,
                keyset_size,
            });
        }
        Self::check_apk(&public_input.apk)?;
        if !Self::proof_points_in_g1(proof) {
            return Ok(false);
        }
        let (challenges, mut fsrng) = self.restore_challenges(
            public_input,
            proof,
            <CountingEvaluations<OC::ScalarField> as VerifierProtocol<IC, OC, S>>::POLYS_OPENED_AT_ZETA
        );
        let evals_at_zeta = utils::lagrange_evaluations(challenges.zeta, &self.domain);
        let count = OC::ScalarField::from(public_input.count as u32);

        let openings_valid = self.validate_evaluations(
            proof,
            &proof.register_evaluations,
            &challenges,
            &mut fsrng,
            &evals_at_zeta,
        );

        let apk = public_input.apk;
        let constraint_polynomial_evals = proof
            .register_evaluations
            .evaluate_constraint_polynomials::<IC, OC>(apk, count, &evals_at_zeta)?;
        let w = utils::horner_field(&constraint_polynomial_evals, challenges.phi);
        Ok(openings_valid
            && (proof.r_zeta_omega + w == proof.q_zeta * evals_at_zeta.vanishing_polynomial))
    }

    fn keyset_size(&self) -> usize {
        self.pks_comm.keyset_size as usize
    }

    /// One bit per key, exactly. The verifier folds the bitmask into a polynomial over the
    /// domain, so extra bits would not be ignored: they would land on padding rows, or wrap
    /// around onto real keys and report signers under the wrong index.
    fn check_bitmask_length(&self, bitmask: &crate::Bitmask) -> Result<(), ApkError> {
        if bitmask.size() != self.keyset_size() {
            return Err(ApkError::BitmaskLengthMismatch {
                bitmask: bitmask.size(),
                keyset: self.keyset_size(),
            });
        }
        Ok(())
    }

    /// The identity as an aggregate key: with no signers selected, the accumulator never moves
    /// off its seed, so a valid proof for `apk = 0` exists. BLS verification against the
    /// identity key then accepts the identity signature on any message. A light client's quorum
    /// check rules that out, but not every caller runs one.
    ///
    /// Nothing else about `apk` is checked here: a verifying proof implies it is in G1. See
    /// [`AccountablePublicInput`].
    fn check_apk(apk: &IC::Affine) -> Result<(), ApkError> {
        use ark_ec::AffineRepr;
        if apk.is_zero() {
            return Err(ApkError::InvalidPublicInput("apk is the identity"));
        }
        Ok(())
    }

    /// Whether every curve point in `proof` is in the outer curve's G1: the register and
    /// quotient commitments and the opening proofs.
    ///
    /// Points outside G1 must not reach the pairing, whose soundness argument holds on G1 only.
    /// The old verifier checked the two points it fed the pairing, after batching the openings;
    /// behind the [`PCS`] interface those are out of reach, so each point is checked instead.
    /// That is also the stronger check: a check on a random combination only, with small
    /// factors in the cofactor, can be passed by grinding the challenges until the components
    /// outside G1 cancel.
    fn proof_points_in_g1<E, C, AC>(
        proof: &Proof<OC::ScalarField, E, C, AC, S::C, S::Proof>,
    ) -> bool
    where
        E: RegisterEvaluations<OC::ScalarField>,
        C: RegisterCommitments<OC::Affine>,
        AC: RegisterCommitments<OC::Affine>,
    {
        let t_subgroup = start_timer!(|| "subgroup checks");
        let mut points = proof.register_commitments.as_vec();
        points.extend(proof.additional_commitments.as_vec());
        points.push(proof.q_comm.to_affine());
        points.extend(proof.w_at_zeta_proof.points());
        points.extend(proof.r_at_zeta_omega_proof.points());
        let in_g1 = points.iter().all(OC::is_in_prime_subgroup);
        end_timer!(t_subgroup);
        in_g1
    }

    fn validate_evaluations<E, C, AC, P>(
        &self,
        proof: &Proof<OC::ScalarField, E, C, AC, S::C, S::Proof>,
        protocol: &P,
        challenges: &Challenges<OC::ScalarField>,
        fsrng: &mut TranscriptRng,
        evals_at_zeta: &LagrangeEvaluations<OC::ScalarField>,
    ) -> bool
    where
        E: RegisterEvaluations<OC::ScalarField>,
        C: RegisterCommitments<OC::Affine>,
        AC: RegisterCommitments<OC::Affine>,
        P: VerifierProtocol<IC, OC, S, C1 = C, C2 = AC>,
    {
        let t_pcs = start_timer!(|| "PCS verification");

        // Reconstruct the commitment to the linearization polynomial
        let t_r_comm = start_timer!(|| "linearization polynomial commitment");
        let r_comm = protocol
            .restore_commitment_to_linearization_polynomial(
                challenges.phi,
                evals_at_zeta.zeta_minus_omega_inv,
                &proof.register_commitments,
                &proof.additional_commitments,
            )
            .into_affine();
        end_timer!(t_r_comm);

        // Aggregate the commitments to be opened at ζ
        let t_aggregate_claims = start_timer!(|| "aggregate evaluation claims at zeta");
        let mut commitment_points = vec![
            self.pks_comm.pks_comm.0.to_affine(),
            self.pks_comm.pks_comm.1.to_affine(),
        ];
        commitment_points.extend(proof.register_commitments.as_vec());
        commitment_points.extend(proof.additional_commitments.as_vec());
        commitment_points.push(proof.q_comm.to_affine());

        let mut register_evals = proof.register_evaluations.as_vec();
        register_evals.push(proof.q_zeta);

        // Both lengths are fixed by the proof type, and `nus` was drawn for exactly that many.
        debug_assert_eq!(commitment_points.len(), challenges.nus.len());
        debug_assert_eq!(register_evals.len(), challenges.nus.len());

        let (w_comm_affine, w_at_zeta) =
            aggregate_claims_multiexp(commitment_points, register_evals, &challenges.nus);
        end_timer!(t_aggregate_claims);

        // Batch verify the two opening proofs
        let t_batch_opening = start_timer!(|| "batched PCS opening verification");

        // Convert affine points back to commitments
        let w_comm = S::C::from_affine(w_comm_affine);
        let r_comm_wrapped = S::C::from_affine(r_comm);

        // Prepare vectors for batch verification
        let commitments = vec![w_comm, r_comm_wrapped];
        let points = vec![challenges.zeta, evals_at_zeta.zeta_omega];
        let values = vec![w_at_zeta, proof.r_zeta_omega];
        let proofs = vec![
            proof.w_at_zeta_proof.clone(),
            proof.r_at_zeta_omega_proof.clone(),
        ];

        let verified = S::batch_verify(
            &self.verifier_key.prepare(),
            commitments,
            points,
            values,
            proofs,
            fsrng, // Use the transcript RNG for randomness
        )
        .is_ok();

        end_timer!(t_batch_opening);
        end_timer!(t_pcs);
        // Returned rather than asserted: proofs come from untrusted sources, so a failed
        // opening check is an invalid proof to reject, not a reason to abort the process.
        verified
    }

    fn restore_challenges<E, C, AC>(
        &self,
        public_input: &impl PublicInput<IC>,
        proof: &Proof<OC::ScalarField, E, C, AC, S::C, S::Proof>,
        batch_size: usize,
    ) -> (Challenges<OC::ScalarField>, TranscriptRng)
    where
        E: RegisterEvaluations<OC::ScalarField>,
        C: RegisterCommitments<OC::Affine>,
        AC: RegisterCommitments<OC::Affine>,
    {
        let mut transcript = self.preprocessed_transcript.clone();

        <Transcript as ApkTranscript<OC::ScalarField>>::append_public_input(
            &mut transcript,
            public_input,
        );
        <Transcript as ApkTranscript<OC::ScalarField>>::append_register_commitments(
            &mut transcript,
            &proof.register_commitments,
        );
        let r = <Transcript as ApkTranscript<OC::ScalarField>>::get_bitmask_aggregation_challenge(
            &mut transcript,
        );
        <Transcript as ApkTranscript<OC::ScalarField>>::append_2nd_round_register_commitments(
            &mut transcript,
            &proof.additional_commitments,
        );
        let phi =
            <Transcript as ApkTranscript<OC::ScalarField>>::get_constraints_aggregation_challenge(
                &mut transcript,
            );
        <Transcript as ApkTranscript<OC::ScalarField>>::append_quotient_commitment(
            &mut transcript,
            &proof.q_comm,
        );
        let zeta =
            <Transcript as ApkTranscript<OC::ScalarField>>::get_evaluation_point(&mut transcript);
        <Transcript as ApkTranscript<OC::ScalarField>>::append_evaluations(
            &mut transcript,
            &proof.register_evaluations,
            &proof.q_zeta,
            &proof.r_zeta_omega,
        );
        let nus =
            <Transcript as ApkTranscript<OC::ScalarField>>::get_opening_aggregation_challenges(
                &mut transcript,
                batch_size,
            );
        <Transcript as ApkTranscript<OC::ScalarField>>::append_opening_proofs(
            &mut transcript,
            &proof.w_at_zeta_proof,
            &proof.r_at_zeta_omega_proof,
        );

        (
            Challenges { r, phi, zeta, nus },
            fiat_shamir_rng(&mut transcript),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::setup::InsecureSetup;
    use crate::test_helpers::random_pks;
    use crate::{
        AccountablePublicInputOf, AccumulatorSeed, Apk, ApkConfig, ApkError, Bitmask,
        Bls12_377Config, Bls12_381Config, CommitmentExt, CountingPublicInputOf, PcsParamsOf,
        ScalarOf, VerifierOf,
    };
    use ark_ec::{CurveGroup, PrimeGroup};
    use ark_std::test_rng;
    use w3f_pcs::pcs::{PcsParams, PCS};

    fn check_rejects_malformed_inputs<C>()
    where
        C: ApkConfig,
        ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
        <C::Pcs as PCS<ScalarOf<C>>>::C:
            CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
        PcsParamsOf<C>: Clone,
        InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
    {
        let rng = &mut test_rng();
        let n = 10;
        let params = Apk::<C>::setup(&mut InsecureSetup::new(rng), n).unwrap();
        let (keyset, commitment) = Apk::<C>::commit_keyset(&params, random_pks(n, rng)).unwrap();
        let domain_size = commitment.domain_size as usize;
        assert_eq!(commitment.keyset_size, n as u64);

        // The commitment: at least one key, and fewer keys than the domain has rows.
        for bad in [0, commitment.domain_size, commitment.domain_size + 1] {
            let mut comm = commitment.clone();
            comm.keyset_size = bad;
            assert!(matches!(
                VerifierOf::<C>::try_new(params.raw_vk(), comm, C::transcript()),
                Err(ApkError::InvalidKeysetCommitment { .. })
            ));
        }

        let mut bits = vec![false; n];
        bits[2] = true;
        let (proof, public_input) = Apk::<C>::prove(
            &params,
            keyset.clone(),
            &commitment,
            Bitmask::from_bits(&bits),
        )
        .unwrap();
        assert_eq!(
            Apk::<C>::verify(&params, commitment.clone(), &public_input, &proof),
            Ok(true)
        );

        // Bit `domain_size + 2` evaluates at the same point as bit 2, so without the length
        // check this bitmask would name a signer who does not exist.
        let mut aliased = vec![false; 2 * domain_size];
        aliased[domain_size + 2] = true;
        for bits in [aliased, vec![true; n - 1], vec![true; n + 1]] {
            let forged = AccountablePublicInputOf::<C> {
                apk: public_input.apk,
                bitmask: Bitmask::from_bits(&bits),
            };
            assert_eq!(
                Apk::<C>::verify(&params, commitment.clone(), &forged, &proof),
                Err(ApkError::BitmaskLengthMismatch {
                    bitmask: bits.len(),
                    keyset: n
                })
            );
        }

        // An apk with h + apk = 0 has no affine form to evaluate the constraints at.
        let forged = AccountablePublicInputOf::<C> {
            apk: (-C::InnerCurve::accumulator_seed().into()).into_affine(),
            bitmask: Bitmask::from_bits(&bits),
        };
        assert_eq!(
            Apk::<C>::verify(&params, commitment.clone(), &forged, &proof),
            Err(ApkError::InvalidPublicInput(
                "apk is the negated accumulator seed"
            ))
        );

        // Counting: between 1 and the keyset size.
        let (proof, public_input) =
            Apk::<C>::prove_counting(&params, keyset, &commitment, Bitmask::from_bits(&bits))
                .unwrap();
        for count in [0, n + 1] {
            let forged = CountingPublicInputOf::<C> {
                apk: public_input.apk,
                count,
            };
            assert_eq!(
                Apk::<C>::verify_counting(&params, commitment.clone(), &forged, &proof),
                Err(ApkError::CountOutOfRange {
                    count,
                    keyset_size: n
                })
            );
        }
        assert_eq!(
            Apk::<C>::verify_counting(&params, commitment, &public_input, &proof),
            Ok(true)
        );
    }

    #[test]
    fn rejects_malformed_inputs() {
        check_rejects_malformed_inputs::<Bls12_377Config>();
        check_rejects_malformed_inputs::<Bls12_381Config>();
    }

    /// `bad_outer` is on the outer curve but outside G1; `off_curve_apk` is not on the inner
    /// curve at all.
    fn check_rejects_points_outside_g1<C>(
        bad_outer: <C::OuterCurve as CurveGroup>::Affine,
        off_curve_apk: <C::InnerCurve as CurveGroup>::Affine,
    ) where
        C: ApkConfig,
        C::Pcs: PCS<ScalarOf<C>, Proof = <C::OuterCurve as CurveGroup>::Affine>,
        ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
        <C::Pcs as PCS<ScalarOf<C>>>::C:
            CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
        PcsParamsOf<C>: Clone,
        InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
    {
        use crate::PrimeSubgroup;
        use ark_ec::AffineRepr;
        type Comm<C> = <<C as ApkConfig>::Pcs as PCS<ScalarOf<C>>>::C;

        assert!(!C::OuterCurve::is_in_prime_subgroup(&bad_outer));
        let rng = &mut test_rng();
        let n = 10;
        let params = Apk::<C>::setup(&mut InsecureSetup::new(rng), n).unwrap();
        let (keyset, commitment) = Apk::<C>::commit_keyset(&params, random_pks(n, rng)).unwrap();
        let verify = |pi: &AccountablePublicInputOf<C>, proof: &_| {
            Apk::<C>::verify(&params, commitment.clone(), pi, proof)
        };

        // A keyset commitment outside G1 is refused before any proof is looked at.
        for i in 0..2 {
            let mut comm = commitment.clone();
            let bad = Comm::<C>::from_affine(bad_outer);
            if i == 0 {
                comm.pks_comm.0 = bad;
            } else {
                comm.pks_comm.1 = bad;
            }
            assert!(matches!(
                VerifierOf::<C>::try_new(params.raw_vk(), comm, C::transcript()),
                Err(ApkError::KeysetCommitmentNotInG1)
            ));
        }

        // Every curve point of a proof, in turn, swapped for one outside G1.
        let bits = Bitmask::from_bits(&[true, false].repeat(n / 2));
        let (mut proof, public_input) =
            Apk::<C>::prove(&params, keyset.clone(), &commitment, bits.clone()).unwrap();
        assert_eq!(verify(&public_input, &proof), Ok(true));
        macro_rules! swapped_out_fails {
            ($field:expr, $bad:expr) => {{
                let original = core::mem::replace(&mut $field, $bad);
                assert_eq!(verify(&public_input, &proof), Ok(false));
                $field = original;
            }};
        }
        swapped_out_fails!(proof.register_commitments.0, bad_outer);
        swapped_out_fails!(proof.register_commitments.1, bad_outer);
        swapped_out_fails!(proof.q_comm, Comm::<C>::from_affine(bad_outer));
        swapped_out_fails!(proof.w_at_zeta_proof, bad_outer);
        swapped_out_fails!(proof.r_at_zeta_omega_proof, bad_outer);
        assert_eq!(verify(&public_input, &proof), Ok(true));

        // An apk off by a point outside G1, or not on the curve at all, is not checked for
        // either, and needs not be: the proof does not verify for it.
        let h = C::InnerCurve::accumulator_seed();
        for apk in [(public_input.apk + h).into_affine(), off_curve_apk] {
            let forged = AccountablePublicInputOf::<C> {
                apk,
                bitmask: bits.clone(),
            };
            assert_eq!(verify(&forged, &proof), Ok(false));
        }

        // The identity as apk: an empty signer set, which a valid proof exists for.
        let empty = AccountablePublicInputOf::<C> {
            apk: <C::InnerCurve as CurveGroup>::Affine::zero(),
            bitmask: Bitmask::from_bits(&vec![false; n]),
        };
        assert_eq!(
            verify(&empty, &proof),
            Err(ApkError::InvalidPublicInput("apk is the identity"))
        );
        let (proof, public_input) =
            Apk::<C>::prove_counting(&params, keyset, &commitment, bits).unwrap();
        let forged = CountingPublicInputOf::<C> {
            apk: <C::InnerCurve as CurveGroup>::Affine::zero(),
            count: public_input.count,
        };
        assert_eq!(
            Apk::<C>::verify_counting(&params, commitment, &forged, &proof),
            Err(ApkError::InvalidPublicInput("apk is the identity"))
        );
    }

    #[test]
    fn rejects_points_outside_g1() {
        use crate::test_helpers::point_of_order;
        let rng = &mut test_rng();

        let t = point_of_order::<ark_bw6_761::g1::Config, _>(2, rng);
        check_rejects_points_outside_g1::<Bls12_377Config>(
            (ark_bw6_761::G1Projective::generator() + t).into_affine(),
            ark_bls12_377::G1Affine::new_unchecked(1u8.into(), 1u8.into()),
        );

        // The order-3 point (0, 1): what the naive BW6-767 check would have let through.
        let t = ark_bw6_767::G1Affine::new_unchecked(0u8.into(), 1u8.into());
        check_rejects_points_outside_g1::<Bls12_381Config>(
            (ark_bw6_767::G1Projective::generator() + t).into_affine(),
            ark_bls12_381::G1Affine::new_unchecked(1u8.into(), 1u8.into()),
        );
    }
}
