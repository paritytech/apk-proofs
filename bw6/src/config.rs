//! Curve configurations, and the API the crate is meant to be used through.
//!
//! An APK proof is parameterised by four things that are not independent: an inner BLS12 curve
//! carrying the signatures, an outer BW6 curve the proof is computed over, a polynomial
//! commitment scheme on the outer curve, and a set of evaluation domains for the outer curve's
//! scalar field. Only a handful of combinations are coherent — the outer curve's scalar field
//! has to be the inner curve's base field, and the domains have to match that field's structure.
//!
//! [`ApkConfig`] bundles them, so choosing a curve is naming one type rather than four, and
//! combinations that do not make sense cannot be spelled. Everything below the config —
//! [`crate::Keyset`], [`crate::Prover`], [`crate::Verifier`], the PIOP — carries a single type
//! parameter, and [`Apk`] hides even that:
//!
//! ```no_run
//! use apk_proofs::{Apk381, Bitmask};
//! # let pks: Vec<ark_bls12_381::G1Projective> = unimplemented!();
//! # let bitmask: Bitmask = unimplemented!();
//! # use apk_proofs::setup::InsecureSetup;
//! let rng = &mut ark_std::test_rng();
//! // Insecure source: samples the trapdoor locally, so tests only (`test-utils`). Production
//! // parameters must come from a trusted-setup ceremony; see `crate::setup`.
//! let setup = Apk381::setup(&mut InsecureSetup::new(rng), 1000).unwrap(); // validators, not log2(anything)
//! let (keyset, commitment) = Apk381::commit_keyset(&setup, pks)?;
//! let (proof, public_input) = Apk381::prove(&setup, keyset, &commitment, bitmask)?;
//! assert!(Apk381::verify(&setup, commitment, &public_input, &proof)?);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use ark_ec::pairing::Pairing;
use ark_ec::CurveGroup;
use w3f_pcs::pcs::PCS;

use crate::domain::DomainSet;

/// The scalar field of a configuration's outer curve, which is also the base field of its inner
/// curve. All polynomial arithmetic happens here.
pub type ScalarOf<C> = <<C as ApkConfig>::OuterCurve as ark_ec::PrimeGroup>::ScalarField;

/// A coherent choice of curves, commitment scheme and evaluation domains.
pub trait ApkConfig: 'static + Sized {
    /// Distinguishes configurations in the Fiat-Shamir transcript: it is the label
    /// [`ApkConfig::transcript`] opens with.
    const NAME: &'static str;

    /// Where the BLS public keys live.
    /// Must seed the accumulator outside G1; see [`crate::AccumulatorSeed`].
    type InnerCurve: crate::AccumulatorSeed;

    /// The pairing `InnerCurve` belongs to.
    ///
    /// The proof system itself only ever touches G1, but a caller signing with these keys needs
    /// G2 and the pairing — so naming it here keeps [`crate::bls`] reachable from a config
    /// rather than forcing callers back to a concrete curve. See `examples/`.
    type InnerPairing: Pairing<G1 = Self::InnerCurve>;

    /// The group the commitment scheme works in. Its scalar field is `InnerCurve`'s base field,
    /// which is what lets inner-curve coordinates be manipulated as native field elements.
    ///
    /// Both sides compute here: the prover commits and opens, and the verifier checks every
    /// curve point of a proof against [`crate::PrimeSubgroup`] before it reaches the pairing,
    /// and combines commitments linearly. The subgroup test is on the verifier's path, so it
    /// should be a fast one.
    type OuterCurve: CurveGroup + crate::PrimeSubgroup;

    /// Commitment scheme over the outer curve. Its opening proofs must list their curve points,
    /// for the verifier to check them; see [`crate::OpeningProofPoints`].
    type Pcs: PCS<
        ScalarOf<Self>,
        Proof: crate::OpeningProofPoints<<Self::OuterCurve as CurveGroup>::Affine>,
    >;

    /// The triple of evaluation domains — `n, 2n, 4n` on APK-377, `n, 2n, 6n` on APK-381 — and
    /// the shift that goes with them. This is the component that cannot be shared between
    /// configurations: see [`crate::domain`].
    type Domains: DomainSet<ScalarOf<Self>>;

    /// A fresh Fiat-Shamir transcript for this configuration.
    ///
    /// Prover and verifier must start from byte-identical transcripts, so both obtain it here
    /// rather than spelling out a label at each call site. Labelling with [`ApkConfig::NAME`]
    /// also means a proof made under one configuration never replays under another.
    fn transcript() -> merlin::Transcript {
        merlin::Transcript::new(Self::NAME.as_bytes())
    }
}

/// The prover for a configuration.
pub type ProverOf<C> = crate::Prover<
    <C as ApkConfig>::InnerCurve,
    <C as ApkConfig>::OuterCurve,
    <C as ApkConfig>::Pcs,
    <C as ApkConfig>::Domains,
>;

/// The verifier for a configuration.
pub type VerifierOf<C> = crate::Verifier<
    <C as ApkConfig>::InnerCurve,
    <C as ApkConfig>::OuterCurve,
    <C as ApkConfig>::Pcs,
    <C as ApkConfig>::Domains,
>;

/// The keyset type for a configuration.
pub type KeysetOf<C> = crate::Keyset<
    <C as ApkConfig>::InnerCurve,
    <C as ApkConfig>::OuterCurve,
    <C as ApkConfig>::Domains,
>;

/// The keyset commitment type for a configuration.
pub type KeysetCommitmentOf<C> =
    crate::KeysetCommitment<ScalarOf<C>, <<C as ApkConfig>::Pcs as PCS<ScalarOf<C>>>::C>;

/// The 'basic' scheme's proof type for a configuration.
pub type SimpleProofOf<C> = crate::SimpleProof<
    ScalarOf<C>,
    <<C as ApkConfig>::OuterCurve as CurveGroup>::Affine,
    <<C as ApkConfig>::Pcs as PCS<ScalarOf<C>>>::C,
    <<C as ApkConfig>::Pcs as PCS<ScalarOf<C>>>::Proof,
>;
/// The 'packed' scheme's proof type for a configuration.
pub type PackedProofOf<C> = crate::PackedProof<
    ScalarOf<C>,
    <<C as ApkConfig>::OuterCurve as CurveGroup>::Affine,
    <<C as ApkConfig>::Pcs as PCS<ScalarOf<C>>>::C,
    <<C as ApkConfig>::Pcs as PCS<ScalarOf<C>>>::Proof,
>;
/// The 'counting' scheme's proof type for a configuration.
pub type CountingProofOf<C> = crate::CountingProof<
    ScalarOf<C>,
    <<C as ApkConfig>::OuterCurve as CurveGroup>::Affine,
    <<C as ApkConfig>::Pcs as PCS<ScalarOf<C>>>::C,
    <<C as ApkConfig>::Pcs as PCS<ScalarOf<C>>>::Proof,
>;
/// The 'basic' and 'packed' schemes' public input for a configuration.
pub type AccountablePublicInputOf<C> = crate::AccountablePublicInput<<C as ApkConfig>::InnerCurve>;
/// The 'counting' scheme's public input for a configuration.
pub type CountingPublicInputOf<C> = crate::CountingPublicInput<<C as ApkConfig>::InnerCurve>;
/// The commitment-scheme parameters (for KZG, the SRS) for a configuration.
pub type PcsParamsOf<C> = <<C as ApkConfig>::Pcs as PCS<ScalarOf<C>>>::Params;

/// BLS12-377 signatures, proofs over BW6-761, KZG commitments, radix-2 domains.
///
/// BW6-761's scalar field has two-adicity 46, so power-of-two domains exist for every size the
/// prover could want and arkworks' radix-2 machinery applies unchanged. The domain triple is
/// `n, 2n, 4n`.
pub struct Bls12_377Config;

impl ApkConfig for Bls12_377Config {
    const NAME: &'static str = "apk-bls12-377-bw6-761";
    type InnerCurve = crate::instances::bls12_377_bw6_761::InnerCurve;
    type InnerPairing = crate::instances::bls12_377_bw6_761::InnerPairing;
    type OuterCurve = crate::instances::bls12_377_bw6_761::OuterCurve;
    type Pcs = crate::instances::bls12_377_bw6_761::kzg::Pcs;
    type Domains = crate::instances::bls12_377_bw6_761::Domains761;
}

/// BLS12-381 signatures, proofs over BW6-767, KZG commitments, mixed-radix domains.
///
/// BW6-767's scalar field has two-adicity **1**. No power-of-two domain beyond size 2 exists
/// there and no domain size is divisible by 4, so sizes are divisors of
/// `2 * 3^2 * 11 * 23 * 47 * 10177` and are transformed by Cooley-Tukey, with Rader for the
/// large factor. The domain triple is `n, 2n, 6n`; see
/// [`crate::instances::bls12_381_bw6_767::APK381_DOMAIN_SIZES`].
///
/// The `packed` scheme is unavailable here — it needs `256 | n`. Asking for it is a compile
/// error; see [`crate::SupportsPackedScheme`].
pub struct Bls12_381Config;

impl ApkConfig for Bls12_381Config {
    const NAME: &'static str = "apk-bls12-381-bw6-767";
    type InnerCurve = crate::instances::bls12_381_bw6_767::InnerCurve;
    type InnerPairing = crate::instances::bls12_381_bw6_767::InnerPairing;
    type OuterCurve = crate::instances::bls12_381_bw6_767::OuterCurve;
    type Pcs = crate::instances::bls12_381_bw6_767::kzg::Pcs;
    type Domains = crate::instances::bls12_381_bw6_767::Domains767;
}

/// The APK proof system for a configuration.
///
/// The entry point when you want the protocol rather than its parts: `Apk::<Bls12_381Config>`
/// picks the curves, the commitment scheme and the domains together. Use
/// [`ProverOf`]/[`VerifierOf`] directly when you need to reuse a prover across many proofs.
pub struct Apk<C: ApkConfig>(core::marker::PhantomData<C>);

/// APK proofs over BLS12-377 / BW6-761.
pub type Apk377 = Apk<Bls12_377Config>;

/// APK proofs over BLS12-381 / BW6-767.
pub type Apk381 = Apk<Bls12_381Config>;

impl<C: ApkConfig> Apk<C>
where
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        crate::CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
{
    /// Commitment-scheme parameters for a validator set of `keyset_size` keys, from `source`.
    ///
    /// `keyset_size` is the number of validators, not a logarithm and not a domain size: which
    /// domain that implies is the configuration's business, and on APK-381 it is not a power of
    /// two. Read it back with [`Apk::domain_size`].
    ///
    /// Whether the result is fit for production is the source's business, not this function's.
    /// See [`crate::setup`].
    pub fn setup<P: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>>(
        source: &mut P,
        keyset_size: usize,
    ) -> Result<PcsParamsOf<C>, crate::setup::SetupError<P::Error>> {
        crate::setup::params_for_keyset::<ScalarOf<C>, C::Pcs, C::Domains, P>(source, keyset_size)
    }

    /// The trace domain size this configuration will use for `keyset_size` validators.
    ///
    /// Always at least `keyset_size + 1` — the extra row holds the affine-addition accumulator's
    /// initial value — and rounded up to a size the field supports.
    pub fn domain_size(keyset_size: usize) -> Result<usize, crate::DomainError> {
        use crate::domain::FftDomain;
        Ok(
            <C::Domains as DomainSet<ScalarOf<C>>>::for_min_size(keyset_size + 1)?
                .base()
                .size(),
        )
    }

    /// Interpolates and commits to a signer set.
    ///
    /// The domain size follows from `pks.len()`; the commitment records it so the verifier
    /// reconstructs the same one.
    ///
    /// Fails if a key is the identity or outside G1, if no domain holds this many keys, or if
    /// `params` are too small to commit to them.
    pub fn commit_keyset(
        params: &PcsParamsOf<C>,
        pks: Vec<C::InnerCurve>,
    ) -> Result<(KeysetOf<C>, KeysetCommitmentOf<C>), crate::ApkError> {
        use w3f_pcs::pcs::PcsParams;
        let keyset = KeysetOf::<C>::new(pks)?;
        let commitment = keyset.commit::<C::Pcs>(&params.ck())?;
        Ok((keyset, commitment))
    }

    fn prover(
        params: &PcsParamsOf<C>,
        keyset: KeysetOf<C>,
        commitment: &KeysetCommitmentOf<C>,
    ) -> Result<ProverOf<C>, crate::ApkError> {
        ProverOf::<C>::new(keyset, commitment, params.clone(), C::transcript())
    }

    fn verifier(
        params: &PcsParamsOf<C>,
        commitment: KeysetCommitmentOf<C>,
    ) -> Result<VerifierOf<C>, crate::ApkError> {
        use w3f_pcs::pcs::PcsParams;
        VerifierOf::<C>::try_new(params.raw_vk(), commitment, C::transcript())
    }

    /// Proves that `bitmask` selects the signers whose aggregate key the returned public input
    /// names.
    pub fn prove(
        params: &PcsParamsOf<C>,
        keyset: KeysetOf<C>,
        commitment: &KeysetCommitmentOf<C>,
        bitmask: crate::Bitmask,
    ) -> Result<(SimpleProofOf<C>, AccountablePublicInputOf<C>), crate::ApkError> {
        Self::prover(params, keyset, commitment)?.prove_simple(bitmask)
    }

    /// Checks a proof against a claimed aggregate key and bitmask.
    ///
    /// `Ok(false)` for a proof that does not verify. `Err` for inputs rejected before the proof
    /// is looked at: a commitment whose sizes this configuration cannot use, or a public input
    /// that is malformed — a bitmask of the wrong length, or an unusable aggregate key.
    pub fn verify(
        params: &PcsParamsOf<C>,
        commitment: KeysetCommitmentOf<C>,
        public_input: &AccountablePublicInputOf<C>,
        proof: &SimpleProofOf<C>,
    ) -> Result<bool, crate::ApkError> {
        Self::verifier(params, commitment)?.verify_simple(public_input, proof)
    }

    /// Proves only how many validators signed, not which ones.
    pub fn prove_counting(
        params: &PcsParamsOf<C>,
        keyset: KeysetOf<C>,
        commitment: &KeysetCommitmentOf<C>,
        bitmask: crate::Bitmask,
    ) -> Result<(CountingProofOf<C>, CountingPublicInputOf<C>), crate::ApkError> {
        Self::prover(params, keyset, commitment)?.prove_counting(bitmask)
    }

    /// Checks a counting proof against a claimed aggregate key and signer count.
    ///
    /// As [`Apk::verify`]; a count of zero or above the keyset size is an `Err`.
    pub fn verify_counting(
        params: &PcsParamsOf<C>,
        commitment: KeysetCommitmentOf<C>,
        public_input: &CountingPublicInputOf<C>,
        proof: &CountingProofOf<C>,
    ) -> Result<bool, crate::ApkError> {
        Self::verifier(params, commitment)?.verify_counting(public_input, proof)
    }
}

/// The `packed` scheme, available only where the domains can supply sizes divisible by 256.
///
/// [`Bls12_381Config`] does not satisfy this and cannot be made to: `q - 1` for BW6-767 carries
/// a single factor of 2, so no domain size there is even a multiple of 4. Calling
/// `Apk::<Bls12_381Config>::prove_packed` is a compile error rather than a runtime panic.
impl<C: ApkConfig> Apk<C>
where
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        crate::CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    C::Domains: crate::SupportsPackedScheme,
{
    /// As [`Apk::prove`], producing a 'packed' proof.
    pub fn prove_packed(
        params: &PcsParamsOf<C>,
        keyset: KeysetOf<C>,
        commitment: &KeysetCommitmentOf<C>,
        bitmask: crate::Bitmask,
    ) -> Result<(PackedProofOf<C>, AccountablePublicInputOf<C>), crate::ApkError> {
        Self::prover(params, keyset, commitment)?.prove_packed(bitmask)
    }

    /// As [`Apk::verify`], for a 'packed' proof.
    pub fn verify_packed(
        params: &PcsParamsOf<C>,
        commitment: KeysetCommitmentOf<C>,
        public_input: &AccountablePublicInputOf<C>,
        proof: &PackedProofOf<C>,
    ) -> Result<bool, crate::ApkError> {
        Self::verifier(params, commitment)?.verify_packed(public_input, proof)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::FftDomain;
    use ark_std::{test_rng, UniformRand};

    /// The API the crate is meant to be used through: name a curve, get a proof.
    /// The same five lines run on either configuration.
    #[test]
    fn apk_facade_roundtrips_on_both_configurations() {
        fn roundtrip<C: ApkConfig>(n: usize) -> bool
        where
            ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
            <C::Pcs as PCS<ScalarOf<C>>>::C: crate::CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine>
                + Clone,
            PcsParamsOf<C>: Clone,
            crate::setup::InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
        {
            let rng = &mut test_rng();
            // Insecure source: fine for a test.
            let params = Apk::<C>::setup(&mut crate::setup::InsecureSetup::new(rng), n).unwrap();
            let pks: Vec<C::InnerCurve> = (0..n).map(|_| C::InnerCurve::rand(rng)).collect();

            let (keyset, commitment) = Apk::<C>::commit_keyset(&params, pks).unwrap();
            let bitmask = crate::Bitmask::from_bits(&vec![true; n]);
            let (proof, public_input) =
                Apk::<C>::prove(&params, keyset, &commitment, bitmask).unwrap();
            Apk::<C>::verify(&params, commitment, &public_input, &proof).unwrap()
        }

        assert!(roundtrip::<Bls12_377Config>(255));
        assert!(roundtrip::<Bls12_381Config>(252));
    }

    /// `setup` takes a validator count. Nothing in the API asks the caller for a log, or for a
    /// domain size: which domain a validator count implies is the configuration's business.
    #[test]
    fn setup_is_sized_by_validator_count_not_by_domain_size() {
        // 252 validators need 253 rows, which is 11 * 23 exactly.
        assert_eq!(Apk381::domain_size(252).unwrap(), 253);
        // 1000 rounds up to 23 * 47; 1500 to 3 * 11 * 47. Neither is a power of two.
        assert_eq!(Apk381::domain_size(1000).unwrap(), 1081);
        assert_eq!(Apk381::domain_size(1500).unwrap(), 1551);
        // The same call on APK-377 rounds up to a power of two instead.
        assert_eq!(Apk377::domain_size(1000).unwrap(), 1024);
        assert_eq!(Apk377::domain_size(1500).unwrap(), 2048);
    }

    /// Selecting a curve is naming one type, and the domain strategy follows from it rather than
    /// being chosen by the caller.
    #[test]
    fn a_config_determines_its_domain_strategy() {
        type D377 = <Bls12_377Config as ApkConfig>::Domains;
        type D381 = <Bls12_381Config as ApkConfig>::Domains;

        // Radix-2: powers of two, expanded 2x and 4x.
        let d377 = <D377 as DomainSet<_>>::for_min_size(200).unwrap();
        assert_eq!(d377.base().size(), 256);
        assert_eq!(d377.large().size(), 4 * 256);

        // Mixed-radix: divisors of 3 * 11 * 23 * 47 * 10177, expanded 2x and 6x. 256 is not one
        // of them, and no size here is even a multiple of 4.
        let d381 = <D381 as DomainSet<_>>::for_min_size(200).unwrap();
        assert_eq!(d381.base().size(), 253); // 11 * 23
        assert_eq!(d381.large().size(), 6 * 253);
        assert_ne!(
            d381.base().size() % 4,
            0,
            "no BW6-767 domain is a multiple of 4"
        );
    }

    #[test]
    fn configs_are_named_distinctly() {
        assert_ne!(Bls12_377Config::NAME, Bls12_381Config::NAME);
    }
}
