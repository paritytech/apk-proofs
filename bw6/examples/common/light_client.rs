//! A light-client simulation, written once against [`ApkConfig`] and run by both examples.
//!
//! This is the crate's primary intended use case: communication-efficient light clients for
//! blockchains.
//!
//! A blockchain is modelled as a set of validators responsible for signing chain events. The
//! validator set changes in periods called 'eras'; within an era only a fraction of the set is
//! assumed malicious or unresponsive.
//!
//! A light client is a resource-constrained client — a mobile app, or better, a smart contract —
//! interested in some chain events but unable to follow the chain itself. It relies on an
//! untrusted helper node that supplies cryptographic proofs of the events it asks for.
//!
//! Such a proof could be a collection of signatures on the event from the relevant validator set,
//! but that would require the client to know every validator's public key, which is inefficient.
//! Knowing the aggregate public key of the set does not help either, since some individual
//! signatures may be missing.
//!
//! This crate provides succinct proofs that a public key is the aggregate public key of a
//! *subset* of the validator set. The whole set is identified by a short commitment, the subset
//! by a bitmask. That is effectively an accountable subset signature whose public key is the
//! commitment.
//!
//! The fundamental event is the validator set change. Given a recent commitment, the client can
//! process proofs of any other event — block finality, say — in the same way.
//!
//! The client's state starts at a commitment `C0` to the genesis validator set. On an era change
//! the helper provides `(C1, asig0, apk0, b0, p0)`: the new commitment, an aggregate signature by
//! a subset of the previous era's validators on `C1`, that subset's aggregate public key, the
//! bitmask identifying it, and a proof that `apk0` really is the aggregate key of the subset `b0`
//! of the set committed to by `C0`. The client verifies the proof, verifies the signature against
//! `apk0`, checks the bitmask makes quorum, and advances to `C1`.
//!
//! Nothing below names a curve. `run::<Bls12_377Config>` and `run::<Bls12_381Config>` are the
//! same code; the config decides the curves, the commitment scheme and the evaluation domains.

#![allow(dead_code)]

use std::cell::RefCell;

use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup};
use ark_serialize::CanonicalSerialize;
use ark_std::{end_timer, start_timer};
use rand::Rng;
use w3f_pcs::pcs::PCS;

use apk_proofs::bls::{PublicKey, SecretKey, Signature};
use apk_proofs::config::{
    AccountablePublicInputOf, Apk, ApkConfig, KeysetCommitmentOf, PcsParamsOf, ScalarOf,
    SimpleProofOf,
};
use apk_proofs::{hash_to_curve, ApkError, Bitmask, CommitmentExt};

type Inner<C> = <C as ApkConfig>::InnerPairing;
type G2Of<C> = <Inner<C> as Pairing>::G2;

pub struct Validator<C: ApkConfig>(SecretKey<Inner<C>>);

impl<C: ApkConfig> Clone for Validator<C> {
    fn clone(&self) -> Self {
        Validator(SecretKey::from_scalar(*self.0.as_ref()))
    }
}

pub struct Approval<C: ApkConfig> {
    comm: KeysetCommitmentOf<C>,
    sig: Signature<Inner<C>>,
    pk: PublicKey<Inner<C>>,
}

// Computing the commitment to a new validator set is expensive, and in a real deployment every
// validator does it in parallel. That is modelled by computing it once and sharing it.
thread_local! {
    static CACHE: RefCell<Option<Box<dyn std::any::Any>>> = RefCell::new(None);
}

fn new_era() {
    CACHE.with(|cell| {
        cell.replace(None);
    });
}

fn shared_keyset_commitment<C, F>(f: F) -> KeysetCommitmentOf<C>
where
    C: ApkConfig,
    KeysetCommitmentOf<C>: Clone + 'static,
    F: FnOnce() -> KeysetCommitmentOf<C>,
{
    CACHE.with(|cell| {
        let mut cell = cell.borrow_mut();
        if let Some(cached) = cell
            .as_ref()
            .and_then(|b| b.downcast_ref::<KeysetCommitmentOf<C>>())
        {
            return cached.clone();
        }
        let fresh = f();
        *cell = Some(Box::new(fresh.clone()));
        fresh
    })
}

impl<C> Validator<C>
where
    C: ApkConfig,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    apk_proofs::setup::InsecureSetup: apk_proofs::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
    KeysetCommitmentOf<C>: 'static,
{
    fn new<R: Rng>(rng: &mut R) -> Self {
        Self(SecretKey::new(rng))
    }

    fn public_key(&self) -> PublicKey<Inner<C>> {
        PublicKey::from(&self.0)
    }

    fn approve(&self, new_validator_set: &ValidatorSet<C>, params: &PcsParamsOf<C>) -> Approval<C> {
        let new_validator_set_commitment = shared_keyset_commitment::<C, _>(|| {
            Apk::<C>::commit_keyset(params, new_validator_set.raw_public_keys())
                .expect("validators only approve well-formed validator sets")
                .1
        });
        let message = hash_commitment::<C>(&new_validator_set_commitment);
        Approval {
            comm: new_validator_set_commitment,
            sig: self.0.sign(&message),
            pk: self.public_key(),
        }
    }
}

pub struct ValidatorSet<C: ApkConfig> {
    validators: Vec<Validator<C>>,
    quorum: usize,
}

impl<C: ApkConfig> Clone for ValidatorSet<C> {
    fn clone(&self) -> Self {
        Self {
            validators: self.validators.to_vec(),
            quorum: self.quorum,
        }
    }
}

impl<C> ValidatorSet<C>
where
    C: ApkConfig,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    apk_proofs::setup::InsecureSetup: apk_proofs::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
    KeysetCommitmentOf<C>: 'static,
{
    fn new<R: Rng>(size: usize, quorum: usize, rng: &mut R) -> Self {
        let validators = (0..size).map(|_| Validator::new(rng)).collect();
        Self { validators, quorum }
    }

    fn public_keys(&self) -> Vec<PublicKey<Inner<C>>> {
        self.validators.iter().map(|v| v.public_key()).collect()
    }

    fn raw_public_keys(&self) -> Vec<C::InnerCurve> {
        self.public_keys()
            .iter()
            .map(|pk| pk.as_group().clone())
            .collect()
    }

    fn size(&self) -> usize {
        self.validators.len()
    }

    fn rotate<R: Rng>(
        &self,
        params: &PcsParamsOf<C>,
        rng: &mut R,
    ) -> (ValidatorSet<C>, Vec<Approval<C>>) {
        new_era();
        let new_validator_set = ValidatorSet::new(self.size(), self.quorum, rng);

        let t_approval = start_timer!(|| {
            format!(
            "Each (honest) validator computes the commitment to the new validator set of size {} and signs the commitment",
            new_validator_set.size()
        )
        });

        let approvals = self
            .validators
            .iter()
            .filter(|_| rng.gen_bool(0.8))
            .map(|v| v.approve(&new_validator_set, params))
            .collect();

        end_timer!(t_approval);
        println!();

        (new_validator_set, approvals)
    }
}

fn hash_commitment<C: ApkConfig>(commitment: &KeysetCommitmentOf<C>) -> G2Of<C> {
    let mut buf = vec![0u8; commitment.compressed_size()];
    commitment.serialize_compressed(&mut buf[..]).unwrap();
    hash_to_curve(&buf)
}

pub struct LightClient<C: ApkConfig> {
    params: PcsParamsOf<C>,
    current_validator_set_commitment: KeysetCommitmentOf<C>,
    quorum: usize,
}

impl<C> LightClient<C>
where
    C: ApkConfig,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    apk_proofs::setup::InsecureSetup: apk_proofs::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
    KeysetCommitmentOf<C>: 'static,
{
    fn init(
        params: PcsParamsOf<C>,
        genesis_keyset_commitment: KeysetCommitmentOf<C>,
        quorum: usize,
    ) -> Self {
        Self {
            params,
            current_validator_set_commitment: genesis_keyset_commitment,
            quorum,
        }
    }

    fn verify_aggregates(
        &mut self,
        public_input: AccountablePublicInputOf<C>,
        proof: &SimpleProofOf<C>,
        aggregate_signature: &Signature<Inner<C>>,
        new_validator_set_commitment: KeysetCommitmentOf<C>,
    ) {
        let n_signers = public_input.bitmask.count_ones();
        let t_verification = start_timer!(|| format!(
            "Light client verifies light client proof for {} signers",
            n_signers
        ));

        let t_apk = start_timer!(|| "apk proof verification");
        assert!(Apk::<C>::verify(
            &self.params,
            self.current_validator_set_commitment.clone(),
            &public_input,
            proof,
        )
        .expect("the commitment and public input must be well-formed"));
        end_timer!(t_apk);

        let t_bls = start_timer!(|| "aggregate BLS signature verification");
        let aggregate_public_key = PublicKey::<Inner<C>>::from_group(public_input.apk.into_group());
        let message = hash_commitment::<C>(&new_validator_set_commitment);
        assert!(aggregate_public_key.verify(aggregate_signature, &message));
        end_timer!(t_bls);

        assert!(
            n_signers >= self.quorum,
            "{} signers don't make the quorum of {}",
            n_signers,
            self.quorum
        );

        self.current_validator_set_commitment = new_validator_set_commitment;

        end_timer!(t_verification);
    }
}

pub struct TrustlessHelper<C: ApkConfig> {
    params: PcsParamsOf<C>,
    current_validator_set: ValidatorSet<C>,
    current_validator_set_commitment: KeysetCommitmentOf<C>,
}

impl<C> TrustlessHelper<C>
where
    C: ApkConfig,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    apk_proofs::setup::InsecureSetup: apk_proofs::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
    KeysetCommitmentOf<C>: 'static,
{
    fn new(
        genesis_validator_set: ValidatorSet<C>,
        genesis_validator_set_commitment: KeysetCommitmentOf<C>,
        params: PcsParamsOf<C>,
    ) -> Self {
        Self {
            params,
            current_validator_set: genesis_validator_set,
            current_validator_set_commitment: genesis_validator_set_commitment,
        }
    }

    fn aggregate_approvals(
        &mut self,
        new_validator_set: ValidatorSet<C>,
        approvals: Vec<Approval<C>>,
    ) -> Result<
        (
            AccountablePublicInputOf<C>,
            SimpleProofOf<C>,
            Signature<Inner<C>>,
            KeysetCommitmentOf<C>,
        ),
        ApkError,
    > {
        let t_approval = start_timer!(|| {
            format!(
            "Helper aggregates {} individual signatures on the same commitment and generates accountable light client proof",
            approvals.len()
        )
        });

        let new_validator_set_commitment =
            approvals.first().ok_or(ApkError::NoSigners)?.comm.clone();
        // Compared as group elements rather than through a HashSet: `PublicKey<E>`'s derived
        // `Hash` would demand `E: Hash`, which no pairing implements.
        let actual_signers: Vec<C::InnerCurve> =
            approvals.iter().map(|a| *a.pk.as_group()).collect();
        let actual_signers_bitmask = self
            .current_validator_set
            .raw_public_keys()
            .iter()
            .map(|pk| actual_signers.contains(pk))
            .collect::<Vec<_>>();

        // The keyset is rebuilt from the current validator set rather than cached, because the
        // proof must be against exactly the set the light client has committed to.
        let (keyset, _) =
            Apk::<C>::commit_keyset(&self.params, self.current_validator_set.raw_public_keys())?;
        let (proof, public_input) = Apk::<C>::prove(
            &self.params,
            keyset,
            &self.current_validator_set_commitment,
            Bitmask::from_bits(&actual_signers_bitmask),
        )?;

        let signatures = approvals.iter().map(|a| &a.sig);
        let aggregate_signature = Signature::aggregate(signatures);

        self.current_validator_set = new_validator_set;
        self.current_validator_set_commitment = new_validator_set_commitment.clone();

        end_timer!(t_approval);
        println!();

        Ok((
            public_input,
            proof,
            aggregate_signature,
            new_validator_set_commitment,
        ))
    }
}

/// Runs `n_eras` validator-set rotations for a set of `validator_set_size` validators.
///
/// `validator_set_size` is a validator count, not a domain size and not a logarithm: the
/// configuration decides which evaluation domain that implies, and on APK-381 it is neither a
/// power of two nor anything the caller could usefully guess.
pub fn run<C, R>(validator_set_size: usize, n_eras: usize, rng: &mut R)
where
    C: ApkConfig,
    R: Rng,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    apk_proofs::setup::InsecureSetup: apk_proofs::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
    KeysetCommitmentOf<C>: 'static,
{
    let domain_size = Apk::<C>::domain_size(validator_set_size)
        .expect("no evaluation domain large enough for this validator set");
    println!(
        "Running {} with {} validators for {} eras.",
        C::NAME,
        validator_set_size,
        n_eras
    );
    println!(
        "The configuration picked an evaluation domain of {} points for {} validators.\n",
        domain_size, validator_set_size
    );

    let t_setup = start_timer!(|| format!(
        "Generating PCS params to support {} signers",
        validator_set_size
    ));
    // Insecure source (trapdoor sampled locally): tests only.
    let params = Apk::<C>::setup(
        &mut apk_proofs::setup::InsecureSetup::new(rng),
        validator_set_size,
    )
    .unwrap();
    end_timer!(t_setup);

    let quorum = ((validator_set_size * 2) / 3) + 1;
    println!(
        "\nGenesis: validator set size = {}, quorum = {}\n",
        validator_set_size, quorum
    );

    let t_genesis = start_timer!(|| format!(
        "Computing commitment to the set of initial {} validators",
        validator_set_size
    ));
    let genesis_validator_set = ValidatorSet::<C>::new(validator_set_size, quorum, rng);
    let (_, genesis_validator_set_commitment) =
        Apk::<C>::commit_keyset(&params, genesis_validator_set.raw_public_keys())
            .expect("genesis validator set is well-formed");
    end_timer!(t_genesis);

    let mut helper = TrustlessHelper::<C>::new(
        genesis_validator_set.clone(),
        genesis_validator_set_commitment.clone(),
        params.clone(),
    );
    let mut light_client = LightClient::<C>::init(params, genesis_validator_set_commitment, quorum);

    let mut current_validator_set = genesis_validator_set;

    for era in 1..=n_eras {
        println!("\nEra {}\n", era);
        let (new_validator_set, approvals) =
            current_validator_set.rotate(&light_client.params, rng);

        // A helper runs this once per block and must outlive a bad one: log the failure and
        // wait for the next round rather than taking the process down.
        let (public_input, proof, aggregate_signature, new_validator_set_commitment) =
            match helper.aggregate_approvals(new_validator_set.clone(), approvals) {
                Ok(out) => out,
                Err(e) => {
                    eprintln!("era {}: skipping, could not produce a proof: {}", era, e);
                    continue;
                }
            };

        light_client.verify_aggregates(
            public_input,
            &proof,
            &aggregate_signature,
            new_validator_set_commitment,
        );

        current_validator_set = new_validator_set;
    }
}

/// Parses `VALIDATORS [N_ERAS]` from the command line.
///
/// `example` is the example target, used only in the hint printed when `VALIDATORS` is missing;
/// which configuration runs is fixed by the calling `main`.
pub fn parse_args(example: &str, default_validators: usize, default_eras: usize) -> (usize, usize) {
    let mut args = std::env::args();
    args.next();

    let validators: usize = match args.next() {
        Some(a) => a.parse().expect("invalid VALIDATORS"),
        None => {
            println!(
                "VALIDATORS not given, using {}. To choose the validator set size and the \
                 number of eras, pass them after `--`, e.g.\n  \
                 cargo run --release --features \"parallel print-trace\" --example {} -- 1023 2",
                default_validators, example
            );
            default_validators
        }
    };
    let eras: usize = args
        .next()
        .map(|a| a.parse().expect("invalid N_ERAS"))
        .unwrap_or(default_eras);

    (validators, eras)
}
