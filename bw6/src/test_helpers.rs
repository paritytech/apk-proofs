//! Shared test bodies, written once against [`ApkConfig`] and run on both configurations.
//!
//! This module is the working demonstration that the two curves share a protocol: apart from the
//! proof sizes (BW6-767 group elements are one byte longer) and the absence of `packed` on
//! APK-381, every function below is called twice with nothing but the config type changed.

use ark_ec::{
    short_weierstrass::{Affine, Projective, SWCurveConfig},
    AffineRepr, CurveConfig, CurveGroup, PrimeGroup,
};
use ark_ff::{FftField, One, PrimeField, Zero};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use ark_std::{end_timer, start_timer, test_rng};
use ark_std::{rand::Rng, UniformRand};
use w3f_pcs::pcs::{PcsParams, PCS};

use crate::config::{
    AccountablePublicInputOf, Apk, ApkConfig, Bls12_377Config, Bls12_381Config, CountingProofOf,
    CountingPublicInputOf, KeysetCommitmentOf, KeysetOf, PackedProofOf, PcsParamsOf, ProverOf,
    ScalarOf, SimpleProofOf, VerifierOf,
};
use crate::domain::{DomainSet, FftDomain};
use crate::{Bitmask, CommitmentExt, PublicInput};

/// A uniformly random point of `E(Fq)`, not only of G1.
pub fn random_curve_point<P: SWCurveConfig, R: Rng>(rng: &mut R) -> Affine<P> {
    loop {
        let x = P::BaseField::rand(rng);
        if let Some(p) = Affine::<P>::get_point_from_x_unchecked(x, false) {
            return p;
        }
    }
}
/// Divides little-endian limbs by a small `d`, returning the quotient and the remainder.
pub fn div_small(limbs: &[u64], d: u64) -> (Vec<u64>, u64) {
    let mut q = vec![0u64; limbs.len()];
    let mut rem: u128 = 0;
    for i in (0..limbs.len()).rev() {
        let cur = (rem << 64) | limbs[i] as u128;
        q[i] = (cur / d as u128) as u64;
        rem = cur % d as u128;
    }
    (q, rem as u64)
}

/// A point of order exactly `l`, for a prime `l` dividing the cofactor.
///
/// Clears the G1 component with `[r]`, then every factor of the cofactor but `l`, then
/// multiplies by `l` until one more step would reach the identity.
pub fn point_of_order<P: SWCurveConfig, R: Rng>(l: u64, rng: &mut R) -> Projective<P> {
    let mut h_prime = P::COFACTOR.to_vec();
    loop {
        let (q, rem) = div_small(&h_prime, l);
        if rem != 0 {
            break;
        }
        h_prime = q;
    }
    assert_ne!(
        h_prime,
        P::COFACTOR.to_vec(),
        "{} does not divide the cofactor",
        l
    );
    for _ in 0..16 {
        let p = random_curve_point::<P, _>(rng).into_group();
        let mut t = p
            .mul_bigint(<P as CurveConfig>::ScalarField::MODULUS)
            .mul_bigint(&h_prime);
        if t.is_zero() {
            continue;
        }
        while !t.mul_bigint([l]).is_zero() {
            t = t.mul_bigint([l]);
        }
        return t;
    }
    panic!("no point of order {} found", l);
}
/// Checks `test` against G1, the identity, points of each small order `l` dividing the
/// cofactor, those points shifted by a G1 point, and random points of `E(Fq)`.
pub fn check_exact<P: SWCurveConfig>(test: impl Fn(&Projective<P>) -> bool, small_orders: &[u64]) {
    let rng = &mut test_rng();
    let g = Projective::<P>::generator();
    assert!(test(&g));
    assert!(test(&Projective::<P>::zero()));
    for _ in 0..4 {
        assert!(test(&Projective::<P>::rand(rng)));
    }
    for &l in small_orders {
        let t = point_of_order::<P, _>(l, rng);
        assert!(!t.is_zero() && t.mul_bigint([l]).is_zero());
        assert!(!test(&t), "a point of order {} passes", l);
        assert!(
            !test(&(t + g)),
            "G1 shifted by a point of order {} passes",
            l
        );
    }
    for _ in 0..8 {
        let p = random_curve_point::<P, _>(rng);
        assert!(!p.is_in_correct_subgroup_assuming_on_curve());
        assert!(!test(&p.into_group()));
    }
}

pub(crate) fn _random_bits<R: Rng>(n: usize, density: f64, rng: &mut R) -> Vec<bool> {
    (0..n).map(|_| rng.gen_bool(density)).collect()
}

pub(crate) fn _random_bitmask<R: Rng, C: CurveGroup>(n: usize, rng: &mut R) -> Vec<C::ScalarField> {
    _random_bits(n, 2.0 / 3.0, rng)
        .into_iter()
        .map(|b| {
            if b {
                C::ScalarField::one()
            } else {
                C::ScalarField::zero()
            }
        })
        .collect()
}

pub fn random_pks<R: Rng, C: CurveGroup>(n: usize, rng: &mut R) -> Vec<C> {
    (0..n).map(|_| C::rand(rng)).collect()
}

// The `where` clause on each generic body below is the same four lines every time: the config,
// the two arithmetic facts the protocol relies on (the outer scalar field is the inner base
// field, and it supports FFTs), and cloneable PCS parameters. It is repeated rather than aliased
// because trait aliases are not stable.

fn keyset_and_commitment<C>(
    params: &PcsParamsOf<C>,
    keyset_size: usize,
) -> (KeysetOf<C>, KeysetCommitmentOf<C>)
where
    C: ApkConfig,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    crate::setup::InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
{
    let rng = &mut test_rng();
    let pks = random_pks::<_, C::InnerCurve>(keyset_size, rng);
    let t = start_timer!(|| "signer set commitment");
    let out = Apk::<C>::commit_keyset(params, pks).unwrap();
    end_timer!(t);
    out
}

/// Proves and verifies one proof of whatever scheme `prove`/`verify` name, and checks that it
/// survives a serialisation round trip at exactly `proof_size` bytes.
fn _test_prove_verify<C, ProofT, PI, P, V>(
    params: PcsParamsOf<C>,
    prove: P,
    verify: V,
    keyset_size: usize,
    proof_size: usize,
) where
    C: ApkConfig,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    crate::setup::InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
    ProofT: CanonicalSerialize + CanonicalDeserialize,
    PI: PublicInput<C::InnerCurve>,
    P: Fn(ProverOf<C>, Bitmask) -> (ProofT, PI),
    V: Fn(&VerifierOf<C>, &PI, &ProofT) -> bool,
{
    let rng = &mut test_rng();
    let (keyset, pks_comm) = keyset_and_commitment::<C>(&params, keyset_size);

    let t_prover_new = start_timer!(|| "prover precomputation");
    let prover = ProverOf::<C>::new(keyset, &pks_comm, params.clone(), C::transcript()).unwrap();
    end_timer!(t_prover_new);

    let verifier = VerifierOf::<C>::try_new(params.raw_vk(), pks_comm, C::transcript()).unwrap();

    let bits = (0..keyset_size)
        .map(|_| rng.gen_bool(2.0 / 3.0))
        .collect::<Vec<_>>();
    let b = Bitmask::from_bits(&bits);

    let prove_ = start_timer!(|| "prove");
    let (proof, public_input) = prove(prover, b.clone());
    end_timer!(prove_);

    let mut serialized_proof = vec![0; proof.compressed_size()];
    proof
        .serialize_compressed(&mut serialized_proof[..])
        .unwrap();
    let deserialized_proof = ProofT::deserialize_compressed(&serialized_proof[..]).unwrap();

    assert_eq!(proof.compressed_size(), proof_size);

    let verify_ = start_timer!(|| "verify");
    let valid = verify(&verifier, &public_input, &deserialized_proof);
    end_timer!(verify_);

    assert!(valid);
}

// ---------------------------------------------------------------------------------------------
// Proof element sizes
//
// Field elements match at 48 bytes across both curves (BW6-767's scalar field is 381 bits,
// BW6-761's 377). Compressed G1 points do not: BW6-767's base field is 767 bits, which fills 96
// bytes to within one spare bit, leaving no room for the two flags arkworks needs for the
// infinity marker and the y-sign, so the encoding spills into a 97th byte. BW6-761's 761-bit
// base field leaves seven spare bits, so its compressed points (the x-coordinate plus the two
// flags) fit in 96.
// ---------------------------------------------------------------------------------------------

/// The compressed size of a proof with `commitments` group elements of `group_bytes` each and
/// `field_elements` field elements of 48 bytes each (both scalar fields fit in 384 bits).
///
/// arkworks writes a field element with no flags. Its size is ceil(modulus_bits / 8):
/// ceil(377 / 8) = 48 bytes, with 7 spare bits
/// ceil(381 / 8) = 48 bytes, with 3 spare bits
/// Both fit in 384 bits = 48 bytes, hence the hard-coded * 48 for both setups.
fn proof_size(commitments: usize, field_elements: usize, group_bytes: usize) -> usize {
    commitments * group_bytes + field_elements * 48
}
// BIT_SIZE = 2, Group bytes size = ceil((field_bits + BIT_SIZE) / 8)
const GROUP_BYTES_761: usize = 96; // x-coordinate, flags in its spare bits
const GROUP_BYTES_767: usize = 97;

// ---------------------------------------------------------------------------------------------
// Positive tests, one body per scheme
// ---------------------------------------------------------------------------------------------

fn test_simple_scheme_for<C>(keyset_size: usize, group_bytes: usize)
where
    C: ApkConfig,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    crate::setup::InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
{
    _test_prove_verify::<C, SimpleProofOf<C>, AccountablePublicInputOf<C>, _, _>(
        // Insecure source (trapdoor sampled locally): tests only.
        Apk::<C>::setup(
            &mut crate::setup::InsecureSetup::new(&mut test_rng()),
            keyset_size,
        )
        .unwrap(),
        |prover, bitmask| prover.prove_simple(bitmask).unwrap(),
        |verifier, public_input, proof| verifier.verify_simple(public_input, proof).unwrap(),
        keyset_size,
        proof_size(5, 6, group_bytes),
    );
}

fn test_counting_scheme_for<C>(keyset_size: usize, group_bytes: usize)
where
    C: ApkConfig,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    crate::setup::InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
{
    _test_prove_verify::<C, CountingProofOf<C>, CountingPublicInputOf<C>, _, _>(
        // Insecure source (trapdoor sampled locally): tests only.
        Apk::<C>::setup(
            &mut crate::setup::InsecureSetup::new(&mut test_rng()),
            keyset_size,
        )
        .unwrap(),
        |prover, bitmask| prover.prove_counting(bitmask).unwrap(),
        |verifier, public_input, proof| verifier.verify_counting(public_input, proof).unwrap(),
        keyset_size,
        proof_size(7, 8, group_bytes),
    );
}

pub fn test_simple_scheme_377(keyset_size: usize) {
    test_simple_scheme_for::<Bls12_377Config>(keyset_size, GROUP_BYTES_761);
}

pub fn test_counting_scheme_377(keyset_size: usize) {
    test_counting_scheme_for::<Bls12_377Config>(keyset_size, GROUP_BYTES_761);
}

/// `packed` exists only on APK-377, so it is not written against a generic config: the bound
/// that would let it be is exactly the one [`Bls12_381Config`] does not satisfy.
pub fn test_packed_scheme_377(keyset_size: usize) {
    type C = Bls12_377Config;
    _test_prove_verify::<C, PackedProofOf<C>, AccountablePublicInputOf<C>, _, _>(
        // Insecure source (trapdoor sampled locally): tests only.
        Apk::<C>::setup(
            &mut crate::setup::InsecureSetup::new(&mut test_rng()),
            keyset_size,
        )
        .unwrap(),
        |prover, bitmask| prover.prove_packed(bitmask).unwrap(),
        |verifier, public_input, proof| verifier.verify_packed(public_input, proof).unwrap(),
        keyset_size,
        proof_size(8, 9, GROUP_BYTES_761),
    );
}

// ---------------------------------------------------------------------------------------------
// APK-381: the same bodies, a different config.
//
// `packed` is absent on purpose — it hardcodes 256-bit bitmask chunks and so needs `256 | n`,
// which no BW6-767 domain satisfies.
// ---------------------------------------------------------------------------------------------

pub fn test_simple_scheme_381(keyset_size: usize) {
    test_simple_scheme_for::<Bls12_381Config>(keyset_size, GROUP_BYTES_767);
}

pub fn test_counting_scheme_381(keyset_size: usize) {
    test_counting_scheme_for::<Bls12_381Config>(keyset_size, GROUP_BYTES_767);
}

/// `prove_packed` must not exist on APK-381. This function is never called; it is here so the
/// gate is visible, and so that deleting `SupportsPackedScheme` would be caught by the doc test
/// below rather than silently re-enabling a scheme that cannot work.
///
/// ```compile_fail
/// use apk_proofs::{Apk381, Bitmask};
/// fn f(params: &apk_proofs::config::PcsParamsOf<apk_proofs::Bls12_381Config>,
///      keyset: apk_proofs::config::KeysetOf<apk_proofs::Bls12_381Config>,
///      comm: &apk_proofs::config::KeysetCommitmentOf<apk_proofs::Bls12_381Config>,
///      b: Bitmask) {
///     let _ = Apk381::prove_packed(params, keyset, comm, b);
/// }
/// ```
///
/// The same call compiles for APK-377:
/// ```
/// use apk_proofs::{Apk377, Bitmask};
/// fn f(params: &apk_proofs::config::PcsParamsOf<apk_proofs::Bls12_377Config>,
///      keyset: apk_proofs::config::KeysetOf<apk_proofs::Bls12_377Config>,
///      comm: &apk_proofs::config::KeysetCommitmentOf<apk_proofs::Bls12_377Config>,
///      b: Bitmask) {
///     let _ = Apk377::prove_packed(params, keyset, comm, b);
/// }
/// ```
pub fn _packed_scheme_is_gated_to_apk_377() {}

// ---------------------------------------------------------------------------------------------
// Negative tests
//
// Everything above checks that an honest proof verifies. These check that dishonest ones do not,
// which is the property that actually matters and the one a refactor is most likely to break
// silently: a verifier that ignores part of its input still passes every positive test.
// ---------------------------------------------------------------------------------------------

/// Builds an honest prover/verifier pair, then asserts that verification fails once the claim
/// no longer matches the proof.
fn _test_rejects_tampering<C>(keyset_size: usize)
where
    C: ApkConfig,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    crate::setup::InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
{
    let rng = &mut test_rng();
    // Insecure source (trapdoor sampled locally): tests only.
    let params = Apk::<C>::setup(&mut crate::setup::InsecureSetup::new(rng), keyset_size).unwrap();
    let (keyset, pks_comm) = keyset_and_commitment::<C>(&params, keyset_size);

    let bits: Vec<bool> = (0..keyset_size).map(|_| rng.gen_bool(2.0 / 3.0)).collect();
    let (proof, public_input) =
        Apk::<C>::prove(&params, keyset, &pks_comm, Bitmask::from_bits(&bits)).unwrap();

    // Sanity: the honest claim verifies, so any failure below is caused by the tampering.
    assert!(
        Apk::<C>::verify(&params, pks_comm.clone(), &public_input, &proof).unwrap(),
        "the honest proof must verify"
    );

    // Flipping a bit changes the claimed signer set. The bitmask feeds the transcript, so this
    // moves every challenge, and it also breaks the aggregate-key relation.
    let mut flipped = bits.clone();
    flipped[0] = !flipped[0];
    let tampered =
        AccountablePublicInputOf::<C>::new(&public_input.apk, &Bitmask::from_bits(&flipped));
    assert!(
        !Apk::<C>::verify(&params, pks_comm.clone(), &tampered, &proof).unwrap(),
        "a proof must not verify against a different bitmask"
    );

    // Claiming a different aggregate key for the same bitmask.
    let wrong_apk = AccountablePublicInputOf::<C>::new(
        &(public_input.apk.into_group() + C::InnerCurve::generator()).into_affine(),
        &Bitmask::from_bits(&bits),
    );
    assert!(
        !Apk::<C>::verify(&params, pks_comm.clone(), &wrong_apk, &proof).unwrap(),
        "a proof must not verify against a different aggregate key"
    );

    // A verifier that disagrees with the prover about the domain must reject. This is what the
    // domain size and generator in the transcript are for: without them the two sides could
    // silently be working over different domains. The doubled size has to be a realisable one,
    // or `verify` would reject it as unconstructible before ever checking the proof.
    let mut wrong_domain_comm = pks_comm;
    wrong_domain_comm.domain_size =
        <C::Domains as DomainSet<ScalarOf<C>>>::for_min_size(2 * (keyset_size + 1))
            .expect("a larger domain must exist")
            .base()
            .size() as u64;
    assert!(
        !Apk::<C>::verify(&params, wrong_domain_comm, &public_input, &proof).unwrap(),
        "a proof must not verify against a verifier using a different domain"
    );
}

/// A keyset commitment arrives from the chain, so its domain size is untrusted. Neither a size
/// the field cannot realise exactly, nor one no field could realise, may panic the verifier.
fn _test_rejects_bad_domain_size<C>(keyset_size: usize, unrealisable: u64)
where
    C: ApkConfig,
    ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + FftField,
    <C::Pcs as PCS<ScalarOf<C>>>::C:
        CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
    PcsParamsOf<C>: Clone,
    crate::setup::InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
{
    let rng = &mut test_rng();
    // Insecure source (trapdoor sampled locally): tests only.
    let params = Apk::<C>::setup(&mut crate::setup::InsecureSetup::new(rng), keyset_size).unwrap();
    let (_keyset, pks_comm) = keyset_and_commitment::<C>(&params, keyset_size);

    let mut not_exact = pks_comm.clone();
    not_exact.domain_size = unrealisable;
    assert!(
        matches!(
            VerifierOf::<C>::try_new(params.raw_vk(), not_exact, C::transcript()),
            Err(crate::ApkError::Domain(crate::DomainError::NotExact { .. }))
        ),
        "a domain size the field cannot realise exactly must be rejected"
    );

    let mut too_large = pks_comm;
    too_large.domain_size = u64::MAX;
    assert!(
        matches!(
            VerifierOf::<C>::try_new(params.raw_vk(), too_large, C::transcript()),
            Err(crate::ApkError::Domain(crate::DomainError::TooLarge { .. }))
        ),
        "an unrealisable domain size must be rejected"
    );
}

pub fn test_rejects_tampering_377(keyset_size: usize) {
    _test_rejects_tampering::<Bls12_377Config>(keyset_size);
}

pub fn test_rejects_tampering_381(keyset_size: usize) {
    _test_rejects_tampering::<Bls12_381Config>(keyset_size);
}

pub fn test_rejects_bad_domain_size_377(keyset_size: usize) {
    // Not a power of two, so no radix-2 domain has exactly this size.
    _test_rejects_bad_domain_size::<Bls12_377Config>(keyset_size, 300);
}

pub fn test_rejects_bad_domain_size_381(keyset_size: usize) {
    // 254 is not a divisor of 3 * 11 * 23 * 47 * 10177, so it is not a table entry.
    _test_rejects_bad_domain_size::<Bls12_381Config>(keyset_size, 254);
}

/// An SRS generated for one keyset size, used with a larger one.
///
/// Nothing ties `Apk::setup`'s argument to the length of the vector handed to `commit_keyset`,
/// so this is a caller mistake the crate has to report well. It used to surface as a panic
/// on an opaque `Err(())` inside a commit closure, long after `commit_keyset` had cheerfully
/// succeeded: the keyset polynomials are degree `n - 1` and fit, while the quotient at `3n - 3`
/// does not.
pub fn test_undersized_srs_is_reported_at_prover_construction() {
    let rng = &mut test_rng();
    type C = Bls12_377Config;

    // 255 keys need a domain of 256; 300 need 512, and so an SRS twice as large.
    // Insecure source (trapdoor sampled locally): tests only.
    let params = Apk::<C>::setup(&mut crate::setup::InsecureSetup::new(rng), 255).unwrap();
    let pks = random_pks::<_, <C as ApkConfig>::InnerCurve>(300, rng);

    // The keyset polynomials (degree 511) still fit, so committing succeeds; the quotient
    // (degree 1533) does not, and the prover has to say so before doing any work.
    let (keyset, commitment) = Apk::<C>::commit_keyset(&params, pks).unwrap();
    let bitmask = Bitmask::from_bits(&vec![true; 300]);
    let err = Apk::<C>::prove(&params, keyset, &commitment, bitmask)
        .err()
        .expect("proving against an undersized SRS must not succeed");

    assert_eq!(
        err,
        crate::ApkError::SrsTooSmall {
            domain_size: 512,
            required_degree: 1533,
            available_degree: 765,
        }
    );
    assert!(err.to_string().contains("SRS too small"));
}

/// One function body, both curves — the same four calls, with only the config type changed.
pub fn test_config_driven_api() {
    fn roundtrip<C: ApkConfig>(keyset_size: usize) -> bool
    where
        ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + FftField,
        <C::Pcs as PCS<ScalarOf<C>>>::C:
            CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine> + Clone,
        PcsParamsOf<C>: Clone,
        crate::setup::InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
    {
        let rng = &mut test_rng();
        // Insecure source (trapdoor sampled locally): tests only.
        let params =
            Apk::<C>::setup(&mut crate::setup::InsecureSetup::new(rng), keyset_size).unwrap();
        let pks = random_pks::<_, C::InnerCurve>(keyset_size, rng);
        let (keyset, commitment) = Apk::<C>::commit_keyset(&params, pks).unwrap();
        let bitmask = Bitmask::from_bits(&vec![true; keyset_size]);
        let (proof, public_input) = Apk::<C>::prove(&params, keyset, &commitment, bitmask).unwrap();
        Apk::<C>::verify(&params, commitment, &public_input, &proof).unwrap()
    }

    assert!(roundtrip::<Bls12_377Config>(255));
    assert!(roundtrip::<Bls12_381Config>(252));
}
