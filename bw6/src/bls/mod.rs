//! Generic BLS signature implementations over a pairing-friendly curve

use ark_ec::pairing::Pairing;
use ark_ec::{AffineRepr, CurveGroup, PrimeGroup};
use ark_ff::{UniformRand, Zero};
use ark_serialize::*;
use rand::Rng;
use std::borrow::Borrow;
use std::ops::Neg;

#[derive(Clone, Debug)]
pub struct Signature<E: Pairing>(E::G2);

impl<E: Pairing> AsRef<E::G2> for Signature<E> {
    fn as_ref(&self) -> &E::G2 {
        &self.0
    }
}

impl<E: Pairing> Signature<E> {
    /// Create a signature from a G2 element
    pub fn from_group(g2: E::G2) -> Self {
        Signature(g2)
    }
    /// Aggregate multiple signatures into a single signature
    pub fn aggregate<S: Borrow<Self>>(signatures: impl IntoIterator<Item = S>) -> Self {
        Signature(signatures.into_iter().map(|s| s.borrow().0).sum::<E::G2>())
    }
}

/// BLS secret key
#[derive(Clone, Debug, CanonicalSerialize, CanonicalDeserialize)]
pub struct SecretKey<E: Pairing>(pub E::ScalarField);

impl<E: Pairing> AsRef<E::ScalarField> for SecretKey<E> {
    fn as_ref(&self) -> &E::ScalarField {
        &self.0
    }
}

impl<E: Pairing> SecretKey<E> {
    /// Generate a new random secret key
    pub fn new<R: Rng>(rng: &mut R) -> Self {
        SecretKey(E::ScalarField::rand(rng))
    }

    /// Create a secret key from a scalar
    pub fn from_scalar(scalar: E::ScalarField) -> Self {
        SecretKey(scalar)
    }

    /// Sign a message (represented as a G2 element)
    ///
    /// The signature is computed as: sig = sk * message
    ///
    /// # Note
    ///
    /// In practice, the message should be hashed to a G2 point using a
    /// hash-to-curve function.
    ///
    /// # Example
    ///
    /// ```
    /// use apk_proofs::bls::{PublicKey, SecretKey};
    /// use ark_bls12_381::{Bls12_381, G2Projective};
    ///
    /// let rng = &mut ark_std::test_rng();
    /// let sk = SecretKey::<Bls12_381>::new(rng);
    /// // A stand-in for hash-to-G2: `hash_to_curve` samples a point from a seeded RNG.
    /// let message: G2Projective = apk_proofs::hash_to_curve(b"Hello, world!");
    ///
    /// let signature = sk.sign(&message);
    /// assert!(PublicKey::from(&sk).verify(&signature, &message));
    /// ```
    pub fn sign(&self, message: &E::G2) -> Signature<E> {
        Signature(*message * self.as_ref())
    }
}

/// BLS public key
///
/// This is an element of G1.
#[derive(Clone, Debug, Eq, PartialEq, Hash, CanonicalSerialize, CanonicalDeserialize)]
pub struct PublicKey<E: Pairing>(E::G1);

impl<E: Pairing> From<&SecretKey<E>> for PublicKey<E> {
    fn from(sk: &SecretKey<E>) -> Self {
        PublicKey(E::G1::generator() * sk.0)
    }
}

impl<E: Pairing> PublicKey<E> {
    /// Create a public key from a G1 element
    pub fn from_group(g1: E::G1) -> Self {
        PublicKey(g1)
    }

    /// Returns the underlying G1 element
    pub fn as_group(&self) -> &E::G1 {
        &self.0
    }

    /// Aggregate multiple public keys into a single public key.
    ///
    /// This is a simple sum of G1 points.
    ///
    /// # Example
    ///
    /// ```
    /// use apk_proofs::bls::{PublicKey, SecretKey, Signature};
    /// use ark_bls12_381::{Bls12_381, G2Projective};
    ///
    /// let rng = &mut ark_std::test_rng();
    /// let sk1 = SecretKey::<Bls12_381>::new(rng);
    /// let sk2 = SecretKey::<Bls12_381>::new(rng);
    /// let message: G2Projective = apk_proofs::hash_to_curve(b"Hello, world!");
    ///
    /// let pk1 = PublicKey::from(&sk1);
    /// let pk2 = PublicKey::from(&sk2);
    /// let aggregate_pk = PublicKey::aggregate([pk1, pk2]);
    ///
    /// // The aggregate key verifies the aggregate of the signatures on the same message.
    /// let aggregate_sig = Signature::aggregate([sk1.sign(&message), sk2.sign(&message)]);
    /// assert!(aggregate_pk.verify(&aggregate_sig, &message));
    /// ```
    pub fn aggregate<P: Borrow<Self>>(public_keys: impl IntoIterator<Item = P>) -> Self {
        PublicKey(public_keys.into_iter().map(|p| p.borrow().0).sum::<E::G1>())
    }

    /// Verify a signature against the public key and the message
    ///
    /// Checks the pairing equation:
    /// ```text
    /// e(G1, sig) = e(pk, message)
    /// OR
    /// e(-G1, sig) * e(pk, message) = 1
    /// ```
    ///
    /// # Example
    ///
    /// ```
    /// use apk_proofs::bls::{PublicKey, SecretKey};
    /// use ark_bls12_381::{Bls12_381, G2Projective};
    ///
    /// let rng = &mut ark_std::test_rng();
    /// let sk = SecretKey::<Bls12_381>::new(rng);
    /// let message: G2Projective = apk_proofs::hash_to_curve(b"Hello, world!");
    /// let signature = sk.sign(&message);
    /// let pk = PublicKey::from(&sk);
    /// assert!(pk.verify(&signature, &message));
    ///
    /// // ...and not on any other message.
    /// let other: G2Projective = apk_proofs::hash_to_curve(b"Goodbye, world!");
    /// assert!(!pk.verify(&signature, &other));
    /// ```
    pub fn verify(&self, signature: &Signature<E>, message: &E::G2) -> bool {
        E::multi_pairing(
            [
                E::G1Affine::generator().into_group().neg().into_affine(),
                self.0.into_affine(),
            ],
            [signature.as_ref().into_affine(), message.into_affine()],
        )
        .is_zero()
    }
}

#[cfg(test)]
mod tests {
    use ark_bls12_377::{Bls12_377, G2Projective};
    use ark_std::test_rng;

    use super::*;

    #[test]
    fn test_apk() {
        let rng = &mut test_rng();
        let message = G2Projective::rand(rng);

        let sks = (0..10)
            .map(|_| SecretKey::<Bls12_377>::new(rng))
            .collect::<Vec<_>>();
        let pks = sks.iter().map(PublicKey::from).collect::<Vec<_>>();
        let sigs = sks.iter().map(|sk| sk.sign(&message)).collect::<Vec<_>>();
        pks.iter()
            .zip(sigs.iter())
            .for_each(|(pk, sig)| assert!(pk.verify(sig, &message)));

        let apk = PublicKey::aggregate(pks);
        let asig = Signature::aggregate(sigs);
        assert!(apk.verify(&asig, &message));
    }
}
