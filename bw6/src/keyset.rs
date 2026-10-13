//! Signer sets and their commitment.

use crate::domain::{DomainSet, FftDomain};
use crate::domains::{Domains, Evals};
use crate::hash_to_curve;
use crate::{ApkError, PrimeSubgroup, PublicKeyFault};
use ark_ec::AffineRepr;
use ark_ec::CurveGroup;
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, DenseUVPolynomial};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use std::marker::PhantomData;
use w3f_pcs::pcs::Commitment;
use w3f_pcs::pcs::{CommitterKey, PCS};

/// What the verifier knows a signer set by: commitments to the key coordinates, and the two
/// sizes it needs to rebuild the domain and bound the bitmask.
///
/// Formally: let `pks` be
/// a vector of public keys with `commit(pks) == KeysetCommitment::pks_comm`, and let
/// `domain_size := KeysetCommitment::domain_size`,
/// `keyset_size := KeysetCommitment::keyset_size`. Then the verifier needs to trust that:
/// 1. `pks` is well-formed:
///    - `pks.len() == domain_size`;
///    - `pks[i]` lies in the inner curve's G1 for `i = 0,...,domain_size-2`;
///    - for the real keys `pks[i]`, `i = 0,...,keyset_size-1`, there exist proofs of
///      possession, and for the padding `pks[i]`, `i = keyset_size,...,domain_size-2`, the
///      discrete log is not known.
///
///    `pks[domain_size-1]` is not part of the relation (not constrained) and could be anything.
///    [`Keyset::new`] pads every row from `keyset_size` on, the last included, with the same
///    point, `hash_to_curve(b"apk-proofs")`.
/// 2. the coordinate vectors of `pks` are interpolated over the order-`domain_size` subgroup
///    of `F*`, i.e. the domain `DomainSet::<F>::base_for_exact_size(domain_size)` returns.
///
/// In light client protocols the commitment is to the upcoming validator set, signed by the
/// current validator set. An honest validator checks the proofs of possession, interpolates with
/// the right padding over the right domain, computes the commitment using the right parameters,
/// and then signs it. The verifier checks the signatures and trusts that the properties hold
/// under an honest-supermajority(2/3 honest validators) assumption. As every honest validator computes the same
/// commitment, the verifier needs to check only the aggregate signature.
///
/// The commitment type is generic over PCS implementations. To extract the underlying curve
/// point, go through [`crate::CommitmentExt::to_affine`] rather than reaching into the concrete
/// type: w3f-pcs's KZG, for instance, wraps the point as `WrappedAffine(pub C::Affine)`, but
/// nothing here depends on that.
#[derive(Clone, Default, Debug, PartialEq, Eq, CanonicalSerialize, CanonicalDeserialize)]
pub struct KeysetCommitment<F, C>
where
    F: PrimeField,
    C: Commitment<F>,
{
    /// Per-coordinate commitments to public key polynomials
    pub pks_comm: (C, C),
    /// Size of the domain used to interpolate the vectors above.
    pub domain_size: u64,
    /// Number of real keys; the rest of the domain, but for its last row, is padding. The
    /// verifier holds a bitmask to exactly this many bits: a longer one would put bits on the
    /// padding rows, or wrap past the domain onto the real keys, since row `domain_size + i`
    /// evaluates at the same point as row `i`.
    pub keyset_size: u64,
    _m: PhantomData<F>,
}

/// A signer set as the prover holds it: the keys, their padded interpolations over the base
/// domain, and the domains themselves.
#[derive(Clone)]
pub struct Keyset<IC, OC, D>
where
    IC: CurveGroup,
    OC: CurveGroup,
    // OC::ScalarField: From<IC::BaseField>,
    D: DomainSet<OC::ScalarField>,
{
    // Actual public keys, no padding.
    pub(crate) pks: Vec<IC>,
    // Interpolations of the coordinate vectors of the public key vector WITH padding.
    pub(crate) pks_polys: [DensePolynomial<OC::ScalarField>; 2],
    // The domains used to compute the interpolations above, and to expand them.
    pub(crate) domains: D,
    // Polynomials above, evaluated over the large domain (at least 4n - 3 points).
    // Filled in by `amplify`, which `Prover::new` calls; used to populate the AIR execution trace.
    pub(crate) pks_evals_x4: Option<[Evals<OC::ScalarField>; 2]>,
}

impl<IC, OC, D> Keyset<IC, OC, D>
where
    IC: CurveGroup,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField>,
    D: DomainSet<OC::ScalarField>,
{
    /// Pads and interpolates a signer set.
    ///
    /// The keys typically come from the chain, so they are checked: each must be a non-identity
    /// point of the prime-order subgroup G1. The check costs about a scalar multiplication per
    /// key, paid once per validator set rather than once per proof.
    pub fn new(pks: Vec<IC>) -> Result<Self, ApkError>
    where
        IC: PrimeSubgroup,
    {
        // One row more than keys: the affine-addition accumulator starts from its seed.
        // Done before key validation: this is a lookup, validation is per key.
        let min_domain_size = pks.len() + 1;
        let domains = D::for_min_size(min_domain_size)?;
        let domain = domains.base();

        let mut padded_pks = pks.clone();
        // a point with unknown discrete log
        let padding_pk = hash_to_curve::<IC>(b"apk-proofs");
        padded_pks.resize(domain.size(), padding_pk);

        // convert into affine coordinates to commit
        let affine_pks = IC::normalize_batch(&padded_pks);
        for (index, pk) in affine_pks[..pks.len()].iter().enumerate() {
            let fault = if pk.is_zero() {
                PublicKeyFault::Identity
            } else if !IC::is_in_prime_subgroup(pk) {
                PublicKeyFault::NotInSubgroup
            } else {
                continue;
            };
            return Err(ApkError::InvalidPublicKey { index, fault });
        }

        let mut pks_x = Vec::with_capacity(affine_pks.len());
        let mut pks_y = Vec::with_capacity(affine_pks.len());
        for affine_point in &affine_pks {
            #[allow(clippy::expect_used, reason = "invariant argued in the message")]
            let (x, y) = affine_point
                .xy()
                .expect("invariant: real keys were checked above, padding is a random G1 element");
            pks_x.push((x).into());
            pks_y.push((y).into());
        }
        let pks_x_poly = DensePolynomial::from_coefficients_vec(domain.interpolate(&pks_x));
        let pks_y_poly = DensePolynomial::from_coefficients_vec(domain.interpolate(&pks_y));
        Ok(Self {
            pks,
            domains,
            pks_polys: [pks_x_poly, pks_y_poly],
            pks_evals_x4: None,
        })
    }

    /// The public keys, without padding.
    pub fn pks(&self) -> &[IC] {
        &self.pks
    }

    /// Actual number of signers, without including the padding.
    pub fn size(&self) -> usize {
        self.pks.len()
    }

    /// The trace domain the public keys were interpolated over.
    pub fn domain(&self) -> &D::Domain {
        self.domains.base()
    }

    /// Evaluates the key polynomials over the large domain.
    pub fn amplify(&mut self) {
        let domains = Domains::<OC::ScalarField, D>::from_set(self.domains.clone());
        let pks_evals_x4 = self
            .pks_polys
            .clone()
            .map(|z| domains.amplify_polynomial(&z));
        self.pks_evals_x4 = Some(pks_evals_x4);
    }

    /// Commits to the two key polynomials, recording the domain size and key count alongside.
    ///
    /// Fails if `committer_key` does not reach degree `n - 1`.
    pub fn commit<S>(
        &self,
        committer_key: &S::CK,
    ) -> Result<KeysetCommitment<OC::ScalarField, S::C>, ApkError>
    where
        S: PCS<OC::ScalarField>,
    {
        let domain_size = self.domain().size();
        // The keyset polynomials have degree `n - 1`.
        if domain_size > committer_key.max_degree() + 1 {
            return Err(ApkError::SrsTooSmall {
                domain_size,
                required_degree: domain_size - 1,
                available_degree: committer_key.max_degree(),
            });
        }
        let pks_x_comm = S::commit(committer_key, &self.pks_polys[0])
            .map_err(|_| ApkError::Pcs("commit to keyset x-coordinates"))?;
        let pks_y_comm = S::commit(committer_key, &self.pks_polys[1])
            .map_err(|_| ApkError::Pcs("commit to keyset y-coordinates"))?;
        Ok(KeysetCommitment {
            pks_comm: (pks_x_comm, pks_y_comm),
            domain_size: domain_size as u64,
            keyset_size: self.size() as u64,
            _m: PhantomData::default(),
        })
    }

    /// The sum of the keys `bitmask` selects: the aggregate public key a proof will claim.
    pub fn aggregate(&self, bitmask: &[bool]) -> Result<IC, ApkError> {
        if bitmask.len() != self.size() {
            return Err(ApkError::BitmaskLengthMismatch {
                bitmask: bitmask.len(),
                keyset: self.size(),
            });
        }
        Ok(bitmask
            .iter()
            .zip(self.pks.iter())
            .filter(|(b, _p)| **b)
            .map(|(_b, p)| p)
            .sum())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::setup::InsecureSetup;
    use crate::test_helpers::random_pks;
    use crate::{
        Apk, ApkConfig, Bls12_377Config, Bls12_381Config, KeysetOf, PcsParamsOf, ScalarOf,
    };
    use ark_ec::short_weierstrass::SWCurveConfig;
    use ark_std::{test_rng, Zero};
    use w3f_pcs::pcs::{PcsParams, PCS};

    fn check_rejects_bad_keys<C, P>()
    where
        C: ApkConfig,
        C::InnerCurve: From<ark_ec::short_weierstrass::Affine<P>>,
        P: SWCurveConfig,
        ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField>,
    {
        let rng = &mut test_rng();

        let mut pks = random_pks::<_, C::InnerCurve>(10, rng);
        pks[3] = C::InnerCurve::zero();
        assert_eq!(
            KeysetOf::<C>::new(pks).err(),
            Some(ApkError::InvalidPublicKey {
                index: 3,
                fault: PublicKeyFault::Identity
            })
        );

        // On the curve, but of order 3: the accumulator seed itself is such a point.
        let mut pks = random_pks::<_, C::InnerCurve>(10, rng);
        pks[5] = crate::point_in_g1_complement::<P>().into();
        assert_eq!(
            KeysetOf::<C>::new(pks).err(),
            Some(ApkError::InvalidPublicKey {
                index: 5,
                fault: PublicKeyFault::NotInSubgroup
            })
        );
    }

    #[test]
    fn rejects_identity_and_non_subgroup_keys() {
        check_rejects_bad_keys::<Bls12_377Config, ark_bls12_377::g1::Config>();
        check_rejects_bad_keys::<Bls12_381Config, ark_bls12_381::g1::Config>();
    }

    fn check_commit_reports_small_srs<C>()
    where
        C: ApkConfig,
        ScalarOf<C>: From<<C::InnerCurve as CurveGroup>::BaseField> + ark_ff::FftField,
        <C::Pcs as PCS<ScalarOf<C>>>::C: crate::CommitmentExt<ScalarOf<C>, Affine = <C::OuterCurve as CurveGroup>::Affine>
            + Clone,
        PcsParamsOf<C>: Clone,
        InsecureSetup: crate::setup::PcsSetup<ScalarOf<C>, C::Pcs>,
    {
        let rng = &mut test_rng();
        // Parameters for 2 keys, a keyset of 100.
        let params = Apk::<C>::setup(&mut InsecureSetup::new(rng), 2).unwrap();
        let keyset = KeysetOf::<C>::new(random_pks(100, rng)).unwrap();
        match keyset.commit::<C::Pcs>(&params.ck()) {
            Err(ApkError::SrsTooSmall { domain_size, .. }) => {
                assert_eq!(domain_size, keyset.domain().size())
            }
            other => panic!("expected SrsTooSmall, got {:?}", other.err()),
        }
    }

    #[test]
    fn commit_reports_small_srs() {
        check_commit_reports_small_srs::<Bls12_377Config>();
        check_commit_reports_small_srs::<Bls12_381Config>();
    }

    #[test]
    fn aggregate_rejects_bitmask_of_wrong_length() {
        let rng = &mut test_rng();
        let keyset = KeysetOf::<Bls12_381Config>::new(random_pks(10, rng)).unwrap();
        assert_eq!(
            keyset.aggregate(&[true; 9]).err(),
            Some(ApkError::BitmaskLengthMismatch {
                bitmask: 9,
                keyset: 10
            })
        );
    }
}
