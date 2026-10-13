//! Commitment-scheme parameters (the SRS) for a validator set.
//!
//! Where the parameters come from is abstracted as [`PcsSetup`], and the sizing — which domain a
//! validator set implies, and so which degree the parameters must reach — is done once, here,
//! for every source.
//!
//! The only source so far is [`InsecureSetup`], which samples the KZG trapdoor `tau` locally.
//! Whoever knows `tau` can open any commitment to any value, i.e. forge every proof, so it is
//! compiled only for this crate's tests and under the `test-utils` feature. Production
//! parameters have to come from a trusted-setup ceremony, through a source still to be designed
//! that implements the same trait.

use ark_ff::PrimeField;
use w3f_pcs::pcs::{CommitterKey, PcsParams, PCS};

use crate::domain::{DomainError, DomainSet, FftDomain};

/// A source of commitment-scheme parameters for the scheme `S` over the field `F`.
pub trait PcsSetup<F: PrimeField, S: PCS<F>> {
    /// Why this source could not produce parameters.
    type Error: core::fmt::Debug;

    /// Parameters supporting polynomials up to `max_degree`.
    fn params(&mut self, max_degree: usize) -> Result<S::Params, Self::Error>;
}

/// Why parameters for a validator set could not be produced.
#[derive(Debug)]
pub enum SetupError<E> {
    /// No evaluation domain exists for the requested size.
    NoDomain(DomainError),
    /// The parameter source failed.
    Source(E),
}

impl<E: core::fmt::Display> core::fmt::Display for SetupError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            SetupError::NoDomain(e) => write!(f, "{}", e),
            SetupError::Source(e) => write!(f, "{}", e),
        }
    }
}

/// Parameters for a keyset of `keyset_size` public keys, from `source`.
pub fn params_for_keyset<F, S, D, P>(
    source: &mut P,
    keyset_size: usize,
) -> Result<S::Params, SetupError<P::Error>>
where
    F: PrimeField,
    S: PCS<F>,
    D: DomainSet<F>,
    P: PcsSetup<F, S>,
{
    // The additional slot is occupied by the affine addition accumulator's initial value.
    params_for_domain::<F, S, D, P>(source, keyset_size + 1)
}

/// Parameters sufficient for a domain of at least `min_domain_size`, from `source`.
///
/// The domain set, not the caller, decides the realized size: `for_min_size` rounds up to the
/// next size the field supports, and the SRS is sized against that. Whether such a domain exists
/// at all is the domain implementation's business — for a radix-2 field that is a two-adicity
/// question, for BW6-767 it is a question of which divisors of `q - 1` are reachable — so the
/// failure surfaces from `for_min_size` rather than from a two-adicity assertion here.
pub fn params_for_domain<F, S, D, P>(
    source: &mut P,
    min_domain_size: usize,
) -> Result<S::Params, SetupError<P::Error>>
where
    F: PrimeField,
    S: PCS<F>,
    D: DomainSet<F>,
    P: PcsSetup<F, S>,
{
    let domain_size = D::for_min_size(min_domain_size)
        .map_err(SetupError::NoDomain)?
        .base()
        .size();
    source
        .params(highest_degree_to_commit(domain_size))
        .map_err(SetupError::Source)
}

/// **Insecure, test-only.** A parameter source that samples the trapdoor locally.
///
/// Parameters from here are only as secret as the rng that made them: whoever holds it can forge
/// every proof. See the [module docs](self).
#[cfg(any(test, feature = "test-utils"))]
pub struct InsecureSetup {
    rng: ark_std::rand::rngs::StdRng,
}

#[cfg(any(test, feature = "test-utils"))]
impl InsecureSetup {
    /// Seeds the source's own rng from `rng`, so the caller keeps using theirs afterwards.
    pub fn new<R: ark_std::rand::RngCore + ?Sized>(rng: &mut R) -> Self {
        use ark_std::rand::SeedableRng;
        let mut seed = [0u8; 32];
        rng.fill_bytes(&mut seed);
        InsecureSetup {
            rng: ark_std::rand::rngs::StdRng::from_seed(seed),
        }
    }
}

/// KZG over any pairing, BW6-761 and BW6-767 alike.
///
/// Not a call to upstream's `KZG::setup`: that delegates to `URS::from_trapdoor`, which asserts
/// `n <= 2^TWO_ADICITY`, and BW6-767's scalar field has two-adicity 1. The assertion is a policy
/// guard, not a mathematical constraint — generating a URS is powers of tau and one batch
/// multiplication, and the KZG operations this crate uses are multi-scalar multiplications and
/// pairings; only w3f-pcs's optional Lagrangian committer key needs a radix-2 domain, and this
/// crate never asks for one. So this reimplements `from_trapdoor` minus the assertion.
#[cfg(any(test, feature = "test-utils"))]
impl<E: ark_ec::pairing::Pairing> PcsSetup<E::ScalarField, w3f_pcs::pcs::kzg::KZG<E>>
    for InsecureSetup
{
    type Error = core::convert::Infallible;

    fn params(&mut self, max_degree: usize) -> Result<w3f_pcs::pcs::kzg::urs::URS<E>, Self::Error> {
        use ark_ec::ScalarMul;
        use ark_ff::One;

        let n1 = max_degree + 1;
        let n2 = 2;
        let (tau, g1, g2) = w3f_pcs::pcs::kzg::urs::URS::<E>::random_params(&mut self.rng);

        let mut powers_of_tau = Vec::with_capacity(n1.max(n2));
        let mut power = E::ScalarField::one();
        for _ in 0..n1.max(n2) {
            powers_of_tau.push(power);
            power *= tau;
        }

        Ok(w3f_pcs::pcs::kzg::urs::URS {
            powers_in_g1: g1.batch_mul(&powers_of_tau[..n1]),
            powers_in_g2: g2.batch_mul(&powers_of_tau[..n2]),
        })
    }
}

/// The highest polynomial degree the prover needs to commit to.
///
/// That is the quotient `q = aggregate_constraint_polynomial / vanishing_polynomial`. The
/// highest constraint degree is `4n - 3`, so `deg(q) = 3n - 3`.
pub(crate) fn highest_degree_to_commit(domain_size: usize) -> usize {
    3 * domain_size - 3
}

/// Verify that PCS parameters are sufficient for a given domain size
pub fn params_fit<S, F>(params: &S::Params, domain_size: usize) -> bool
where
    F: PrimeField,
    S: PCS<F>,
{
    highest_degree_to_commit(domain_size) <= params.ck().max_degree()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::Radix2DomainSet;
    use ark_bw6_761::{Fr, BW6_761};
    use ark_std::test_rng;
    use w3f_pcs::pcs::kzg::KZG;

    type TestKzg = KZG<BW6_761>;
    type TestDomain = Radix2DomainSet<Fr>;

    #[test]
    fn test_generate_for_domain() {
        let rng = &mut test_rng();
        let domain_size = 256;

        let params = params_for_domain::<Fr, TestKzg, TestDomain, _>(
            &mut InsecureSetup::new(rng),
            domain_size,
        )
        .unwrap();

        assert!(params_fit::<TestKzg, Fr>(&params, domain_size));
    }

    /// The SRS must be sized against the domain the prover will actually build, not the size
    /// that was asked for.
    #[test]
    fn test_generate_for_domain_rounds_up() {
        let rng = &mut test_rng();

        let params =
            params_for_domain::<Fr, TestKzg, TestDomain, _>(&mut InsecureSetup::new(rng), 200)
                .unwrap();

        // 200 is not a valid radix-2 size; the prover will get a domain of 256.
        assert!(params_fit::<TestKzg, Fr>(&params, 256));
    }

    #[test]
    fn test_generate_for_keyset() {
        let rng = &mut test_rng();
        let keyset_size = 100;

        let params = params_for_keyset::<Fr, TestKzg, TestDomain, _>(
            &mut InsecureSetup::new(rng),
            keyset_size,
        )
        .unwrap();

        // keyset_size + 1 (for the accumulator), rounded up to a power of two.
        let required_domain_size = (keyset_size + 1).next_power_of_two();
        assert!(params_fit::<TestKzg, Fr>(&params, required_domain_size));
    }

    /// The SRS is sized from the domain the prover will actually build, never from the
    /// caller's request. That matters most where a field's subgroup orders are sparse: on
    /// BW6-767 a keyset of 3243 lands on a domain of 11891, and an SRS cut to the request
    /// would be a third of the size the quotient needs.
    ///
    /// Demonstrated with a deliberately gappy table rather than the real one, so the property
    /// is pinned independently of how APK-381's sizes happen to be spaced today.
    #[test]
    fn srs_follows_the_realised_domain_not_the_request() {
        use crate::domain::{DomainSizes, DomainTriple, SmoothDomainSet};
        type Fr767 = ark_bw6_767::Fr;

        struct Gappy;
        impl DomainSizes<Fr767> for Gappy {
            fn triple_for(min_size: usize) -> Option<DomainTriple> {
                // Legal but sparse: every entry is still a subgroup order of this field.
                // `.iter().copied()` rather than `.into_iter()`: this crate is on edition
                // 2018, where an array's `into_iter` still yields references.
                [1usize, 3, 1081]
                    .iter()
                    .copied()
                    .find(|&n| n >= min_size)
                    .map(|n| DomainTriple::new(n, 2 * n, 6 * n))
            }
        }

        let realised = SmoothDomainSet::<Fr767, Gappy>::for_min_size(31)
            .unwrap()
            .base()
            .size();
        assert_eq!(realised, 1081, "31 rows round up across the gap");
        assert_eq!(highest_degree_to_commit(realised), 3 * 1081 - 3);
        // What sizing off the request would have given: 36x too small.
        assert!(highest_degree_to_commit(31) < highest_degree_to_commit(realised) / 30);
    }

    /// BW6-761's scalar field has two-adicity 46, so a domain of 2^50 does not exist. The
    /// failure is now a typed error from the domain rather than an assertion in this module.
    #[test]
    fn test_insufficient_adicity() {
        use crate::domain::{DomainError, DomainSet};
        assert_eq!(
            TestDomain::for_min_size(2usize.pow(50)).err(),
            Some(DomainError::TooLarge {
                requested: 2usize.pow(50)
            })
        );
    }
}
