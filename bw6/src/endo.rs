//! Endomorphism-based G1 membership test for BW6 curves, used for BW6-761.
//!
//! BW6-761 G1 has an efficient endomorphism `φ: (x, y) ↦ (ωx, y)`, `ω` a primitive cube root of
//! unity in the base field, acting on G1 as multiplication by an eigenvalue `λ`. Membership in
//! G1 then reduces to a few multiplications by the (64-bit) BLS12-377 seed `u` instead of one by
//! the 377-bit group order `r`.

use ark_ec::{
    bw6::{BW6Config, G1Projective},
    AdditiveGroup,
};
use ark_ff::{BitIteratorBE, Zero};
use std::ops::AddAssign;

// See https://github.com/celo-org/zexe/blob/master/algebra/src/bw6_761/curves/g1.rs#L37-L71
// and also https://github.com/celo-org/zexe/blob/master/scripts/glv_lattice_basis/src/lib.rs

/// Double-and-add by the non-negative integer whose little-endian limbs are `u`.
fn mul_by_u<C: BW6Config>(p: &G1Projective<C>, u: &[u64]) -> G1Projective<C> {
    let mut res = G1Projective::<C>::zero();
    for i in BitIteratorBE::without_leading_zeros(u) {
        res.double_in_place();
        if i {
            res.add_assign(p)
        }
    }
    res
}

/// `φ(p) = (ωx, y)`, computed on projective coordinates.
///
/// `omega` is trusted to be the cube root of unity whose eigenvalue [`subgroup_check`]'s
/// formula assumes; for BW6-761 that pairing is pinned by the tests below.
fn glv_endomorphism_proj<C: BW6Config>(p: &G1Projective<C>, omega: C::Fp) -> G1Projective<C> {
    G1Projective::<C>::new_unchecked(p.x * omega, p.y, p.z)
}

/// Whether `p` is in G1, given that it is on the curve: checks
/// `[u + 1]p + φ([u^3 - u^2 + 1]p) = 0`.
///
/// `(u + 1, u^3 - u^2 + 1)` is a short vector of the lattice of pairs `(a, b)` with
/// `a + λb ≡ 0 (mod r)`, so the combination vanishes on G1. That it vanishes nowhere else is a
/// property of BW6-761's cofactor, checked against every small cofactor order by
/// `instances::bls12_377_bw6_761::prime_subgroup`.
/// The same test as gnark-crypto's `G1Jac.IsInSubGroup` for BW6-761,
/// <https://github.com/Consensys/gnark-crypto/blob/fffb11080f9fcf9be636edd955860680e37e90f3/ecc/bw6-761/g1.go#L637-L660>.
/// Also refer to [ePrint 2020/351](https://eprint.iacr.org/2020/351.pdf), section 3.1
///
/// `omega` and `u` must be BW6-761's [`OMEGA`](crate::instances::bls12_377_bw6_761::OMEGA) and
/// [`U`](crate::instances::bls12_377_bw6_761::U): the formula is specific to that curve.
pub fn subgroup_check<C: BW6Config>(p: &G1Projective<C>, omega: C::Fp, u: &[u64]) -> bool {
    let up = mul_by_u::<C>(p, u);
    let u2p = mul_by_u::<C>(&up, u);
    let u3p = mul_by_u::<C>(&u2p, u);

    (up + p + glv_endomorphism_proj::<C>(&(u3p - u2p + p), omega)).is_zero()
}

#[cfg(test)]
mod tests {
    use crate::instances::bls12_377_bw6_761::{LAMBDA, OMEGA, U};
    use ark_bw6_761::{Config, Fq, G1Affine};
    use ark_ec::{AffineRepr, CurveGroup};
    use ark_ff::{Field, One};
    use ark_std::{test_rng, UniformRand};

    use super::*;

    fn glv_endomorphism_in_place(p: &mut G1Affine) {
        let x = &mut p.x;
        *x *= &OMEGA;
    }

    fn glv_endomorphism(p: &G1Affine) -> G1Affine {
        let mut p = p.clone();
        glv_endomorphism_in_place(&mut p);
        p
    }

    #[test]
    pub fn test_omega() {
        assert!(OMEGA.pow([3]).is_one());
    }

    #[test]
    pub fn test_endo() {
        let rng = &mut test_rng();

        let p1 = ark_bw6_761::G1Projective::rand(rng).into_affine();
        let mut p2 = p1.clone();

        assert_eq!(glv_endomorphism(&p1), p1 * LAMBDA);
        glv_endomorphism_in_place(&mut p2);
        assert_eq!(p2, p1 * LAMBDA);
    }

    #[test]
    pub fn test_endo_proj() {
        let rng = &mut test_rng();

        let p = ark_bw6_761::G1Projective::rand(rng);

        assert_eq!(glv_endomorphism_proj::<Config>(&p, OMEGA), p * LAMBDA);
    }

    #[test]
    pub fn test_subgroup_check() {
        let rng = &mut test_rng();

        let p = ark_bw6_761::G1Projective::rand(rng);

        assert!(subgroup_check::<Config>(&p, OMEGA, U));

        let point_not_in_g1 = loop {
            let x = Fq::rand(rng);
            let p = G1Affine::get_point_from_x_unchecked(x, false);
            if p.is_some() && !p.unwrap().is_in_correct_subgroup_assuming_on_curve() {
                break p.unwrap();
            }
        };

        assert!(point_not_in_g1.is_on_curve());
        assert!(!point_not_in_g1.is_in_correct_subgroup_assuming_on_curve());
        assert!(!subgroup_check::<Config>(
            &point_not_in_g1.into_group(),
            OMEGA,
            U
        ));
    }
}
