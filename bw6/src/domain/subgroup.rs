//! Generators of the multiplicative subgroups of a prime field.

use ark_ff::PrimeField;
use ark_std::convert::TryInto;
use num_bigint::BigUint;
use num_integer::Integer;

/// The generator of the order-`n` subgroup of `F*`, or `None` when no such subgroup exists.
///
/// `F*` is cyclic of order `q - 1`, so it has a subgroup of order `n` exactly when `n | q - 1`,
/// and `GENERATOR^((q-1)/n)` generates it. This is the general form of what arkworks' radix-2
/// domains do with `TWO_ADIC_ROOT_OF_UNITY`, and agrees with them on power-of-two `n`: that root
/// is itself `GENERATOR^((q-1)/2^TWO_ADICITY)`.
pub fn subgroup_generator<F: PrimeField>(n: usize) -> Option<F> {
    if n == 0 {
        return None;
    }
    let group_order: BigUint = Into::<BigUint>::into(F::MODULUS) - 1u8;
    let (cofactor, remainder) = group_order.div_rem(&BigUint::from(n));
    if remainder != BigUint::from(0u32) {
        return None;
    }
    let cofactor: F::BigInt = cofactor.try_into().ok()?;
    Some(F::GENERATOR.pow(cofactor))
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Field;
    use ark_std::One;

    // BW6-767's scalar field (= BLS12-381's base field): two-adicity 1.
    type Fr767 = ark_bw6_767::Fr;

    #[test]
    fn subgroup_generator_has_the_claimed_order() {
        for n in [2usize, 3, 9, 11, 23, 47, 1551, 9306] {
            let w = subgroup_generator::<Fr767>(n).unwrap();
            assert!(w.pow([n as u64]).is_one(), "w^{} != 1", n);
            // For prime n, w != 1 is enough for w to have order exactly n.
            if n > 1 && [2usize, 3, 11, 23, 47].contains(&n) {
                assert!(!w.is_one(), "generator collapsed to 1 for n = {}", n);
            }
        }
    }
}
