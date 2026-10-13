use ark_ff::{BitIteratorLE, PrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

const BITS_IN_LIMB: usize = 64;

/// A bitmask that can be encoded as a Vec of field elements.
#[derive(Clone, CanonicalSerialize, CanonicalDeserialize)]
pub struct Bitmask {
    limbs: Vec<u64>,
    padding_size: usize,
}

impl AsRef<[u64]> for Bitmask {
    fn as_ref(&self) -> &[u64] {
        &self.limbs
    }
}

// u64 limbs are used for internal representation for compatibility with ark_ff::BigInteger.
// Everything is little-endian, meaning that the least significant bit/limb is the leftmost.
// As a consequence 64-bit chunks are reversed to form equivalent u64 limbs.
// The highest limb of the bitmask is right padded with 0s.
impl Bitmask {
    pub fn from_bits(bits: &[bool]) -> Self {
        // repr = bitmask + padding
        let bitmask_size = bits.len();
        let limbs_required = bitmask_size.div_ceil(BITS_IN_LIMB);
        let repr_size = BITS_IN_LIMB * limbs_required;
        let padding_size = repr_size - bitmask_size;

        let mut repr_bits = Vec::with_capacity(repr_size);
        repr_bits.extend_from_slice(bits);
        repr_bits.extend_from_slice(&vec![false; padding_size]);

        let limbs_bits_iter = repr_bits.chunks_exact(BITS_IN_LIMB);
        assert_eq!(limbs_bits_iter.remainder().len(), 0);
        let limbs = limbs_bits_iter.map(bits_to_limb).collect();
        Self {
            limbs,
            padding_size,
        }
    }

    pub fn to_bits(&self) -> Vec<bool> {
        let mut bits: Vec<_> = BitIteratorLE::new(self).collect();
        bits.truncate(self.size()); // remove padding
        bits
    }

    pub fn to_bits_as_field_elements<F: PrimeField>(&self) -> Vec<F> {
        self.to_bits()
            .iter()
            .map(|b| if *b { F::one() } else { F::zero() })
            .collect()
    }

    pub fn size(&self) -> usize {
        let repr_size = BITS_IN_LIMB * self.limbs.len();
        // Saturating: a bitmask off the wire may claim more padding than it has bits. It then
        // reads as empty, which no keyset size matches, rather than underflowing.
        repr_size.saturating_sub(self.padding_size)
    }

    // TODO: padding
    pub fn count_ones(&self) -> usize {
        self.limbs
            .iter()
            .map(|limb| limb.count_ones() as usize)
            .sum()
    }

    /// Splits the bits into chunks of the specified size and converts the chunks into field elements,
    /// interpreting the lowest bit of the chunk as the least significant (~little-endian).
    /// Panics if the chunk doesn't have the unique representation in the field (chunk size exceeds the field capacity).
    pub fn to_chunks_as_field_elements<F: PrimeField>(&self, limbs_in_chunk: usize) -> Vec<F> {
        let bits_in_chunk = BITS_IN_LIMB * limbs_in_chunk;
        let capacity = (F::MODULUS_BIT_SIZE - 1) as usize;
        assert!(bits_in_chunk <= capacity);
        self.limbs
            .chunks(limbs_in_chunk)
            .map(limbs_to_field_elements::<F>)
            .collect()
    }
}

// Assembles a u64 limb from a little-endian bit slice of at most 64 bits.
fn bits_to_limb(bits: &[bool]) -> u64 {
    bits.iter().rev().fold(0u64, |limb_acc, next_bit| {
        (limb_acc << 1) ^ (*next_bit as u64)
    })
}

fn limbs_to_field_elements<F: PrimeField>(limbs: &[u64]) -> F {
    let mut repr = F::BigInt::default();
    let repr_limbs = repr.as_mut();
    assert!(repr_limbs.len() >= limbs.len());
    repr_limbs.iter_mut().zip(limbs).for_each(|(a, b)| *a = *b);
    #[allow(clippy::expect_used, reason = "invariant argued in the message")]
    let chunk = F::from_bigint(repr).expect(
        "invariant: callers assert the chunk has fewer bits than the modulus, so it is below it",
    );
    chunk
}

#[cfg(test)]
mod tests {
    use ark_bls12_381::Fq;
    use ark_bw6_761::Fr;
    use ark_ff::Field;
    use ark_std::convert::TryInto;
    use ark_std::test_rng;

    use crate::test_helpers::_random_bits;

    use super::*;

    pub fn _test_from_bits_to_bits(size: usize) {
        let bits = _random_bits(size, 1.0 / 2.0, &mut test_rng());
        let bitmask = Bitmask::from_bits(&bits);
        let limbs_in_bitmask = size.div_ceil(BITS_IN_LIMB);
        assert_eq!(bitmask.limbs.len(), limbs_in_bitmask);
        assert_eq!(bitmask.padding_size, BITS_IN_LIMB * limbs_in_bitmask - size);
        assert_eq!(bitmask.size(), size);
        assert_eq!(bits, bitmask.to_bits())
    }

    #[test]
    pub fn test_bitmask_from_bits_to_bits() {
        _test_from_bits_to_bits(0);
        _test_from_bits_to_bits(63);
        _test_from_bits_to_bits(64);
        _test_from_bits_to_bits(65);
        _test_from_bits_to_bits(1024);
    }

    #[test]
    pub fn test_bitmask_limbs_to_field_elements() {
        let bits = _random_bits(1024, 1.0 / 2.0, &mut test_rng());
        let bitmask = Bitmask::from_bits(&bits);
        let chunks = bitmask.to_chunks_as_field_elements::<Fr>(1);
        assert_eq!(
            chunks,
            bitmask
                .limbs
                .iter()
                .cloned()
                .map(Fr::from)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    pub fn test_bitmask_bigger_chunks_to_field_elements() {
        let bits = _random_bits(1024, 1.0 / 2.0, &mut test_rng());
        let bitmask = Bitmask::from_bits(&bits);
        let chunks = bitmask.to_chunks_as_field_elements::<Fr>(4);
        assert_eq!(chunks.len(), 4);
    }

    pub fn _test_to_field_element(set_bit_positions: Vec<usize>, expected: u128) {
        let mut bits = vec![false; 128];
        assert!(set_bit_positions.len() <= 128);
        for i in set_bit_positions {
            bits[i] = true;
        }

        let bitmask = Bitmask::from_bits(&bits);
        let chunks = bitmask.to_chunks_as_field_elements::<Fr>(2);
        assert_eq!(chunks.len(), 1);
        assert_eq!(chunks[0], Fr::from(expected));
    }

    #[test]
    pub fn test_to_field_elements() {
        _test_to_field_element(vec![0], 1);
        _test_to_field_element(vec![63], 2u128.pow(63));
        _test_to_field_element(vec![0, 64], 2u128.pow(64) + 1);
        _test_to_field_element(vec![126], 2u128.pow(126));
    }

    // ark_bls12_381::Fq::MODULUS_BIT_SIZE is 381, so the assertion in
    // `to_chunks_as_field_elements` caps chunks at 380 bits (`MODULUS_BIT_SIZE - 1`).
    // Since chunk sizes are quantized in units of `BITS_IN_LIMB` (64), the largest
    // chunk that satisfies the assertion is 5 limbs = 320 bits; 6 limbs = 384 bits
    // exceeds the capacity and must panic before it can overflow the field modulus p.

    #[test]
    pub fn test_max_capacity_chunk_unwraps_successfully() {
        let capacity: usize = (Fq::MODULUS_BIT_SIZE - 1).try_into().unwrap();
        let limbs_in_chunk = capacity / BITS_IN_LIMB;
        let bits_in_chunk = BITS_IN_LIMB * limbs_in_chunk;
        assert!(bits_in_chunk <= capacity);

        let bits = vec![true; bits_in_chunk];
        let bitmask = Bitmask::from_bits(&bits);
        let chunks = bitmask.to_chunks_as_field_elements::<Fq>(limbs_in_chunk);

        assert_eq!(chunks.len(), 1);
        let expected = Fq::from(2u8).pow([bits_in_chunk as u64]) - Fq::from(1u8);
        assert_eq!(chunks[0], expected);
    }

    #[test]
    #[should_panic]
    pub fn test_chunk_exceeding_capacity_panics() {
        let capacity: usize = (Fq::MODULUS_BIT_SIZE - 1).try_into().unwrap();
        let limbs_in_chunk = capacity / BITS_IN_LIMB + 1;
        let bits_in_chunk = BITS_IN_LIMB * limbs_in_chunk;
        assert!(bits_in_chunk > capacity);

        let bits = vec![true; bits_in_chunk];
        let bitmask = Bitmask::from_bits(&bits);
        let _ = bitmask.to_chunks_as_field_elements::<Fq>(limbs_in_chunk);
    }
}
