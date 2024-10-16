//! Random number generator support

use super::Uint;
use crate::{Int, Random, RandomBits, RandomBitsError};
use rand_core::CryptoRngCore;

impl<const LIMBS: usize> Random for Int<LIMBS> {
    /// Generate a cryptographically secure random [`Int`].
    fn random(rng: &mut impl CryptoRngCore) -> Self {
        Self(Uint::random(rng))
    }
}

impl<const LIMBS: usize> RandomBits for Int<LIMBS> {
    fn try_random_bits(
        rng: &mut impl CryptoRngCore,
        bit_length: u32,
    ) -> Result<Self, RandomBitsError> {
        Self::try_random_bits_with_precision(rng, bit_length, Self::BITS)
    }

    fn try_random_bits_with_precision(
        rng: &mut impl CryptoRngCore,
        bit_length: u32,
        bits_precision: u32,
    ) -> Result<Self, RandomBitsError> {
        Uint::try_random_bits_with_precision(rng, bit_length, bits_precision)
            .map(|val| Self(val))
    }
}
