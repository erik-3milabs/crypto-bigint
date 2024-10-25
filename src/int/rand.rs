//! Random number generator support

use rand_core::CryptoRngCore;

use crate::{Int, Random};

use super::Uint;

impl<const LIMBS: usize> Random for Int<LIMBS> {
    /// Generate a cryptographically secure random [`Int`].
    fn random(rng: &mut impl CryptoRngCore) -> Self {
        Self(Uint::random(rng))
    }
}