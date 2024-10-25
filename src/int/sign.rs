use subtle::CtOption;
use crate::{CtChoice, Int, Uint, Word};

impl<const LIMBS: usize> Int<LIMBS> {
    /// Returns the word of most significant [`Limb`].
    /// For the generative case where the number of limbs is zero,
    /// zeroed word is returned (which is semantically correct).
    /// This method leaks the limb length of the value, which is also OK.
    const fn most_significant_word(&self) -> Word {
        if Self::LIMBS == 0 {
            0
        } else {
            self.0.to_words()[LIMBS - 1]
        }
    }

    /// Construct new [`Int`] from an absolute value and sign.
    /// Returns `None` when the absolute value does not fit in an [`Int<LIMBS>`].
    pub fn new_from_abs_sign(
        abs: Uint<LIMBS>,
        is_negative: CtChoice,
    ) -> CtOption<Self> {
        let magnitude = Self(abs).conditional_wrapping_neg(is_negative);
        let fits = Uint::ct_lte(&abs, &Int::MAX.0).or(is_negative.and(Uint::ct_eq(&abs, &Int::MIN.0)));
        CtOption::new(magnitude, fits.into())
    }

    /// Whether this [`Int`] is negative, as a `CtChoice`.
    pub const fn is_negative(&self) -> CtChoice {
        CtChoice::from_word_msb(self.most_significant_word())
    }

    /// Whether this [`Int`] is positive, as a `CtChoice`.
    pub const fn is_positive(&self) -> CtChoice {
        self.is_negative().not().and(self.ct_is_nonzero())
    }

    /// The sign and magnitude of this [`Int`].
    pub const fn abs_sign(&self) -> (Uint<LIMBS>, CtChoice) {
        let sign = self.is_negative();
        // Note: this negate_if is safe to use, since we are negating based on self.is_negative()
        let abs = self.conditional_wrapping_neg(sign);
        (abs.0, sign)
    }

    /// The magnitude of this [`Int`].
    pub const fn abs(&self) -> Uint<LIMBS> {
        self.abs_sign().0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::I128;

    #[test]
    fn is_negative() {
        assert_eq!(I128::MIN.is_negative().to_u8(), CtChoice::TRUE.to_u8());
        assert_eq!(I128::MINUS_ONE.is_negative().to_u8(), CtChoice::TRUE.to_u8());
        assert_eq!(I128::ZERO.is_negative().to_u8(), CtChoice::FALSE.to_u8());
        assert_eq!(I128::ONE.is_negative().to_u8(), CtChoice::FALSE.to_u8());
        assert_eq!(I128::MAX.is_negative().to_u8(), CtChoice::FALSE.to_u8());

        let random_negative = I128::from_be_hex("91113333555577779999BBBBDDDDFFFF");
        assert_eq!(random_negative.is_negative().to_u8(), CtChoice::TRUE.to_u8());

        let random_positive = I128::from_be_hex("71113333555577779999BBBBDDDDFFFF");
        assert_eq!(random_positive.is_negative().to_u8(), CtChoice::FALSE.to_u8());
    }

    #[test]
    fn is_positive() {
        assert_eq!(I128::MIN.is_positive().to_u8(), CtChoice::FALSE.to_u8());
        assert_eq!(I128::MINUS_ONE.is_positive().to_u8(), CtChoice::FALSE.to_u8());
        assert_eq!(I128::ZERO.is_positive().to_u8(), CtChoice::FALSE.to_u8());
        assert_eq!(I128::ONE.is_positive().to_u8(), CtChoice::TRUE.to_u8());
        assert_eq!(I128::MAX.is_positive().to_u8(), CtChoice::TRUE.to_u8());

        let random_negative = I128::from_be_hex("deadbeefcafebabedeadbeefcafebabe");
        assert_eq!(random_negative.is_positive().to_u8(), CtChoice::FALSE.to_u8());

        let random_positive = I128::from_be_hex("0badc0dedeadc0decafebabedeadcafe");
        assert_eq!(random_positive.is_positive().to_u8(), CtChoice::TRUE.to_u8());
    }
}
