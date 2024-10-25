//! [`Int`] bitwise right shift operations.

use core::ops::{Shr, ShrAssign};
use subtle::{CtOption};
use crate::{CtChoice, Int, Uint};

impl<const LIMBS: usize> Int<LIMBS> {
    /// Computes `self >> shift`.
    ///
    /// Note, this is signed shift right, i.e., the value shifted in on the left is equal to
    /// the most significant bit.
    ///
    /// Panics if `shift >= Self::BITS`.
    pub fn shr(&self, shift: usize) -> Self {
        self.overflowing_shr(shift)
            .expect("`shift` within the bit size of the integer")
    }

    /// Computes `self >> shift` in variable time.
    ///
    /// Note, this is signed shift right, i.e., the value shifted in on the left is equal to
    /// the most significant bit.
    ///
    /// Panics if `shift >= Self::BITS`.
    pub fn shr_vartime(&self, shift: usize) -> Self {
        self.overflowing_shr_vartime(shift)
            .expect("`shift` within the bit size of the integer")
    }

    /// Computes `self >> shift`.
    ///
    /// Note, this is signed shift right, i.e., the value shifted in on the left is equal to
    /// the most significant bit.
    ///
    /// Returns `None` if `shift >= Self::BITS`.
    pub fn overflowing_shr(&self, shift: usize) -> CtOption<Self> {
        // `floor(log2(BITS - 1))` is the number of bits in the representation of `shift`
        // (which lies in range `0 <= shift < BITS`).
        let shift_bits = usize::BITS - (Self::BITS - 1).leading_zeros();
        let overflow = CtChoice::from_usize_lt(shift, Self::BITS).not();
        let shift = shift % Self::BITS;
        let mut result = *self;
        let mut i = 0;
        while i < shift_bits {
            let bit = CtChoice::from_usize_lsb((shift >> i) & 1);
            result = Int::select(
                &result,
                &result
                    .overflowing_shr_vartime(1 << i)
                    .expect("shift within range"),
                bit,
            );
            i += 1;
        }

        CtOption::new(result, overflow.not().into())
    }

    /// Computes `self >> shift`.
    ///
    /// Returns `None` if `shift >= Self::BITS`.
    ///
    /// NOTE: this operation is variable time with respect to `shift` *ONLY*.
    ///
    /// When used with a fixed `shift`, this function is constant-time with respect
    /// to `self`.
    #[inline(always)]
    pub fn overflowing_shr_vartime(&self, shift: usize) -> CtOption<Self> {
        if shift >= Self::BITS {
            return CtOption::new(Self::ZERO, CtChoice::FALSE.into());
        }

        // safe to unwrap, due to above check
        let logical_shr = Self(Uint::shr_vartime(&self.0, shift));

        // To turn a logical shr into an arithmetical, the shifted in bits need to match the
        // msb of self.

        let masked_msb = self.bitand(&Self::SIGN_MASK);
        let inverse_shift = shift.saturating_sub(1);
        let shifted_masked_msb = Uint::shr_vartime(&masked_msb.0, inverse_shift);
        let msbs = Self(shifted_masked_msb).negc().0;

        CtOption::new(logical_shr.bitxor(&msbs), CtChoice::TRUE.into())
    }

    /// Computes `self >> shift` in a panic-free manner, returning zero if the shift exceeds the
    /// precision.
    pub fn wrapping_shr(&self, shift: usize) -> Self {
        self.overflowing_shr(shift).unwrap_or(Self::ZERO)
    }

    /// Computes `self >> shift` in variable-time in a panic-free manner, returning zero if the
    /// shift exceeds the precision.
    pub fn wrapping_shr_vartime(&self, shift: usize) -> Self {
        self.overflowing_shr_vartime(shift).unwrap_or(Self::ZERO)
    }
}

impl<const LIMBS: usize> Shr<usize> for Int<LIMBS> {
    type Output = Int<LIMBS>;

    /// NOTE: this operation is variable time with respect to `rhs` *ONLY*.
    ///
    /// When used with a fixed `rhs`, this function is constant-time with respect
    /// to `self`.
    fn shr(self, rhs: usize) -> Int<LIMBS> {
        Int::<LIMBS>::shr(&self, rhs)
    }
}

impl<const LIMBS: usize> Shr<usize> for &Int<LIMBS> {
    type Output = Int<LIMBS>;

    /// NOTE: this operation is variable time with respect to `rhs` *ONLY*.
    ///
    /// When used with a fixed `rhs`, this function is constant-time with respect
    /// to `self`.
    fn shr(self, rhs: usize) -> Int<LIMBS> {
        self.shr(rhs)
    }
}

impl<const LIMBS: usize> ShrAssign<usize> for Int<LIMBS> {
    fn shr_assign(&mut self, rhs: usize) {
        *self = self.shr(rhs);
    }
}


#[cfg(test)]
mod tests {
    use core::ops::Div;

    use crate::I256;

    const N: I256 =
        I256::from_be_hex("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141");

    const N_2: I256 =
        I256::from_be_hex("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF5D576E7357A4501DDFE92F46681B20A0");

    #[test]
    fn shr0() {
        assert_eq!(I256::MAX >> 0, I256::MAX);
        assert_eq!(I256::MIN >> 0, I256::MIN);
    }

    #[test]
    fn shr1() {
        assert_eq!(N >> 1, N_2);
    }

    #[test]
    fn shr5() {
        assert_eq!(
            I256::MAX >> 5,
            I256::MAX.div(I256::from(32).to_nz().unwrap()).unwrap()
        );
        assert_eq!(
            I256::MIN >> 5,
            I256::MIN.div(I256::from(32).to_nz().unwrap()).unwrap()
        );
    }

    #[test]
    fn shr7_vartime() {
        assert_eq!(
            I256::MAX.shr_vartime(7),
            I256::MAX.div(I256::from(128).to_nz().unwrap()).unwrap()
        );
        assert_eq!(
            I256::MIN.shr_vartime(7),
            I256::MIN.div(I256::from(128).to_nz().unwrap()).unwrap()
        );
    }

    #[test]
    fn shr256_const() {
        assert_eq!(N.overflowing_shr(256).is_none().unwrap_u8(), 1);
        assert_eq!(N.overflowing_shr_vartime(256).is_none().unwrap_u8(), 1);
    }

    #[test]
    #[should_panic(expected = "`shift` within the bit size of the integer")]
    fn shr256() {
        let _ = N >> 256;
    }
}
