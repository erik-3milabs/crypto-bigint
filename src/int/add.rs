//! [`Int`] addition operations.

use core::ops::{Add, AddAssign};

use num_traits::WrappingAdd;
use subtle::{Choice, CtOption};

use crate::{Checked, CheckedAdd, ConstChoice, Int, Wrapping};

impl<const LIMBS: usize> Int<LIMBS> {
    /// Perform checked addition.
    pub const fn overflowing_add(&self, rhs: &Self) -> (Self, ConstChoice) {
        // Step 1. add operands
        let res = Self(self.0.wrapping_add(&rhs.0));

        // Step 2. check whether overflow happened.
        // Note:
        // - overflow can only happen when the inputs have the same sign, and then
        // - overflow occurs if and only if the result has the opposite sign of both inputs.
        //
        // We can thus express the overflow flag as: (self.msb == rhs.msb) & (self.msb != res.msb)
        let self_msb = self.is_negative();
        let overflow = self_msb
            .xor(rhs.is_negative())
            .not()
            .and(self_msb.xor(res.is_negative()));

        // Step 3. Construct result
        (res, overflow)
    }

    /// Perform wrapping addition, discarding overflow.
    pub const fn wrapping_add(&self, rhs: &Self) -> Self {
        Self(self.0.wrapping_add(&rhs.0))
    }
}

impl<const LIMBS: usize> Add for Int<LIMBS> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        self.add(&rhs)
    }
}

impl<const LIMBS: usize> Add<&Int<LIMBS>> for Int<LIMBS> {
    type Output = Self;

    fn add(self, rhs: &Self) -> Self {
        CtOption::from(self.checked_add(rhs)).expect("attempted to add with overflow")
    }
}

impl<const LIMBS: usize> AddAssign for Int<LIMBS> {
    fn add_assign(&mut self, other: Self) {
        *self += &other;
    }
}

impl<const LIMBS: usize> AddAssign<&Int<LIMBS>> for Int<LIMBS> {
    fn add_assign(&mut self, other: &Self) {
        *self = *self + other;
    }
}

impl<const LIMBS: usize> AddAssign for Wrapping<Int<LIMBS>> {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl<const LIMBS: usize> AddAssign<&Wrapping<Int<LIMBS>>> for Wrapping<Int<LIMBS>> {
    fn add_assign(&mut self, other: &Self) {
        *self = *self + other;
    }
}

impl<const LIMBS: usize> AddAssign for Checked<Int<LIMBS>> {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl<const LIMBS: usize> AddAssign<&Checked<Int<LIMBS>>> for Checked<Int<LIMBS>> {
    fn add_assign(&mut self, other: &Self) {
        *self = *self + other;
    }
}

impl<const LIMBS: usize> CheckedAdd for Int<LIMBS> {
    fn checked_add(&self, rhs: &Self) -> CtOption<Self> {
        let (val, overflow) = self.overflowing_add(rhs);
        CtOption::new(val, Choice::from(overflow.not()))
    }
}

impl<const LIMBS: usize> WrappingAdd for Int<LIMBS> {
    fn wrapping_add(&self, v: &Self) -> Self {
        self.wrapping_add(v)
    }
}

#[cfg(test)]
mod tests {

    #[cfg(test)]
    mod tests {
        use crate::int::I128;
        use crate::{ConstChoice, U128};

        #[test]
        fn overflowing_add() {
            let min_plus_one = I128 {
                0: I128::MIN.0.wrapping_add(&I128::ONE.0),
            };
            let max_minus_one = I128 {
                0: I128::MAX.0.wrapping_sub(&I128::ONE.0),
            };
            let two = I128 {
                0: U128::from(2u32),
            };

            // lhs = MIN

            let (_val, overflow) = I128::MIN.overflowing_add(&I128::MIN);
            assert_eq!(overflow, ConstChoice::TRUE);

            let (_val, overflow) = I128::MIN.overflowing_add(&I128::MINUS_ONE);
            assert_eq!(overflow, ConstChoice::TRUE);

            let (val, overflow) = I128::MIN.overflowing_add(&I128::ZERO);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::MIN);

            let (val, overflow) = I128::MIN.overflowing_add(&I128::ONE);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, min_plus_one);

            let (val, overflow) = I128::MIN.overflowing_add(&I128::MAX);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::MINUS_ONE);

            // lhs = -1

            let (_val, overflow) = I128::MINUS_ONE.overflowing_add(&I128::MIN);
            assert_eq!(overflow, ConstChoice::TRUE);

            let (val, overflow) = I128::MINUS_ONE.overflowing_add(&I128::MINUS_ONE);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, two.neg().unwrap());

            let (val, overflow) = I128::MINUS_ONE.overflowing_add(&I128::ZERO);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::MINUS_ONE);

            let (val, overflow) = I128::MINUS_ONE.overflowing_add(&I128::ONE);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::ZERO);

            let (val, overflow) = I128::MINUS_ONE.overflowing_add(&I128::MAX);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, max_minus_one);

            // lhs = 0

            let (val, overflow) = I128::ZERO.overflowing_add(&I128::MIN);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::MIN);

            let (val, overflow) = I128::ZERO.overflowing_add(&I128::MINUS_ONE);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::MINUS_ONE);

            let (val, overflow) = I128::ZERO.overflowing_add(&I128::ZERO);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::ZERO);

            let (val, overflow) = I128::ZERO.overflowing_add(&I128::ONE);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::ONE);

            let (val, overflow) = I128::ZERO.overflowing_add(&I128::MAX);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::MAX);

            // lhs = 1

            let (val, overflow) = I128::ONE.overflowing_add(&I128::MIN);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, min_plus_one);

            let (val, overflow) = I128::ONE.overflowing_add(&I128::MINUS_ONE);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::ZERO);

            let (val, overflow) = I128::ONE.overflowing_add(&I128::ZERO);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::ONE);

            let (val, overflow) = I128::ONE.overflowing_add(&I128::ONE);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, two);

            let (_val, overflow) = I128::ONE.overflowing_add(&I128::MAX);
            assert_eq!(overflow, ConstChoice::TRUE);

            // lhs = MAX

            let (val, overflow) = I128::MAX.overflowing_add(&I128::MIN);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::MINUS_ONE);

            let (val, overflow) = I128::MAX.overflowing_add(&I128::MINUS_ONE);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, max_minus_one);

            let (val, overflow) = I128::MAX.overflowing_add(&I128::ZERO);
            assert_eq!(overflow, ConstChoice::FALSE);
            assert_eq!(val, I128::MAX);

            let (_val, overflow) = I128::MAX.overflowing_add(&I128::ONE);
            assert_eq!(overflow, ConstChoice::TRUE);

            let (_val, overflow) = I128::MAX.overflowing_add(&I128::MAX);
            assert_eq!(overflow, ConstChoice::TRUE);
        }
    }
}
