//! [`Int`] comparisons.
//!
//! By default, these are all constant-time and use the `subtle` crate.
#![allow(dead_code)]

use core::cmp::Ordering;

use subtle::{Choice, ConstantTimeEq, ConstantTimeGreater, ConstantTimeLess};

use crate::{ConstChoice, Int, Uint};

impl<const LIMBS: usize> Int<LIMBS> {
    /// Return `b` if `c` is truthy, otherwise return `a`.
    #[inline]
    pub(crate) const fn select(a: &Self, b: &Self, c: ConstChoice) -> Self {
        Self(Uint::select(&a.0, &b.0, c))
    }

    /// Returns the truthy value if `self`!=0 or the falsy value otherwise.
    #[inline]
    pub(crate) const fn is_nonzero(&self) -> ConstChoice {
        Uint::is_nonzero(&self.0)
    }

    /// Returns the truthy value if `self` is odd or the falsy value otherwise.
    pub(crate) const fn is_odd(&self) -> ConstChoice {
        Uint::is_odd(&self.0)
    }

    /// Returns the truthy value if `self == rhs` or the falsy value otherwise.
    #[inline]
    pub(crate) const fn eq(lhs: &Self, rhs: &Self) -> ConstChoice {
        Uint::eq(&lhs.0, &rhs.0)
    }

    /// Returns the truthy value if `self < rhs` and the falsy value otherwise.
    #[inline]
    pub(crate) const fn lt(lhs: &Self, rhs: &Self) -> ConstChoice {
        Uint::lt(&lhs.invert_msb().0, &rhs.invert_msb().0)
    }

    /// Returns the truthy value if `self > rhs` and the falsy value otherwise.
    #[inline]
    pub(crate) const fn gt(lhs: &Self, rhs: &Self) -> ConstChoice {
        Uint::gt(&lhs.invert_msb().0, &rhs.invert_msb().0)
    }

    /// Returns the ordering between `self` and `rhs` as an i8.
    /// Values correspond to the Ordering enum:
    ///   -1 is Less
    ///   0 is Equal
    ///   1 is Greater
    #[inline]
    pub(crate) const fn cmp(lhs: &Self, rhs: &Self) -> i8 {
        Uint::cmp(&lhs.invert_msb().0, &rhs.invert_msb().0)
    }

    /// Returns the Ordering between `self` and `rhs` in variable time.
    pub const fn cmp_vartime(&self, rhs: &Self) -> Ordering {
        self.invert_msb().0.cmp_vartime(&rhs.invert_msb().0)
    }
}

impl<const LIMBS: usize> ConstantTimeEq for Int<LIMBS> {
    #[inline]
    fn ct_eq(&self, other: &Self) -> Choice {
        Int::eq(self, other).into()
    }
}

impl<const LIMBS: usize> ConstantTimeGreater for Int<LIMBS> {
    #[inline]
    fn ct_gt(&self, other: &Self) -> Choice {
        Int::gt(self, other).into()
    }
}

impl<const LIMBS: usize> ConstantTimeLess for Int<LIMBS> {
    #[inline]
    fn ct_lt(&self, other: &Self) -> Choice {
        Int::lt(self, other).into()
    }
}

impl<const LIMBS: usize> Eq for Int<LIMBS> {}

impl<const LIMBS: usize> Ord for Int<LIMBS> {
    fn cmp(&self, other: &Self) -> Ordering {
        let c = Self::cmp(self, other);
        match c {
            -1 => Ordering::Less,
            0 => Ordering::Equal,
            _ => Ordering::Greater,
        }
    }
}

impl<const LIMBS: usize> PartialOrd for Int<LIMBS> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<const LIMBS: usize> PartialEq for Int<LIMBS> {
    fn eq(&self, other: &Self) -> bool {
        self.ct_eq(other).into()
    }
}
