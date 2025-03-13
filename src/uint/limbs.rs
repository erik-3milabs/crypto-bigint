use crate::{ConstChoice, Limb, Uint};

/// Calculate the number of leading zero limbs in the binary representation of this number
#[inline(always)]
pub(crate) const fn leading_zero_limbs(limbs: &[Limb]) -> usize {
    let mut count = 0;
    let mut i = limbs.len();
    let mut nonzero_limb_not_encountered = ConstChoice::TRUE;
    while i > 0 {
        i -= 1;
        let limb = limbs[i];
        nonzero_limb_not_encountered = nonzero_limb_not_encountered.and(limb.is_nonzero().not());
        count += nonzero_limb_not_encountered.if_true_u32(1);
    }

    count as usize
}

/// Calculate the number of leading zero limbs in the binary representation of this number in
/// variable-time with respect to `limbs`.
#[inline(always)]
pub(crate) const fn leading_zero_limbs_vartime(limbs: &[Limb]) -> usize {
    let mut count = 0;
    let mut i = limbs.len();
    while i > 0 {
        i -= 1;
        if limbs[i].0 != 0 {
            return count;
        }
        count += 1;
    }
    count
}

impl<const LIMBS: usize> Uint<LIMBS> {
    /// Calculate the minimal number of limbs needed to represent this number.
    pub const fn limbs(&self) -> usize {
        LIMBS - leading_zero_limbs(&self.limbs)
    }

    /// Calculate the minimal number of limbs needed to represent this number.
    /// Executes in variable-time with respect to `self`.
    pub const fn limbs_vartime(&self) -> usize {
        LIMBS - leading_zero_limbs_vartime(&self.limbs)
    }
}

#[cfg(test)]
mod tests {
    use crate::U256;

    #[test]
    fn test_limbs() {
        assert_eq!(U256::ZERO.limbs(), 0);
        assert_eq!(U256::ONE.limbs(), 1);
        assert_eq!(U256::ONE.shl(127).limbs(), U256::LIMBS / 2);
        assert_eq!(U256::MAX.limbs(), U256::LIMBS);
    }

    #[test]
    fn test_limbs_vartime() {
        assert_eq!(U256::ZERO.limbs_vartime(), 0);
        assert_eq!(U256::ONE.limbs_vartime(), 1);
        assert_eq!(U256::ONE.shl(127).limbs_vartime(), U256::LIMBS / 2);
        assert_eq!(U256::MAX.limbs_vartime(), U256::LIMBS);
    }
}
