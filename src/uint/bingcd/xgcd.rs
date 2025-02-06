use crate::uint::bingcd::matrix::IntMatrix;
use crate::{Odd, Uint};

impl<const LIMBS: usize> Odd<Uint<LIMBS>> {
    /// Constructs a matrix `M` s.t. for `(A, B) = M(a,b)` it holds that
    /// - `gcd(A, B) = gcd(a, b)`, and
    /// - `A.bits() < a.bits()` and/or `B.bits() < b.bits()`.
    ///
    /// Moreover, it returns `log_upper_bound: u32` s.t. each element in `M` lies in the interval
    /// `(-2^log_upper_bound, 2^log_upper_bound]`.
    ///
    /// Assumes `iterations < Uint::<UPDATE_LIMBS>::BITS`.
    ///
    /// The function executes in time variable in `iterations`.
    #[inline]
    pub(super) const fn restricted_extended_gcd_vartime<const UPDATE_LIMBS: usize>(
        &self,
        rhs: &Uint<LIMBS>,
        iterations: u32,
    ) -> (IntMatrix<UPDATE_LIMBS>, u32) {
        debug_assert!(iterations < Uint::<UPDATE_LIMBS>::BITS);
        // (self, rhs) corresponds to (b_, a_) in the Algorithm 1 notation.
        let (mut a, mut b) = (*rhs, *self.as_ref());

        // Compute the update matrix.
        let mut matrix = IntMatrix::UNIT;
        let mut log_upper_bound = 0;
        let mut j = 0;
        while j < iterations {
            j += 1;

            let a_odd = a.is_odd();
            let a_lt_b = Uint::lt(&a, &b);

            // swap if a odd and a < b
            let do_swap = a_odd.and(a_lt_b);
            Uint::conditional_swap(&mut a, &mut b, do_swap);
            matrix.conditional_swap_rows(do_swap);

            // subtract b from a when a is odd
            a = Uint::select(&a, &a.wrapping_sub(&b), a_odd);
            matrix.conditional_subtract_bottom_row_from_top(a_odd);

            // Div a by 2 and double the bottom row of the matrix when a, b ≠ 0.
            let do_apply = a.is_nonzero().and(b.is_nonzero());
            // safe to vartime; shr_vartime is variable in the value of shift only. Since this shift
            // is a public constant, the constant time property of this algorithm is not impacted.
            a = Uint::select(&a, &a.shr_vartime(1), do_apply);
            matrix.conditional_double_bottom_row(do_apply);
            log_upper_bound = do_apply.select_u32(log_upper_bound, log_upper_bound + 1);
        }

        (matrix, log_upper_bound)
    }
}

#[cfg(test)]
mod tests {
    use crate::uint::bingcd::matrix::IntMatrix;
    use crate::{I64, U64};

    #[test]
    fn test_restricted_extended_gcd() {
        let a = U64::from_be_hex("CA048AFA63CD6A1F").to_odd().unwrap();
        let b = U64::from_be_hex("AE693BF7BE8E5566");
        let (matrix, iters) = a.restricted_extended_gcd_vartime(&b, 5);
        assert_eq!(iters, 5);
        assert_eq!(
            matrix,
            IntMatrix::new(I64::from(5), I64::from(-2), I64::from(-4), I64::from(8))
        );
    }

    #[test]
    fn test_restricted_extended_gcd_stops_early() {
        // Stop before max_iters
        let a = U64::from_be_hex("0000000003CD6A1F").to_odd().unwrap();
        let b = U64::from_be_hex("000000000E8E5566");
        let (.., iters) = a.restricted_extended_gcd_vartime::<{ U64::LIMBS }>(&b, 60);
        assert_eq!(iters, 35);
    }
}
