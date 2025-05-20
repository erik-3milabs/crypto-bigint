use crate::{ConstChoice, Limb, Uint};

/// The matrix representation used in the partial extended gcd algorithm.
///
/// Instead of using Ints to represent the internal state, this matrix uses [Uint]s and a
/// [ConstChoice] indicating the pattern of the signs in the matrix. This is because during the
/// execution of the partial extended gcd algorithm, the signs of the elements in the matrix
/// always adhere to one of the following two patterns:
/// ```text
///   true            false
/// [ + - ]          [ - + ]
/// [ - + ]    or    [ + - ]
/// ```
/// The `pattern` attribute indicates which of the two patterns is active.
///
/// Secondly, the determinant of this matrix is always `+/- 1`. The `pattern` attribute moreover
/// indicates whether the sign of this determinant is `1` (`true`) or `-1` (`false`).
#[derive(Debug, PartialEq)]
pub struct PxgcdMatrix<const LIMBS: usize> {
    m00: Uint<LIMBS>,
    m01: Uint<LIMBS>,
    m10: Uint<LIMBS>,
    m11: Uint<LIMBS>,
    pattern: ConstChoice,
}

impl<const LIMBS: usize> PxgcdMatrix<LIMBS> {
    /// The unit / identity matrix
    const UNIT: Self = Self {
        m00: Uint::ONE,
        m01: Uint::ZERO,
        m10: Uint::ZERO,
        m11: Uint::ONE,
        pattern: ConstChoice::TRUE,
    };

    /// Construct the adjugate of this matrix.
    ///
    /// Returns the absolute values of the adjugate matrix elements, as well as a `ConstChoice`
    /// indicating the signs of the elements: all four values are negative whenever it is `truthy`.
    ///
    /// Recall that
    /// ```text
    ///     ( [ m00 m01 ] )   [  m11 -m01 ]
    /// Adj ( [ m10 m11 ] ) = [ -m10  m00 ]
    /// ```
    pub const fn adjugate(
        &self,
    ) -> (
        (Uint<LIMBS>, Uint<LIMBS>, Uint<LIMBS>, Uint<LIMBS>),
        ConstChoice,
    ) {
        ((self.m11, self.m01, self.m10, self.m00), self.pattern.not())
    }

    /// Swap the rows of this matrix if `swap` is truthy. Otherwise, do nothing.
    ///
    /// Note: this operation preserves the fact that `det(self) = +/-1`.
    #[inline]
    fn swap_rows_if(&mut self, swap: ConstChoice) {
        Uint::conditional_swap(&mut self.m00, &mut self.m10, swap);
        Uint::conditional_swap(&mut self.m01, &mut self.m11, swap);
        self.pattern = self.pattern.xor(swap);
    }

    /// Subtract the bottom row of this matrix `2^k` times from the top row if `sub` is set.
    /// Otherwise, do nothing.
    ///
    /// In matrix notation, execute
    /// ```text
    /// [ 1 -2^k ] [ m00 m01 ]
    /// [ 0    1 ] [ m10 m11 ]
    /// ```
    /// if `mul.is_true()`, and nothing otherwise.
    ///
    /// Note: this operation preserves the fact that `det(self) = +/-1`.
    #[inline]
    fn conditional_subtract_bottom_row_2k_times_from_top_row(&mut self, k: u32, sub: ConstChoice) {
        // Note: these are additions (and not subtractions like the function name suggests) because
        // of how the sign-information of this matrix is stored in `pattern`.
        self.m00 = Uint::select(&self.m00, &self.m00.wrapping_add(&self.m10.shl(k)), sub);
        self.m01 = Uint::select(&self.m01, &self.m01.wrapping_add(&self.m11.shl(k)), sub);
    }
}

impl<const LIMBS: usize> Uint<LIMBS> {
    /// Extended GCD-reduce `self` and `rhs` until both can be represented using `threshold` bits.
    pub fn partial_xgcd(
        &self,
        rhs: &Uint<LIMBS>,
        threshold: u32,
    ) -> (Self, Self, PxgcdMatrix<LIMBS>) {
        self.bounded_partial_xgcd(rhs, threshold, Uint::<LIMBS>::BITS)
    }

    /// Extended GCD-reduce `self` and `rhs` until both can be represented using `threshold` bits.
    ///
    /// Executes in variable time with respect to `bits_upper_bound`: an upper bound on the bit size
    /// of `self` and `rhs`. Note that this algorithm becomes constant time when a static upper
    /// bound is passed.
    pub fn bounded_partial_xgcd(
        &self,
        rhs: &Uint<LIMBS>,
        threshold: u32,
        bits_upper_bound: u32,
    ) -> (Self, Self, PxgcdMatrix<LIMBS>) {
        let (mut a, mut b) = (*self, *rhs);
        let (mut a_bits, mut b_bits) = (a.bits(), b.bits());
        let mut matrix = PxgcdMatrix::UNIT;

        // Make sure a >= b
        let a_lt_b = Uint::lt(&a, &b);
        Uint::conditional_swap(&mut a, &mut b, a_lt_b);
        a_lt_b.conditional_swap_u32(&mut a_bits, &mut b_bits);
        matrix.swap_rows_if(a_lt_b);

        // TODO: we think we can prove that 7/4 * (upperbound - threshold) + C iterations will
        //  suffice (as long as the gap between upperbound and threshold is large enough?)
        let iterations = bits_upper_bound
            .saturating_sub(threshold)
            .saturating_mul(2)
            .saturating_sub(1);
        let mut i = 0;
        while i < iterations {
            i += 1;

            debug_assert!(a >= b);
            let loop_is_active = ConstChoice::from_u32_lt(threshold, a_bits);

            // Find the largest k such that a >= b*2^k
            let mut k = a_bits - b_bits;
            let mut b_2k = b.shl(k);
            let overshoot = Uint::gt(&b_2k, &a);
            k -= overshoot.select_u32(0, 1);
            b_2k = Uint::select(&b_2k, &b_2k.shr_vartime(1), overshoot);

            // Subtract b*2^k from a
            a = Uint::select(&a, &a.wrapping_sub(&b_2k), loop_is_active);
            matrix.conditional_subtract_bottom_row_2k_times_from_top_row(k, loop_is_active);
            a_bits = a.bits();

            // Make sure a >= b
            let a_lt_b = Uint::lt(&a, &b);
            Uint::conditional_swap(&mut a, &mut b, a_lt_b);
            a_lt_b.conditional_swap_u32(&mut a_bits, &mut b_bits);
            matrix.swap_rows_if(a_lt_b);
        }

        // Make sure the matrix has a positive determinant
        Uint::conditional_swap(&mut a, &mut b, matrix.pattern.not());
        matrix.swap_rows_if(matrix.pattern.not());

        (a, b, matrix)
    }

    /// XGCD reduce `self` and `rhs` until both can be represented with `threshold` bits.
    ///
    /// Requires that the most significant bit of `self` and `rhs` is NOT set.
    ///
    /// Assumes `self` and `rhs` to be random elements, i.e., at no point in the process will `a`
    /// and `b` differ more than 63 bits in size.
    pub fn partial_xgcd_randomized(
        &self,
        rhs: &Uint<LIMBS>,
        threshold: u32,
    ) -> (Self, Self, PxgcdMatrix<LIMBS>, u32) {
        // TODO: deal with situations where a and b have their top bit set
        assert!(self.as_int().is_negative().not().to_bool_vartime());
        assert!(rhs.as_int().is_negative().not().to_bool_vartime());

        let (mut a, mut b) = (*self, *rhs);
        let mut matrix = PxgcdMatrix::UNIT;

        // Make sure a >= b
        let a_lt_b = Uint::lt(&a, &b);
        Uint::conditional_swap(&mut a, &mut b, a_lt_b);
        matrix.swap_rows_if(a_lt_b);

        // c = b * 2^k
        let mut k = a.bits() - b.bits();
        let mut c = b.shl(k);

        let mut iterations = 0;
        while a.bits() > threshold && b != Uint::ZERO {
            let double_c = c.shl_vartime(1);
            let half_c = c.shr_vartime(1);
            let (a_sub_c, borrow) = a.sbb(&c, Limb::ZERO);

            if double_c <= a {
                c = double_c;
                k = k.saturating_add(1);
            } else {
                let no_underflow = borrow == Limb::ZERO;
                if no_underflow {
                    a = a_sub_c;
                    matrix.conditional_subtract_bottom_row_2k_times_from_top_row(
                        k,
                        ConstChoice::TRUE,
                    );
                }

                c = half_c;
                k = k.saturating_sub(1);
            }

            if c < b {
                assert_eq!(k, 0);
                Self::swap(&mut a, &mut b);
                matrix.swap_rows_if(ConstChoice::TRUE);
                c = b;
                k = 0;
            }

            iterations += 1;
        }

        // Make sure the matrix has a positive determinant
        Uint::conditional_swap(&mut a, &mut b, matrix.pattern.not());
        matrix.swap_rows_if(matrix.pattern.not());

        (a, b, matrix, iterations)
    }
}

#[cfg(test)]
mod tests {
    use crate::modular::partial_xgcd::PxgcdMatrix;
    use crate::{ConstChoice, Uint};

    type Vector<const LIMBS: usize> = (Uint<LIMBS>, Uint<LIMBS>);

    impl<const LIMBS: usize> PxgcdMatrix<LIMBS> {
        /// Construct a [PxgcdMatrix].
        const fn new(
            m00: Uint<LIMBS>,
            m01: Uint<LIMBS>,
            m10: Uint<LIMBS>,
            m11: Uint<LIMBS>,
            pattern: ConstChoice,
        ) -> Self {
            Self {
                m00,
                m01,
                m10,
                m11,
                pattern,
            }
        }

        /// Wrapping apply this matrix to `vector` and return the result.
        fn wrapping_apply(&self, vector: Vector<LIMBS>) -> Vector<LIMBS> {
            let (a, b) = vector;
            (
                self.m00
                    .wrapping_mul(&a)
                    .wrapping_sub(&self.m01.wrapping_mul(&b))
                    .wrapping_neg_if(self.pattern.not()),
                self.m10
                    .wrapping_mul(&a)
                    .wrapping_sub(&self.m11.wrapping_mul(&b))
                    .wrapping_neg_if(self.pattern),
            )
        }
    }

    mod test_pxgcd_matrix {
        use crate::modular::partial_xgcd::PxgcdMatrix;
        use crate::{ConstChoice, Uint, U64};

        const MATRIX: PxgcdMatrix<{ U64::LIMBS }> = PxgcdMatrix::new(
            Uint::from_u64(3u64),
            Uint::from_u64(2u64),
            Uint::ONE,
            Uint::ZERO,
            ConstChoice::TRUE,
        );

        #[test]
        fn test_swap() {
            let mut matrix = MATRIX;
            matrix.swap_rows_if(ConstChoice::FALSE);
            assert_eq!(matrix, MATRIX);

            matrix.swap_rows_if(ConstChoice::TRUE);
            assert_eq!(
                matrix,
                PxgcdMatrix::new(
                    Uint::ONE,
                    Uint::ZERO,
                    Uint::from(3u64),
                    Uint::from(2u64),
                    ConstChoice::FALSE
                )
            );
        }

        #[test]
        fn test_conditional_upper_triangular_left_mul() {
            let mut matrix = MATRIX;
            matrix.conditional_subtract_bottom_row_2k_times_from_top_row(1, ConstChoice::FALSE);
            assert_eq!(matrix, MATRIX);

            matrix.conditional_subtract_bottom_row_2k_times_from_top_row(1, ConstChoice::TRUE);
            assert_eq!(
                matrix,
                PxgcdMatrix::new(
                    Uint::from(5u64),
                    Uint::from(2u64),
                    Uint::ONE,
                    Uint::ZERO,
                    ConstChoice::TRUE
                )
            );
        }

        #[test]
        fn test_apply() {
            let matrix = PxgcdMatrix::new(
                U64::from(8u64),
                U64::from(2u64),
                U64::from(3u64),
                U64::from(7u64),
                ConstChoice::TRUE,
            );
            let (a, b) = (U64::from(55u64), U64::from(61u64));
            let (a, b) = matrix.wrapping_apply((a, b));
            assert_eq!(a, U64::from(318u64));
            assert_eq!(b, U64::from(262u64));
        }

        #[test]
        fn test_adjugate() {
            let mut x = MATRIX;
            let (matrix, sign) = x.adjugate();
            assert_eq!(
                matrix,
                (
                    Uint::ZERO,
                    Uint::from_u64(2u64),
                    Uint::ONE,
                    Uint::from_u64(3u64)
                ),
            );
            assert_eq!(sign, ConstChoice::FALSE);

            x.swap_rows_if(ConstChoice::TRUE);
            let (matrix, sign) = x.adjugate();
            assert_eq!(
                matrix,
                (
                    Uint::from_u64(2u64),
                    Uint::ZERO,
                    Uint::from_u64(3u64),
                    Uint::ONE,
                )
            );
            assert_eq!(sign, ConstChoice::TRUE);
        }
    }

    mod test_pxgcd {
        use crate::modular::partial_xgcd::PxgcdMatrix;
        use crate::{ConstChoice, U1024, U64};

        #[test]
        fn text_partial_xgcd_unit() {
            let (a, b, matrix) = U64::ONE.partial_xgcd(&U64::ZERO, 1);
            assert_eq!(a, U64::ONE);
            assert_eq!(b, U64::ZERO);
            assert_eq!(matrix, PxgcdMatrix::UNIT);
        }

        #[test]
        fn text_partial_xgcd_unit_swapped() {
            let (a, b, matrix) = U64::ZERO.partial_xgcd(&U64::ONE, 1);
            assert_eq!(a, U64::ZERO);
            assert_eq!(b, U64::ONE);
            assert_eq!(matrix, PxgcdMatrix::UNIT);
        }

        #[test]
        fn test_partial_xgcd_non_unitary_elements() {
            let (a, b, matrix) = U64::from(2u64).partial_xgcd(&U64::ONE, 1);
            assert_eq!(a, U64::ZERO);
            assert_eq!(b, U64::ONE);
            assert_eq!(
                matrix,
                PxgcdMatrix::new(
                    U64::ONE,
                    U64::from_u64(2u64),
                    U64::ZERO,
                    U64::ONE,
                    ConstChoice::TRUE
                )
            );
        }

        #[test]
        fn test_partial_xgcd_zero() {
            let threshold = 6;

            let (a, b) = (U64::from(554u64), U64::ZERO);
            let (partial_a, partial_b, matrix) = a.partial_xgcd(&b, threshold);
            assert_eq!(partial_a, U64::from(554u64));
            assert_eq!(partial_b, U64::ZERO);
            assert_eq!(matrix.wrapping_apply((a, b)), (partial_a, partial_b));

            let (partial_a, partial_b, matrix) = b.partial_xgcd(&a, threshold);

            assert_eq!(partial_a, U64::ZERO);
            assert_eq!(partial_b, U64::from(554u64));
            assert_eq!(matrix.wrapping_apply((b, a)), (partial_a, partial_b));
        }

        #[test]
        fn test_partial_xgcd() {
            let threshold = 6;

            let (a, b) = (U64::from(554u64), U64::from(3321u64));
            let (partial_a, partial_b, matrix) = a.partial_xgcd(&b, threshold);

            assert_eq!(partial_a, U64::from(3u64));
            assert_eq!(partial_b, U64::from(23u64));

            assert_eq!(matrix.wrapping_apply((a, b)), (partial_a, partial_b));
        }

        #[test]
        fn test_partial_xgcd_large() {
            let threshold = 512;

            let (a, b) = (U1024::MAX, U1024::ONE.shl(750));
            let (partial_a, partial_b, matrix) = a.partial_xgcd(&b, threshold);

            assert!(partial_a.bits() <= threshold);
            assert!(partial_b.bits() <= threshold);

            assert_eq!(matrix.wrapping_apply((a, b)), (partial_a, partial_b));
        }

        #[test]
        fn test_partial_xgcd_no_threshold_underflow() {
            let threshold = 700;
            let (a, b) = (U64::ONE, U64::ZERO);
            let (partial_a, partial_b, matrix) = a.partial_xgcd(&b, threshold);

            assert!(partial_a.bits() <= threshold);
            assert!(partial_b.bits() <= threshold);

            assert_eq!(matrix.wrapping_apply((a, b)), (partial_a, partial_b));
        }
    }

    mod test_pxgcd_randomized {
        use crate::modular::partial_xgcd::tests::Vector;
        use crate::{Concat, ConstChoice, PxgcdMatrix, Uint, U1024, U64};
        use core::ops::Sub;

        fn test_fast_partial_xgcd_output<const LIMBS: usize, const DOUBLE: usize>(
            input: Vector<LIMBS>,
            threshold: u32,
            output: Vector<LIMBS>,
            matrix: PxgcdMatrix<LIMBS>,
        ) where
            Uint<LIMBS>: Concat<Output = Uint<DOUBLE>>,
        {
            let (res_a, res_b) = output;
            assert!(res_a.bits() <= threshold);
            assert!(res_b.bits() <= threshold);

            assert_eq!(matrix.wrapping_apply(input), output)
        }

        fn pxgcd_randomized_test<const LIMBS: usize, const DOUBLE: usize>(
            input: Vector<LIMBS>,
            threshold: u32,
        ) where
            Uint<LIMBS>: Concat<Output = Uint<DOUBLE>>,
        {
            let (a, b) = input;
            let (res_a, res_b, matrix, iterations) = a.partial_xgcd_randomized(&b, threshold);
            test_fast_partial_xgcd_output((a, b), threshold, (res_a, res_b), matrix);
            assert!(iterations < (Uint::<LIMBS>::BITS - threshold) * 5 / 2);
        }

        #[test]
        fn test_pxgcd_randomized_unit() {
            let (a, b, matrix, ..) = U64::ONE.partial_xgcd_randomized(&U64::ZERO, 1);
            assert_eq!(a, U64::ONE);
            assert_eq!(b, U64::ZERO);
            assert_eq!(matrix, PxgcdMatrix::UNIT);
        }

        #[test]
        fn test_pxgcd_randomized_unit_swapped() {
            let (a, b, matrix, ..) = U64::ZERO.partial_xgcd_randomized(&U64::ONE, 1);
            assert_eq!(a, U64::ZERO);
            assert_eq!(b, U64::ONE);
            assert_eq!(matrix, PxgcdMatrix::UNIT);
        }

        #[test]
        fn test_pxgcd_randomized_non_unitary_elements() {
            let (a, b, matrix, ..) = U64::from(2u64).partial_xgcd_randomized(&U64::ONE, 1);
            assert_eq!(a, U64::ZERO);
            assert_eq!(b, U64::ONE);
            assert_eq!(
                matrix,
                PxgcdMatrix::new(
                    U64::ONE,
                    U64::from_u64(2u64),
                    U64::ZERO,
                    U64::ONE,
                    ConstChoice::TRUE
                )
            );
        }

        #[test]
        fn test_pxgcd_randomized_zero() {
            let threshold = 6;

            let (a, b) = (U64::from(554u64), U64::ZERO);
            let (partial_a, partial_b, matrix, ..) = a.partial_xgcd_randomized(&b, threshold);
            assert_eq!(partial_a, U64::from(554u64));
            assert_eq!(partial_b, U64::ZERO);
            assert_eq!(matrix.wrapping_apply((a, b)), (partial_a, partial_b));

            let (partial_a, partial_b, matrix, ..) = b.partial_xgcd_randomized(&a, threshold);

            assert_eq!(partial_a, U64::ZERO);
            assert_eq!(partial_b, U64::from(554u64));
            assert_eq!(matrix.wrapping_apply((b, a)), (partial_a, partial_b));
        }

        #[test]
        fn test_pxgcd_randomized() {
            let threshold = 6;

            let (a, b) = (U64::from(554u64), U64::from(3321u64));
            let (partial_a, partial_b, matrix, ..) = a.partial_xgcd_randomized(&b, threshold);

            assert_eq!(partial_a, U64::from(3u64));
            assert_eq!(partial_b, U64::from(23u64));

            assert_eq!(matrix.wrapping_apply((a, b)), (partial_a, partial_b));
        }

        #[test]
        fn test_pxgcd_randomized_large() {
            let threshold = 512;

            let (a, b) = (U1024::MAX.shr1(), U1024::ONE.shl(750));
            let (partial_a, partial_b, matrix, ..) = a.partial_xgcd_randomized(&b, threshold);

            assert!(partial_a.bits() <= threshold);
            assert!(partial_b.bits() <= threshold);

            assert_eq!(matrix.wrapping_apply((a, b)), (partial_a, partial_b));
        }

        #[test]
        fn test_pxgcd_randomized_no_threshold_underflow() {
            let threshold = 700;
            let (a, b) = (U64::ONE, U64::ZERO);
            let (partial_a, partial_b, matrix, ..) = a.partial_xgcd_randomized(&b, threshold);

            assert!(partial_a.bits() <= threshold);
            assert!(partial_b.bits() <= threshold);

            assert_eq!(matrix.wrapping_apply((a, b)), (partial_a, partial_b));
        }
    }
}
