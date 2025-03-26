use crate::{ConstChoice, Uint};
use num_traits::WrappingSub;

/// The matrix representation used in the partial extended gcd algorithm.
#[derive(Debug, PartialEq)]
pub struct PxgcdMatrix<const LIMBS: usize> {
    m00: Uint<LIMBS>,
    m01: Uint<LIMBS>,
    m10: Uint<LIMBS>,
    m11: Uint<LIMBS>,
    pattern: ConstChoice,
}

type Vector<const LIMBS: usize> = (Uint<LIMBS>, Uint<LIMBS>);

impl<const LIMBS: usize> PxgcdMatrix<LIMBS> {
    /// The unit / identity matrix
    const UNIT: Self = Self {
        m00: Uint::ONE,
        m01: Uint::ZERO,
        m10: Uint::ZERO,
        m11: Uint::ONE,
        pattern: ConstChoice::TRUE,
    };

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

    /// Swap the rows of this matrix if `swap` is truthy. Otherwise, do nothing.
    fn swap_rows_if(&mut self, swap: ConstChoice) {
        Uint::conditional_swap(&mut self.m00, &mut self.m10, swap);
        Uint::conditional_swap(&mut self.m01, &mut self.m11, swap);
        self.pattern = self.pattern.xor(swap);
    }

    /// Left-multiply this matrix with the upper triangular matrix that has `-2^k` as its top-right
    /// element whenever `mul` is truthy. Otherwise, do nothing. In other words, compute
    /// ```text
    /// [ 1 -2^k ] [ m00 m01 ]
    /// [ 0    1 ] [ m10 m11 ]
    /// ```
    /// if `mul.is_true()`, and nothing otherwise.
    pub(crate) fn conditional_left_mul_by_upper_triangular_matrix_with_negative_2k(
        &mut self,
        k: u32,
        mul: ConstChoice,
    ) {
        // Note: these are additions (and not subtractions like the matrix representation suggests)
        // because of how the sign-information of this matrix is stored in `pattern`.
        self.m00 = Uint::select(&self.m00, &self.m00.wrapping_add(&self.m10.shl(k)), mul);
        self.m01 = Uint::select(&self.m01, &self.m01.wrapping_add(&self.m11.shl(k)), mul);
    }
}

impl<const LIMBS: usize> Uint<LIMBS> {
    /// Extended GCD-reduce `self` and `rhs` until both can be represented using `threshold` bits.
    pub fn partial_xgcd(
        &self,
        rhs: &Uint<LIMBS>,
        threshold: u32,
    ) -> (Self, Self, PxgcdMatrix<LIMBS>) {
        let (mut a, mut b) = (*self, *rhs);
        let mut matrix = PxgcdMatrix::UNIT;

        loop {
            // Make sure a >= b
            let a_lt_b = Uint::lt(&a, &b);
            Uint::conditional_swap(&mut a, &mut b, a_lt_b);
            matrix.swap_rows_if(a_lt_b);

            // loop invariant
            debug_assert!(a >= b);

            let a_bits = a.bits();
            if a_bits <= threshold {
                break;
            }

            // Find the largest k such that a >= b*2^k
            let b_bits = b.bits();
            let mut k = a_bits - b_bits;
            let mut b_2k = b.shl(k);
            if b_2k > a {
                k -= 1;
                b_2k = b_2k.shr_vartime(1)
            }

            // Subtract b*2^k from a
            a = a.wrapping_sub(&b_2k);
            matrix.conditional_left_mul_by_upper_triangular_matrix_with_negative_2k(
                k,
                ConstChoice::TRUE,
            );
        }

        // Make sure the matrix has a positive determinant
        Uint::conditional_swap(&mut a, &mut b, matrix.pattern.not());
        matrix.swap_rows_if(matrix.pattern.not());

        (a, b, matrix)
    }
}

#[cfg(test)]
mod tests {
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
            matrix.conditional_left_mul_by_upper_triangular_matrix_with_negative_2k(
                1,
                ConstChoice::FALSE,
            );
            assert_eq!(matrix, MATRIX);

            matrix.conditional_left_mul_by_upper_triangular_matrix_with_negative_2k(
                1,
                ConstChoice::TRUE,
            );
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
    }

    mod test_pxgcd {
        use crate::modular::partial_xgcd::PxgcdMatrix;
        use crate::{ConstChoice, U64};

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
        fn text_partial_xgcd_non_unitary_elements() {
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
        fn text_partial_xgcd() {
            let threshold = 6;

            let (a, b) = (U64::from(554u64), U64::from(3321u64));
            let (partial_a, partial_b, matrix) = a.partial_xgcd(&b, threshold);

            assert_eq!(partial_a, U64::from(3u64));
            assert_eq!(partial_b, U64::from(23u64));

            assert_eq!(matrix.wrapping_apply((a, b)), (partial_a, partial_b));
        }
    }
}
