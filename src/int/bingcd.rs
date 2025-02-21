//! This module implements (a constant variant of) the Optimized Extended Binary GCD algorithm,
//! which is described by Pornin in "Optimized Binary GCD for Modular Inversion".
//! Ref: <https://eprint.iacr.org/2020/972.pdf>

use crate::modular::bingcd::tools::const_min;
use crate::{ConstChoice, Int, NonZero, Odd, Uint};

pub struct BinXgcdOutput<const LIMBS: usize> {
    gcd: Odd<Uint<LIMBS>>,
    x: Int<LIMBS>,
    y: Int<LIMBS>,
    lhs_on_gcd: Int<LIMBS>,
    rhs_on_gcd: Int<LIMBS>,
}

impl<const LIMBS: usize> BinXgcdOutput<LIMBS> {
    /// Return the quotients `lhs.gcd` and `rhs/gcd`.
    pub const fn quotients(&self) -> (Int<LIMBS>, Int<LIMBS>) {
        (self.lhs_on_gcd, self.rhs_on_gcd)
    }

    /// Return the Bézout coefficients `x` and `y` s.t. `lhs * x + rhs * y = gcd`.
    pub const fn bezout_coefficients(&self) -> (Int<LIMBS>, Int<LIMBS>) {
        (self.x, self.y)
    }
}

impl<const LIMBS: usize> Int<LIMBS> {
    /// Compute the gcd of `self` and `rhs` leveraging the Binary GCD algorithm.
    pub const fn bingcd(&self, rhs: &Self) -> Uint<LIMBS> {
        self.abs().bingcd(&rhs.abs())
    }

    /// Executes the Binary Extended GCD algorithm.
    ///
    /// Given `(self, rhs)`, computes `(g, x, y)`, s.t. `self * x + rhs * y = g = gcd(self, rhs)`.
    pub const fn binxgcd(&self, rhs: &Self) -> (Uint<LIMBS>, Self, Self) {
        // Make sure `self` and `rhs` are nonzero.
        let self_is_nz = self.is_nonzero();
        let self_nz = Int::select(&Int::ONE, self, self_is_nz)
            .to_nz()
            .expect("self is non zero by construction");
        let rhs_is_nz = rhs.is_nonzero();
        let rhs_nz = Int::select(&Int::ONE, rhs, rhs_is_nz)
            .to_nz()
            .expect("rhs is non zero by construction");

        let (gcd, mut x, mut y) = self_nz.binxgcd(&rhs_nz);

        // Correct the gcd in case self and/or rhs was zero
        let gcd = Uint::select(&rhs.abs(), gcd.as_ref(), self_is_nz);
        let gcd = Uint::select(&self.abs(), &gcd, rhs_is_nz);

        // Correct the Bézout coefficients in case self and/or rhs was zero.
        let signum_self = Int::new_from_abs_sign(Uint::ONE, self.is_negative()).expect("+/- 1");
        let signum_rhs = Int::new_from_abs_sign(Uint::ONE, rhs.is_negative()).expect("+/- 1");
        x = Int::select(&Int::ZERO, &x, self_is_nz);
        y = Int::select(&signum_rhs, &y, self_is_nz);
        x = Int::select(&signum_self, &x, rhs_is_nz);
        y = Int::select(&Int::ZERO, &y, rhs_is_nz);

        (gcd, x, y)
    }
}

impl<const LIMBS: usize> NonZero<Int<LIMBS>> {
    /// Compute the gcd of `self` and `rhs` leveraging the Binary GCD algorithm.
    pub const fn bingcd(&self, rhs: &Self) -> NonZero<Uint<LIMBS>> {
        self.abs().bingcd(&rhs.as_ref().abs())
    }

    /// Execute the Binary Extended GCD algorithm.
    ///
    /// Given `(self, rhs)`, computes `(g, x, y)` s.t. `self * x + rhs * y = g = gcd(self, rhs)`.
    pub const fn binxgcd(&self, rhs: &Self) -> (NonZero<Uint<LIMBS>>, Int<LIMBS>, Int<LIMBS>) {
        let (mut lhs, mut rhs) = (*self.as_ref(), *rhs.as_ref());

        // Leverage the property that gcd(2^k * a, 2^k *b) = 2^k * gcd(a, b)
        let i = lhs.0.trailing_zeros();
        let j = rhs.0.trailing_zeros();
        let k = const_min(i, j);
        lhs = lhs.shr(k);
        rhs = rhs.shr(k);

        // Note: at this point, either lhs or rhs is odd (or both).
        // Swap to make sure lhs is odd.
        let swap = ConstChoice::from_u32_lt(j, i);
        Int::conditional_swap(&mut lhs, &mut rhs, swap);
        let lhs = lhs.to_odd().expect("odd by construction");

        let rhs = rhs.to_nz().expect("non-zero by construction");
        let output = lhs.binxgcd(&rhs);
        let gcd = output.gcd;
        let (mut x, mut y) = output.bezout_coefficients();

        Int::conditional_swap(&mut x, &mut y, swap);

        // Add the factor 2^k to the gcd.
        let gcd = gcd
            .as_ref()
            .shl(k)
            .to_nz()
            .expect("gcd of non-zero element with another element is non-zero");

        (gcd, x, y)
    }
}

impl<const LIMBS: usize> Odd<Int<LIMBS>> {
    /// Compute the gcd of `self` and `rhs` leveraging the Binary GCD algorithm.
    pub const fn bingcd(&self, rhs: &Self) -> Odd<Uint<LIMBS>> {
        self.abs().bingcd(&rhs.as_ref().abs())
    }

    /// Execute the Binary Extended GCD algorithm.
    ///
    /// Given `(self, rhs)`, computes `(g, x, y)` s.t. `self * x + rhs * y = g = gcd(self, rhs)`.
    pub const fn binxgcd(&self, rhs: &NonZero<Int<LIMBS>>) -> BinXgcdOutput<LIMBS> {
        let (abs_lhs, sgn_lhs) = self.abs_sign();
        let (abs_rhs, sgn_rhs) = rhs.abs_sign();

        let output = abs_lhs.binxgcd_nz(&abs_rhs);

        let (mut x, mut y) = output.bezout_coefficients();
        x = x.wrapping_neg_if(sgn_lhs);
        y = y.wrapping_neg_if(sgn_rhs);

        let (abs_lhs_on_gcd, abs_rhs_on_gcd) = output.quotients();
        let lhs_on_gcd = Int::new_from_abs_sign(abs_lhs_on_gcd, sgn_lhs).expect("no overflow");
        let rhs_on_gcd = Int::new_from_abs_sign(abs_rhs_on_gcd, sgn_rhs).expect("no overflow");

        BinXgcdOutput {
            gcd: output.gcd,
            x,
            y,
            lhs_on_gcd,
            rhs_on_gcd,
        }
    }
}

#[cfg(test)]
mod test {

    mod test_int_binxgcd {
        use crate::{
            ConcatMixed, Gcd, Int, Random, Uint, U1024, U128, U192, U2048, U256, U384, U4096, U512,
            U64, U768, U8192,
        };
        use rand_core::OsRng;

        fn int_binxgcd_test<const LIMBS: usize, const DOUBLE: usize>(
            lhs: Int<LIMBS>,
            rhs: Int<LIMBS>,
        ) where
            Uint<LIMBS>: ConcatMixed<Uint<LIMBS>, MixedOutput = Uint<DOUBLE>>,
            Int<LIMBS>: Gcd<Output = Uint<LIMBS>>,
        {
            let gcd = lhs.gcd(&rhs);
            let (xgcd, x, y) = lhs.binxgcd(&rhs);
            assert_eq!(gcd, xgcd);
            let x_lhs = x.widening_mul(&lhs);
            let y_rhs = y.widening_mul(&rhs);
            let prod = x_lhs.wrapping_add(&y_rhs);
            assert_eq!(prod, xgcd.resize().as_int());
        }

        fn int_binxgcd_tests<const LIMBS: usize, const DOUBLE: usize>()
        where
            Uint<LIMBS>: ConcatMixed<Uint<LIMBS>, MixedOutput = Uint<DOUBLE>>,
            Int<LIMBS>: Gcd<Output = Uint<LIMBS>>,
        {
            int_binxgcd_test(Int::MIN, Int::MIN);
            int_binxgcd_test(Int::MIN, Int::MINUS_ONE);
            int_binxgcd_test(Int::MIN, Int::ZERO);
            int_binxgcd_test(Int::MIN, Int::ONE);
            int_binxgcd_test(Int::MIN, Int::MAX);
            int_binxgcd_test(Int::ONE, Int::MIN);
            int_binxgcd_test(Int::ONE, Int::MINUS_ONE);
            int_binxgcd_test(Int::ONE, Int::ZERO);
            int_binxgcd_test(Int::ONE, Int::ONE);
            int_binxgcd_test(Int::ONE, Int::MAX);
            int_binxgcd_test(Int::ZERO, Int::MIN);
            int_binxgcd_test(Int::ZERO, Int::MINUS_ONE);
            int_binxgcd_test(Int::ZERO, Int::ZERO);
            int_binxgcd_test(Int::ZERO, Int::ONE);
            int_binxgcd_test(Int::ZERO, Int::MAX);
            int_binxgcd_test(Int::ONE, Int::MIN);
            int_binxgcd_test(Int::ONE, Int::MINUS_ONE);
            int_binxgcd_test(Int::ONE, Int::ZERO);
            int_binxgcd_test(Int::ONE, Int::ONE);
            int_binxgcd_test(Int::ONE, Int::MAX);
            int_binxgcd_test(Int::MAX, Int::MIN);
            int_binxgcd_test(Int::MAX, Int::MINUS_ONE);
            int_binxgcd_test(Int::MAX, Int::ZERO);
            int_binxgcd_test(Int::MAX, Int::ONE);
            int_binxgcd_test(Int::MAX, Int::MAX);

            for _ in 0..100 {
                let x = Int::random(&mut OsRng);
                let y = Int::random(&mut OsRng);
                int_binxgcd_test(x, y);
            }
        }

        #[test]
        fn test_int_binxgcd() {
            int_binxgcd_tests::<{ U64::LIMBS }, { U128::LIMBS }>();
            int_binxgcd_tests::<{ U128::LIMBS }, { U256::LIMBS }>();
            int_binxgcd_tests::<{ U192::LIMBS }, { U384::LIMBS }>();
            int_binxgcd_tests::<{ U256::LIMBS }, { U512::LIMBS }>();
            int_binxgcd_tests::<{ U384::LIMBS }, { U768::LIMBS }>();
            int_binxgcd_tests::<{ U512::LIMBS }, { U1024::LIMBS }>();
            int_binxgcd_tests::<{ U1024::LIMBS }, { U2048::LIMBS }>();
            int_binxgcd_tests::<{ U2048::LIMBS }, { U4096::LIMBS }>();
            int_binxgcd_tests::<{ U4096::LIMBS }, { U8192::LIMBS }>();
        }
    }

    mod test_nonzero_int_binxgcd {
        use crate::{
            ConcatMixed, Gcd, Int, NonZero, RandomMod, Uint, U1024, U128, U192, U2048, U256, U384,
            U4096, U512, U64, U768, U8192,
        };
        use rand_core::OsRng;

        fn nz_int_binxgcd_test<const LIMBS: usize, const DOUBLE: usize>(
            lhs: NonZero<Int<LIMBS>>,
            rhs: NonZero<Int<LIMBS>>,
        ) where
            Uint<LIMBS>: ConcatMixed<Uint<LIMBS>, MixedOutput = Uint<DOUBLE>>,
            Int<LIMBS>: Gcd<Output = Uint<LIMBS>>,
        {
            let gcd = lhs.gcd(&rhs);
            let (xgcd, x, y) = lhs.binxgcd(&rhs);
            assert_eq!(gcd.to_nz().unwrap(), xgcd);
            assert_eq!(
                x.widening_mul(&lhs).wrapping_add(&y.widening_mul(&rhs)),
                xgcd.resize().as_int()
            );
        }

        fn nz_int_binxgcd_tests<const LIMBS: usize, const DOUBLE: usize>()
        where
            Uint<LIMBS>: ConcatMixed<Uint<LIMBS>, MixedOutput = Uint<DOUBLE>>,
            Int<LIMBS>: Gcd<Output = Uint<LIMBS>>,
        {
            let nz_min = Int::MIN.to_nz().expect("is nz");
            let nz_minus_one = Int::MINUS_ONE.to_nz().expect("is nz");
            let nz_one = Int::ONE.to_nz().expect("is nz");
            let nz_max = Int::MAX.to_nz().expect("is nz");

            nz_int_binxgcd_test(nz_min, nz_min);
            nz_int_binxgcd_test(nz_min, nz_minus_one);
            nz_int_binxgcd_test(nz_min, nz_one);
            nz_int_binxgcd_test(nz_min, nz_max);
            nz_int_binxgcd_test(nz_one, nz_min);
            nz_int_binxgcd_test(nz_one, nz_minus_one);
            nz_int_binxgcd_test(nz_one, nz_one);
            nz_int_binxgcd_test(nz_one, nz_max);
            nz_int_binxgcd_test(nz_max, nz_min);
            nz_int_binxgcd_test(nz_max, nz_minus_one);
            nz_int_binxgcd_test(nz_max, nz_one);
            nz_int_binxgcd_test(nz_max, nz_max);

            let bound = Int::MIN.abs().to_nz().unwrap();
            for _ in 0..100 {
                let x = Uint::random_mod(&mut OsRng, &bound)
                    .as_int()
                    .to_nz()
                    .unwrap();
                let y = Uint::random_mod(&mut OsRng, &bound)
                    .as_int()
                    .to_nz()
                    .unwrap();
                nz_int_binxgcd_test(x, y);
            }
        }

        #[test]
        fn test_nz_int_binxgcd() {
            nz_int_binxgcd_tests::<{ U64::LIMBS }, { U128::LIMBS }>();
            nz_int_binxgcd_tests::<{ U128::LIMBS }, { U256::LIMBS }>();
            nz_int_binxgcd_tests::<{ U192::LIMBS }, { U384::LIMBS }>();
            nz_int_binxgcd_tests::<{ U256::LIMBS }, { U512::LIMBS }>();
            nz_int_binxgcd_tests::<{ U384::LIMBS }, { U768::LIMBS }>();
            nz_int_binxgcd_tests::<{ U512::LIMBS }, { U1024::LIMBS }>();
            nz_int_binxgcd_tests::<{ U1024::LIMBS }, { U2048::LIMBS }>();
            nz_int_binxgcd_tests::<{ U2048::LIMBS }, { U4096::LIMBS }>();
            nz_int_binxgcd_tests::<{ U4096::LIMBS }, { U8192::LIMBS }>();
        }
    }

    mod test_odd_int_binxgcd {
        use crate::{
            ConcatMixed, Int, Random, Uint, U1024, U128, U192, U2048, U256, U384, U4096, U512, U64,
            U768, U8192,
        };
        use rand_core::OsRng;

        fn odd_int_binxgcd_test<const LIMBS: usize, const DOUBLE: usize>(
            lhs: Int<LIMBS>,
            rhs: Int<LIMBS>,
        ) where
            Uint<LIMBS>: ConcatMixed<Uint<LIMBS>, MixedOutput = Uint<DOUBLE>>,
        {
            let gcd = lhs.bingcd(&rhs);
            let output = lhs.to_odd().unwrap().binxgcd(&rhs.to_nz().unwrap());
            assert_eq!(gcd.to_odd().unwrap(), output.gcd);

            // Test quotients
            let (lhs_on_gcd, rhs_on_gcd) = output.quotients();
            assert_eq!(lhs_on_gcd, lhs.div_uint(&gcd.to_nz().unwrap()));
            assert_eq!(rhs_on_gcd, rhs.div_uint(&gcd.to_nz().unwrap()));

            // Test the Bezout coefficients
            let (x, y) = output.bezout_coefficients();
            assert_eq!(
                x.widening_mul(&lhs).wrapping_add(&y.widening_mul(&rhs)),
                gcd.resize().as_int()
            );
        }

        fn odd_int_binxgcd_tests<const LIMBS: usize, const DOUBLE: usize>()
        where
            Uint<LIMBS>: ConcatMixed<Uint<LIMBS>, MixedOutput = Uint<DOUBLE>>,
        {
            let neg_max = Int::MAX.wrapping_neg();
            odd_int_binxgcd_test(neg_max, neg_max);
            odd_int_binxgcd_test(neg_max, Int::MINUS_ONE);
            odd_int_binxgcd_test(neg_max, Int::ONE);
            odd_int_binxgcd_test(neg_max, Int::MAX);
            odd_int_binxgcd_test(Int::ONE, neg_max);
            odd_int_binxgcd_test(Int::ONE, Int::MINUS_ONE);
            odd_int_binxgcd_test(Int::ONE, Int::ONE);
            odd_int_binxgcd_test(Int::ONE, Int::MAX);
            odd_int_binxgcd_test(Int::MAX, neg_max);
            odd_int_binxgcd_test(Int::MAX, Int::MINUS_ONE);
            odd_int_binxgcd_test(Int::MAX, Int::ONE);
            odd_int_binxgcd_test(Int::MAX, Int::MAX);

            for _ in 0..100 {
                let x = Int::<LIMBS>::random(&mut OsRng).bitor(&Int::ONE);
                let y = Int::<LIMBS>::random(&mut OsRng);
                odd_int_binxgcd_test(x, y);
            }
        }

        #[test]
        fn test_odd_int_binxgcd() {
            odd_int_binxgcd_tests::<{ U64::LIMBS }, { U128::LIMBS }>();
            odd_int_binxgcd_tests::<{ U128::LIMBS }, { U256::LIMBS }>();
            odd_int_binxgcd_tests::<{ U192::LIMBS }, { U384::LIMBS }>();
            odd_int_binxgcd_tests::<{ U256::LIMBS }, { U512::LIMBS }>();
            odd_int_binxgcd_tests::<{ U384::LIMBS }, { U768::LIMBS }>();
            odd_int_binxgcd_tests::<{ U512::LIMBS }, { U1024::LIMBS }>();
            odd_int_binxgcd_tests::<{ U1024::LIMBS }, { U2048::LIMBS }>();
            odd_int_binxgcd_tests::<{ U2048::LIMBS }, { U4096::LIMBS }>();
            odd_int_binxgcd_tests::<{ U4096::LIMBS }, { U8192::LIMBS }>();
        }
    }
}
