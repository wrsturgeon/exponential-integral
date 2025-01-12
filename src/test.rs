mod doesnt_crash {
    // Chebyshev approximation can balloon out of control,
    // so this function doesn't need to succeed for all inputs;
    // it needs only to succeed on the inputs we give it.
    /*
    mod chebyshev {
        use {
            crate::chebyshev,
            core::cmp::Ordering,
            quickcheck::TestResult,
            quickcheck_macros::quickcheck,
            sigma_types::{Finite, Sorted, less_than::usize::LessThan},
        };

        #[quickcheck]
        fn eval(
            a: Finite<f64>,
            b: Finite<f64>,
            coefficient_tuple: (f64, f64, f64, f64, f64, f64, f64, f64),
            raw_order: usize,
            x: Finite<f64>,
        ) -> TestResult {
            const N: usize = 8;

            let endpoints = Sorted::new(match a.partial_cmp(&b) {
                Some(Ordering::Less) => [a, b],
                Some(Ordering::Greater) => [b, a],
                Some(Ordering::Equal) | None => return TestResult::discard(),
            });

            let raw_coefficients: [_; N] = coefficient_tuple.into();
            let coefficients = Finite::all(&raw_coefficients);

            let order = LessThan::new(raw_order % N);

            _ = chebyshev::eval(endpoints, coefficients, order, x);

            TestResult::passed()
        }
    }
    */

    mod implementation {

        mod piecewise {
            use {
                crate::{constants, implementation::piecewise::*},
                quickcheck::TestResult,
                quickcheck_macros::quickcheck,
                sigma_types::{Finite, NonZero},
            };

            #[quickcheck]
            fn neg_10(x: NonZero<Finite<f64>>) -> TestResult {
                if **x < constants::NXMAX {
                    return TestResult::discard();
                }
                if **x > -10_f64 {
                    return TestResult::discard();
                }
                _ = le_neg_10(x);
                TestResult::passed()
            }

            #[quickcheck]
            fn neg_4(x: NonZero<Finite<f64>>) -> TestResult {
                if **x <= -10_f64 {
                    return TestResult::discard();
                }
                if **x > -4_f64 {
                    return TestResult::discard();
                }
                _ = le_neg_4(x);
                TestResult::passed()
            }
        }

        use {
            crate::implementation::*,
            quickcheck_macros::quickcheck,
            sigma_types::{Finite, NonZero},
        };

        #[quickcheck]
        fn e1(x: NonZero<Finite<f64>>) {
            _ = E1(x);
        }
    }

    mod with_error {
        use {
            crate::with_error::*,
            quickcheck_macros::quickcheck,
            sigma_types::{Finite, NonZero},
        };

        #[quickcheck]
        fn e1(x: NonZero<Finite<f64>>) {
            _ = E1(x);
        }

        #[quickcheck]
        fn ei(x: NonZero<Finite<f64>>) {
            _ = Ei(x);
        }
    }

    use {
        crate::Ei,
        quickcheck_macros::quickcheck,
        sigma_types::{Finite, NonZero},
    };

    #[quickcheck]
    fn ei(x: NonZero<Finite<f64>>) {
        _ = Ei(x);
    }
}
