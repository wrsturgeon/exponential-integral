mod doesnt_crash {
    mod chebyshev {
        extern crate alloc;

        use {
            crate::chebyshev, alloc::format, quickcheck::TestResult, quickcheck_macros::quickcheck,
        };

        // Chebyshev approximation can balloon out of control,
        // so it doesn't need to succeed for all inputs,
        // but only on those we give it.

        // Test our `const` implementation of `min`
        // again the standard library.
        #[quickcheck]
        fn min(a: usize, b: usize) -> TestResult {
            let c = chebyshev::min(a, b);
            let m = a.min(b);
            if c == m {
                TestResult::passed()
            } else {
                TestResult::error(format!(
                    "`chebyshev::min({a}, {b}) = {c}`, but `{a}.min({b}) = {m}`"
                ))
            }
        }
    }

    mod implementation {

        mod neg {
            use {
                crate::implementation::neg::*,
                quickcheck_macros::quickcheck,
                sigma_types::{Finite, Negative},
            };

            #[quickcheck]
            fn e1(x: Negative<Finite<f64>>) {
                _ = E1(x);
            }
        }

        mod piecewise {
            use {
                crate::{constants, implementation::piecewise::*},
                quickcheck::TestResult,
                quickcheck_macros::quickcheck,
                sigma_types::{Finite, Negative, NonZero, Positive},
            };

            #[quickcheck]
            fn neg_10(x: Negative<Finite<f64>>) -> TestResult {
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
            fn neg_4(x: Negative<Finite<f64>>) -> TestResult {
                if **x <= -10_f64 {
                    return TestResult::discard();
                }
                if **x > -4_f64 {
                    return TestResult::discard();
                }
                _ = le_neg_4(x);
                TestResult::passed()
            }

            #[quickcheck]
            fn neg_1(x: Negative<Finite<f64>>) -> TestResult {
                if **x <= -4_f64 {
                    return TestResult::discard();
                }
                if **x > -1_f64 {
                    return TestResult::discard();
                }
                _ = le_neg_1(x);
                TestResult::passed()
            }

            #[quickcheck]
            fn pos_1(x: NonZero<Finite<f64>>) -> TestResult {
                if **x <= -1_f64 {
                    return TestResult::discard();
                }
                if **x > 1_f64 {
                    return TestResult::discard();
                }
                _ = le_pos_1(x);
                TestResult::passed()
            }

            #[quickcheck]
            fn pos_4(x: Positive<Finite<f64>>) -> TestResult {
                if **x <= 1_f64 {
                    return TestResult::discard();
                }
                if **x > 4_f64 {
                    return TestResult::discard();
                }
                _ = le_pos_4(x);
                TestResult::passed()
            }

            #[quickcheck]
            fn pos_max(x: Positive<Finite<f64>>) -> TestResult {
                if **x <= 4_f64 {
                    return TestResult::discard();
                }
                if **x > constants::XMAX {
                    return TestResult::discard();
                }
                _ = le_pos_max(x);
                TestResult::passed()
            }
        }

        mod pos {
            use {
                crate::implementation::pos::*,
                quickcheck_macros::quickcheck,
                sigma_types::{Finite, Positive},
            };

            #[quickcheck]
            fn e1(x: Positive<Finite<f64>>) {
                _ = E1(x);
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

    use {
        crate::{E1, Ei},
        quickcheck::TestResult,
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

    #[quickcheck]
    fn ei_near_zero(x: NonZero<Finite<f64>>) -> TestResult {
        let Some(smaller) = Finite::try_new(**x / 1_000_000_000_000_f64).and_then(NonZero::try_new)
        else {
            return TestResult::discard();
        };
        _ = Ei(smaller);
        TestResult::passed()
    }
}
