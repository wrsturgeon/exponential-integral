mod doesnt_crash {
    use {
        crate::Ei,
        quickcheck_macros::quickcheck,
        sigma_types::{Finite, NonZero},
    };

    #[quickcheck]
    fn ei(x: NonZero<Finite<f64>>) {
        _ = Ei(x);
    }

    mod chebyshev {
        use {crate::chebyshev::*, quickcheck_macros::quickcheck, sigma_types::Finite};

        #[quickcheck]
        fn series_eval(
            a: Finite<f64>,
            b: Finite<f64>,
            (c1, c2, c3, c4, c5, c6, c7, c8): (f64, f64, f64, f64, f64, f64, f64, f64),
            order: usize,
            x: Finite<f64>,
        ) {
            let series = Series {
                a,
                b,
                coefficients: &[c1, c2, c3, c4, c5, c6, c7, c8],
                order,
            };
            _ = series.eval(x);
        }
    }

    mod implementation {
        use {
            crate::implementation::*,
            quickcheck_macros::quickcheck,
            sigma_types::{Finite, NonZero},
        };

        #[quickcheck]
        fn e1(x: NonZero<Finite<f64>>) {
            _ = E1(x);
        }

        mod piecewise {
            use {
                crate::implementation::piecewise::*,
                quickcheck_macros::quickcheck,
                sigma_types::{Finite, NonZero},
            };

            #[quickcheck]
            fn lt10(x: NonZero<Finite<f64>>) {
                _ = lt_10(x);
            }
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
}
