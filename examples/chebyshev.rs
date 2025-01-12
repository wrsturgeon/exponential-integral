//! Evaluate a Chebyshev polynomial at a point.

#![expect(
    unused_crate_dependencies,
    reason = "examples won't necessarily use each dev-dependency"
)]
#![expect(
    clippy::print_stdout,
    clippy::use_debug,
    reason = "executable, not a library"
)]

use {
    core::cmp::Ordering,
    exponential_integral::chebyshev,
    quickcheck::{Arbitrary, Gen},
    sigma_types::{Finite, Sorted, less_than::usize::LessThan},
};

fn main() {
    const N: usize = 8;

    let mut g = Gen::new(256);

    loop {
        #[expect(clippy::type_complexity, reason = "grow up")]
        let (a, b, coefficient_tuple, raw_order, x): (
            Finite<f64>,
            Finite<f64>,
            (f64, f64, f64, f64, f64, f64, f64, f64),
            usize,
            Finite<f64>,
        ) = Arbitrary::arbitrary(&mut g);

        let endpoints = Sorted::new(match a.partial_cmp(&b) {
            Some(Ordering::Less) => [a, b],
            Some(Ordering::Greater) => [b, a],
            Some(Ordering::Equal) | None => continue,
        });

        let raw_coefficients: [_; N] = coefficient_tuple.into();
        let coefficients = Finite::all(&raw_coefficients);

        #[expect(clippy::integer_division_remainder_used, reason = "not cryptographic")]
        let order = LessThan::new(raw_order % N);

        let evaluated = chebyshev::eval(endpoints, coefficients, order, x);
        println!(
            "chebyshev::eval{:#?} = {evaluated:#?}",
            (endpoints, coefficients, order, x),
        );

        return;
    }
}
