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
    exponential_integral::chebyshev,
    quickcheck::{Arbitrary, Gen},
    sigma_types::Finite,
};

#[cfg(feature = "precision")]
use sigma_types::less_than::usize::LessThan;

fn main() {
    const N: usize = 8;

    let mut g = Gen::new(256);

    #[expect(clippy::type_complexity, reason = "grow up")]
    let (coefficient_tuple, x): ((f64, f64, f64, f64, f64, f64, f64, f64), Finite<f64>) =
        Arbitrary::arbitrary(&mut g);

    #[cfg(feature = "precision")]
    let raw_order = usize::arbitrary(&mut g);

    let raw_coefficients: [_; N] = coefficient_tuple.into();
    let coefficients = Finite::all(&raw_coefficients);

    #[cfg(feature = "precision")]
    #[expect(clippy::integer_division_remainder_used, reason = "not cryptographic")]
    let order = LessThan::new(raw_order % N);

    let evaluated = chebyshev::eval(
        coefficients,
        x,
        #[cfg(feature = "precision")]
        order,
    );
    println!(
        "chebyshev::eval{:#?} = {evaluated:#?}",
        (
            coefficients,
            x,
            #[cfg(feature = "precision")]
            order,
        ),
    );
}
