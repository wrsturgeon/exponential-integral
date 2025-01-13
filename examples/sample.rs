//! Evaluate $\text{Ei}$ at a point.

#![expect(
    unused_crate_dependencies,
    reason = "examples won't necessarily use each dev-dependency"
)]
#![expect(clippy::print_stdout, reason = "executable, not a library")]
#![expect(
    clippy::as_conversions,
    clippy::cast_precision_loss,
    reason = "not worth sacrificing clarity to eliminate a one-in-a-million scenario for a one-off example"
)]

use {
    exponential_integral::Ei,
    quickcheck::{Arbitrary, Gen},
    sigma_types::{Finite, NonZero},
};

/// Generate a value within a range, not inclusive.
#[inline]
#[expect(clippy::single_call_fn, reason = "`loop` and `return` semantics")]
fn in_range(lo: f64, hi: f64, g: &mut Gen) -> NonZero<Finite<f64>> {
    loop {
        let u: usize = Arbitrary::arbitrary(g);
        let on_unit = (u as f64) / (usize::MAX as f64);
        let f = on_unit.mul_add(hi - lo, lo);
        if let Some(checked) = Finite::try_new(f).and_then(NonZero::try_new) {
            return checked;
        }
    }
}

fn main() {
    let mut g = Gen::new(256);
    let x = in_range(-16_f64, 16_f64, &mut g);
    println!("x = {x}");
    let ei = Ei(
        x,
        #[cfg(feature = "precision")]
        Arbitrary::arbitrary(&mut g),
    );
    match ei {
        Ok(ok) => println!("Ei({x}) = {ok}"),
        Err(e) => println!("Ei({x}) = [ERROR: {e}]"),
    }
}
