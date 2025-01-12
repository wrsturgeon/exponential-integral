//! Evaluate $\text{Ei}$ at a point.

#![expect(
    unused_crate_dependencies,
    reason = "examples won't necessarily use each dev-dependency"
)]
#![expect(clippy::print_stdout, reason = "executable, not a library")]

use {
    core::cmp::Ordering,
    exponential_integral::with_error::Ei,
    quickcheck::{Arbitrary, Gen},
    sigma_types::{Finite, NonZero},
};

/// Generate a value within a range, not inclusive.
#[inline]
#[expect(clippy::single_call_fn, reason = "`loop` and `return` semantics")]
fn in_range<T: Arbitrary + PartialOrd>(lo: &T, hi: &T, g: &mut Gen) -> T {
    loop {
        let t = T::arbitrary(g);
        let Some(Ordering::Greater) = t.partial_cmp(lo) else {
            continue;
        };
        let Some(Ordering::Less) = t.partial_cmp(hi) else {
            continue;
        };
        return t;
    }
}

fn main() {
    let x = in_range(
        &NonZero::new(Finite::new(1_f64)),
        &NonZero::new(Finite::new(16_f64)),
        &mut Gen::new(256),
    );
    println!("x = {x}");
    let ei = Ei(x);
    match ei {
        Ok(ok) => println!("Ei({x}) = {ok}"),
        Err(e) => println!("Ei({x}) = [ERROR: {e}]"),
    }
}
