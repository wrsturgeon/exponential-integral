//! The exponential integral, often written $\text{Ei}$,
//! equal to the the integral of an exponentiated input over the input itself:
//! $\text{Ei}(t) = \int_{-\infty}^{t} \frac{ e^{u} }{ u } \text{d}u$
//!
//! Inspired by [GSL's implementation](https://github.com/ampl/gsl/blob/ff49e28bdffb893a1c0f6e3eff151296e0e71f82/specfunc/expint.c#L8).

#![no_std]
#![expect(non_snake_case, reason = "Proper mathematical names")]

pub mod chebyshev;
mod constants;
mod implementation;

#[cfg(test)]
mod test;

pub mod with_error;

use {
    core::fmt,
    sigma_types::{Finite, NonNegative, NonZero},
};

/// An approximate value alongside an estimate of its own approximation error.
/// # Original C code
/// ```c
/// struct gsl_sf_result_struct {
///   double val;
///   double err;
/// };
/// typedef struct gsl_sf_result_struct gsl_sf_result;
/// ```
#[expect(clippy::exhaustive_structs, reason = "Simple structure")]
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub struct Approx {
    /// Estimate of the approximation error for `value`.
    pub error: NonNegative<Finite<f64>>,
    /// Approximate value.
    pub value: Finite<f64>,
}

impl fmt::Display for Approx {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self {
            ref error,
            ref value,
        } = *self;
        write!(f, "{value} +/- {error}")
    }
}

/// An approximate value alongside an estimate of its own approximation error.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub enum Error {
    /// Argument was less than the safe minimum.
    ArgumentTooNegative(f64),
}

impl fmt::Display for Error {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::ArgumentTooNegative(arg) => write!(
                f,
                "Argument too negative: minimum is {}, but {arg} was supplied",
                constants::NXMAX,
            ),
        }
    }
}

/// Original C code:
/// ```c
/// #define EVAL_RESULT(fn) \
///    gsl_sf_result result; \
///    int status = fn; \
///    if (status != GSL_SUCCESS) { \
///      GSL_ERROR_VAL(#fn, status, result.val); \
///    } ; \
///    return result.val;
///
/// // ...
///
/// double gsl_sf_expint_Ei(const double x)
/// {
///   EVAL_RESULT(gsl_sf_expint_Ei_e(x, &result));
/// }
/// ```
///
/// # Errors
/// See `Error`.
#[inline(always)]
pub fn Ei(x: NonZero<Finite<f64>>) -> Result<Finite<f64>, Error> {
    with_error::Ei(x).map(|Approx { error: _, value }| value)
}
