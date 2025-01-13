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

use {
    core::fmt,
    sigma_types::{Finite, Negative, NonZero, Positive},
};

#[cfg(feature = "error")]
use sigma_types::NonNegative;

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
    #[cfg(feature = "error")]
    pub error: NonNegative<Finite<f64>>,
    /// Approximate value.
    pub value: Finite<f64>,
}

impl fmt::Display for Approx {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self {
            #[cfg(feature = "error")]
            ref error,
            ref value,
        } = *self;
        #[cfg(feature = "error")]
        {
            write!(f, "{value} +/- {error}")
        }
        #[cfg(not(feature = "error"))]
        {
            write!(f, "{value}")
        }
    }
}

/// An approximate value alongside an estimate of its own approximation error.
#[non_exhaustive]
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
pub enum Error {
    /// Argument was less than the safe minimum.
    ArgumentTooNegative(Negative<Finite<f64>>),
    /// Argument was less than the safe maximum.
    ArgumentTooPositive(Positive<Finite<f64>>),
}

impl fmt::Display for Error {
    #[inline]
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::ArgumentTooNegative(arg) => write!(
                f,
                "Argument too large (negative): minimum is {}, but {arg} was supplied",
                constants::NXMAX,
            ),
            Self::ArgumentTooPositive(arg) => write!(
                f,
                "Argument too large (positive): maximum is {}, but {arg} was supplied",
                constants::XMAX,
            ),
        }
    }
}

/// # Original C code
/// ```c
/// int gsl_sf_expint_E1_e(const double x, gsl_sf_result * result)
/// {
///   return expint_E1_impl(x, result, 0);
/// }
/// ```
///
/// # Errors
/// See `Error`.
#[inline]
pub fn E1(x: NonZero<Finite<f64>>) -> Result<Approx, Error> {
    implementation::E1(x)
}

/// # Original C code
/// ```c
/// int gsl_sf_expint_Ei_e(const double x, gsl_sf_result * result)
/// {
///   /* CHECK_POINTER(result) */
///
///   {
///     int status = gsl_sf_expint_E1_e(-x, result);
///     result->val = -result->val;
///     return status;
///   }
/// }
/// ```
///
/// # Errors
/// See `Error`.
#[inline(always)]
pub fn Ei(x: NonZero<Finite<f64>>) -> Result<Approx, Error> {
    #![expect(
        clippy::arithmetic_side_effects,
        reason = "property-based testing ensures this never happens"
    )]

    E1(-x).map(|mut approx| {
        approx.value = -approx.value;
        approx
    })
}
