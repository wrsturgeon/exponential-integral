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

pub mod neg {
    //! Inputs less than 0.

    use {
        crate::{Approx, constants, implementation::neg, pos},
        core::fmt,
        sigma_types::{Finite, Negative},
    };

    /// Argument too large (negative): minimum is `constants::NXMAX`, just under -710.
    #[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
    pub struct HugeArgument(pub(crate) Negative<Finite<f64>>);

    impl fmt::Display for HugeArgument {
        #[inline]
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            let Self(ref arg) = *self;
            write!(
                f,
                "Argument too large (negative): minimum is {}, but {arg} was supplied",
                constants::NXMAX,
            )
        }
    }

    /// E1 on inputs less than 0.
    /// # Errors
    /// If `x` is so large that floating-point operations will fail down the line (absolute value of just over 710).
    #[inline]
    pub fn E1(
        x: Negative<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Result<Approx, HugeArgument> {
        neg::E1(
            x,
            #[cfg(feature = "precision")]
            max_precision,
        )
    }

    /// Ei on inputs less than 0.
    /// # Errors
    /// If `x` is so large that floating-point operations will fail down the line (absolute value of just over 710).
    #[inline(always)]
    pub fn Ei(
        x: Negative<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Result<Approx, HugeArgument> {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        pos::E1(
            -x,
            #[cfg(feature = "precision")]
            max_precision,
        )
        .map(|mut approx| {
            approx.value = -approx.value;
            approx
        })
        .map_err(|pos::HugeArgument(arg)| HugeArgument(-arg))
    }
}

pub mod pos {
    //! Inputs greater than 0.

    use {
        crate::{Approx, constants, implementation::pos, neg},
        core::fmt,
        sigma_types::{Finite, Positive},
    };

    /// Argument too large (positive): maximum is `constants::XMAX`, just over 710.
    #[derive(Clone, Copy, Debug, PartialEq, PartialOrd)]
    pub struct HugeArgument(pub(crate) Positive<Finite<f64>>);

    impl fmt::Display for HugeArgument {
        #[inline]
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            let Self(ref arg) = *self;
            write!(
                f,
                "Argument too large (positive): maximum is {}, but {arg} was supplied",
                constants::XMAX,
            )
        }
    }

    /// E1 on inputs less than 0.
    /// # Errors
    /// If `x` is so large that floating-point operations will fail down the line (absolute value of just over 710).
    #[inline]
    pub fn E1(
        x: Positive<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Result<Approx, HugeArgument> {
        pos::E1(
            x,
            #[cfg(feature = "precision")]
            max_precision,
        )
    }

    /// Ei on inputs less than 0.
    /// # Errors
    /// If `x` is so large that floating-point operations will fail down the line (absolute value of just over 710).
    #[inline(always)]
    pub fn Ei(
        x: Positive<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Result<Approx, HugeArgument> {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        neg::E1(
            -x,
            #[cfg(feature = "precision")]
            max_precision,
        )
        .map(|mut approx| {
            approx.value = -approx.value;
            approx
        })
        .map_err(|neg::HugeArgument(arg)| HugeArgument(-arg))
    }
}

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
            Self::ArgumentTooNegative(arg) => fmt::Display::fmt(&neg::HugeArgument(arg), f),
            Self::ArgumentTooPositive(arg) => fmt::Display::fmt(&pos::HugeArgument(arg), f),
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
/// If `x` is so large that floating-point operations will fail down the line (absolute value of just over 710).
#[inline]
pub fn E1(
    x: NonZero<Finite<f64>>,
    #[cfg(feature = "precision")] max_precision: usize,
) -> Result<Approx, Error> {
    implementation::E1(
        x,
        #[cfg(feature = "precision")]
        max_precision,
    )
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
/// If `x` is so large that floating-point operations will fail down the line (absolute value of just over 710).
#[inline(always)]
pub fn Ei(
    x: NonZero<Finite<f64>>,
    #[cfg(feature = "precision")] max_precision: usize,
) -> Result<Approx, Error> {
    #![expect(
        clippy::arithmetic_side_effects,
        reason = "property-based testing ensures this never happens"
    )]

    E1(
        -x,
        #[cfg(feature = "precision")]
        max_precision,
    )
    .map(|mut approx| {
        approx.value = -approx.value;
        approx
    })
}
