//! Approximate values alongside estimates of their own approximation errors.

use {
    crate::{Approx, Error, implementation},
    sigma_types::{Finite, NonZero},
};

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

    E1(-x).map(|Approx { error, value }| Approx {
        error,
        value: -value,
    })
}
