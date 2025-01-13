//! Chebyshev series/polynomial approximation.

use {
    crate::Approx,
    sigma_types::{Finite, Zero as _, less_than::usize::LessThan},
};

#[cfg(feature = "error")]
use {crate::constants, sigma_types::NonNegative};

/// Chebyshev series/polynomial approximation.
/// # Original C code
/// ```c
/// struct cheb_series_struct {
///   double * c;   /* coefficients                */
///   int order;    /* order of expansion          */
///   double a;     /* lower interval point        */
///   double b;     /* upper interval point        */
///   int order_sp; /* effective single precision order */
/// };
/// typedef struct cheb_series_struct cheb_series;
///
/// // ...
///
/// static inline int
/// cheb_eval_e(const cheb_series * cs,
///             const double x,
///             gsl_sf_result * result)
/// {
///   int j;
///   double d  = 0.0;
///   double dd = 0.0;
///
///   double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
///   double two_x = 2.0 * y;
///
///   double e = 0.0;
///
///   for(j = cs->order; j>=1; j--) {
///     double temp = d;
///     d = two_x*d - dd + cs->c[j];
///     e += fabs(two_x*temp) + fabs(dd) + fabs(cs->c[j]);
///     dd = temp;
///   }
///
///   {
///     double temp = d;
///     d = y*d - dd + 0.5 * cs->c[0];
///     e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
///   }
///
///   result->val = d;
///   result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);
///
///   return GSL_SUCCESS;
/// }
/// ```
#[inline]
#[must_use]
pub fn eval<const N_COEFFICIENTS: usize>(
    coefficients: &[Finite<f64>; N_COEFFICIENTS],
    x: Finite<f64>,
    #[cfg(feature = "precision")] order: LessThan<{ N_COEFFICIENTS }>,
) -> Approx {
    #![expect(
        clippy::arithmetic_side_effects,
        reason = "property-based testing ensures this never happens"
    )]

    debug_assert!(
        N_COEFFICIENTS > 0,
        "Chebyshev series without any coefficients",
    );

    let two_x: Finite<f64> = Finite::new(2_f64) * x;

    #[cfg(feature = "error")]
    let mut e = NonNegative::<Finite<f64>>::ZERO;

    let mut d = Finite::<f64>::ZERO;
    let mut dd = Finite::<f64>::ZERO;

    {
        let mut j: LessThan<{ N_COEFFICIENTS }> = {
            #[cfg(feature = "precision")]
            {
                order
            }
            #[cfg(not(feature = "precision"))]
            {
                LessThan::new(const { N_COEFFICIENTS - 1 })
            }
        };
        while *j >= 1 {
            // SAFETY:
            // See the `debug_assert` above.
            let coefficient = *unsafe { coefficients.get_unchecked(*j) };
            let tmp = d;
            d = ((two_x * d) - dd) + coefficient;
            #[cfg(feature = "error")]
            {
                e += NonNegative::new((two_x * tmp).map(f64::abs))
                    + NonNegative::new(dd.map(f64::abs))
                    + NonNegative::new(coefficient.map(f64::abs));
            }
            dd = tmp;

            j.map_mut(|u| *u -= 1);
        }
    }

    {
        #[cfg(feature = "error")]
        let tmp = d;
        // SAFETY:
        // Sigma types ensure validity.
        let coefficient = *unsafe { coefficients.get_unchecked(0) };
        let half_coefficient = coefficient.map(|c| 0.5_f64 * c);
        d = x * d - dd + half_coefficient;
        #[cfg(feature = "error")]
        {
            e += NonNegative::new((x * tmp).map(f64::abs))
                + NonNegative::new(dd.map(f64::abs))
                + NonNegative::new(half_coefficient.map(f64::abs));
        }
    }

    #[cfg(feature = "error")]
    // SAFETY:
    // See `debug_assert`s above.
    let last_coefficient = *unsafe { coefficients.get_unchecked(const { N_COEFFICIENTS - 1 }) };

    Approx {
        value: d,
        #[cfg(feature = "error")]
        error: NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON)) * e
            + NonNegative::new(last_coefficient.map(f64::abs)),
    }
}

/// Compile-time-compatible minimum of two large unsigned integers.
#[inline]
#[cfg_attr(not(test), expect(dead_code, reason = "TODO: REMOVE"))]
#[cfg_attr(test, expect(clippy::single_call_fn, reason = "TODO: REMOVE"))]
pub(crate) const fn min(a: usize, b: usize) -> usize {
    if a.checked_sub(b).is_some() { b } else { a }
}
