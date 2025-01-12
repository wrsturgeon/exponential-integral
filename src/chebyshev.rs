//! Chebyshev series.

use {
    crate::{Approx, constants},
    sigma_types::{Finite, NonNegative, Zero as _},
};

/// Chebyshev series.
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
/// ```
pub(crate) struct Series<'c, const N: usize> {
    /// Lower bound.
    pub(crate) a: Finite<f64>,
    /// Upper bound.
    pub(crate) b: Finite<f64>,
    /// Coefficients.
    pub(crate) coefficients: &'c [f64; N],
    /// Order of expansion.
    pub(crate) order: usize,
    // pub(crate) order_sp: i32,
}

impl<const N: usize> Series<'_, N> {
    /// # Original C code
    /// ```c
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
    ///   double y2 = 2.0 * y;
    ///
    ///   double e = 0.0;
    ///
    ///   for(j = cs->order; j>=1; j--) {
    ///     double temp = d;
    ///     d = y2*d - dd + cs->c[j];
    ///     e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
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
    #[expect(clippy::many_single_char_names, reason = "original names")]
    pub(crate) fn eval(&self, x: Finite<f64>) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        debug_assert!(N > 0, "Chebyshev series without any coefficients");

        debug_assert!(
            self.order < N,
            "Order of expansion out of bounds: {} >= {N}",
            self.order,
        );

        let y: Finite<f64> = (((Finite::new(2_f64) * x) - self.a) - self.b) / (self.b - self.a);
        let y2: Finite<f64> = Finite::new(2_f64) * y;

        let mut e = NonNegative::<Finite<f64>>::ZERO;

        let mut d = Finite::<f64>::ZERO;
        let mut dd = Finite::<f64>::ZERO;

        {
            let mut j = self.order;
            while j >= 1 {
                // SAFETY:
                // See `debug_assert`s above.
                let coefficient = Finite::new(*unsafe { self.coefficients.get_unchecked(j) });
                let tmp = d;
                d = ((y2 * d) - dd) + coefficient;
                e += NonNegative::new((y2 * tmp).map(f64::abs))
                    + NonNegative::new(dd.map(f64::abs))
                    + NonNegative::new(coefficient.map(f64::abs));
                dd = tmp;

                j -= 1;
            }
        }

        {
            let tmp = d;
            // SAFETY:
            // See `debug_assert`s above.
            let coefficient = Finite::new(*unsafe { self.coefficients.get_unchecked(0) });
            let half_coefficient = coefficient.map(|c| 0.5_f64 * c);
            d = y * d - dd + half_coefficient;
            e += NonNegative::new((y * tmp).map(f64::abs))
                + NonNegative::new(dd.map(f64::abs))
                + NonNegative::new(half_coefficient.map(f64::abs));
        }

        // SAFETY:
        // See `debug_assert`s above.
        let last_coefficient = Finite::new(*unsafe { self.coefficients.get_unchecked(self.order) });

        Approx {
            value: d,
            error: NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON)) * e
                + NonNegative::new(last_coefficient.map(f64::abs)),
        }
    }
}
