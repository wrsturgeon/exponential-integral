//! Behind the curtain: actual implementations. May change (but almost surely won't).

/// Specialized approximations to be used on disjoint intervals of the domain,
/// instead of a one-size-fits-all approach.
pub(crate) mod piecewise {
    #![expect(clippy::single_call_fn, reason = "disjoint, so that's kinda the point")]

    use {
        crate::{Approx, chebyshev, constants},
        sigma_types::{Finite, NonNegative, NonZero, One},
    };

    /// Between the minimum input and -10.
    /// # Original C code
    /// ```c
    /// const double s = 1.0/x * ( 0 ? 1.0 : exp(-x) );
    /// gsl_sf_result result_c;
    /// cheb_eval_e(&AE11_cs, 20.0/x+1.0, &result_c);
    /// result->val  = s * (1.0 + result_c.val);
    /// result->err  = s * result_c.err;
    /// result->err += 2.0 * GSL_DBL_EPSILON * (fabs(x) + 1.0) * fabs(result->val);
    /// return GSL_SUCCESS;
    /// ```
    #[inline]
    pub(crate) fn lt_10(x: NonZero<Finite<f64>>) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        let s: Finite<f64> = (Finite::<f64>::ONE / *x) * (-*x).map(libm::exp);

        let ae11 = chebyshev::Series {
            a: Finite::new(-1_f64),
            b: Finite::new(1_f64),
            coefficients: &constants::AE11,
            order: 38,
            // order_sp: 20,
        };
        let Approx {
            value: cheb,
            error: _,
        } = ae11.eval((Finite::new(20_f64) / *x) + One::ONE);

        let value = s * (Finite::ONE + cheb);
        let abs_x: NonNegative<Finite<f64>> = x.map(|f| f.map(f64::abs));
        let abs_value: NonNegative<Finite<f64>> = NonNegative::new(value.map(f64::abs));
        let epsilon = NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON));
        let two = NonNegative::new(Finite::new(2_f64));
        Approx {
            value,
            error: two * epsilon * (abs_x + One::ONE) * abs_value,
        }
    }
}

use {
    crate::{Approx, Error, constants},
    sigma_types::{Finite, NonZero},
};

/// # Original C code
/// Note that `scale` is pinned to `0`.
/// ```c
/// /* implementation for E1, allowing for scaling by exp(x) */
/// static
/// int expint_E1_impl(const double x, gsl_sf_result * result, const int scale)
/// {
///   const double xmaxt = -GSL_LOG_DBL_MIN;      /* XMAXT = -LOG (R1MACH(1)) */
///   const double xmax  = xmaxt - log(xmaxt);    /* XMAX = XMAXT - LOG(XMAXT) */
///
///   /* CHECK_POINTER(result) */
///
///   if(x < -xmax && !scale) {
///       OVERFLOW_ERROR(result);
///   }
///   else if(x <= -10.0) {
///     const double s = 1.0/x * ( scale ? 1.0 : exp(-x) );
///     gsl_sf_result result_c;
///     cheb_eval_e(&AE11_cs, 20.0/x+1.0, &result_c);
///     result->val  = s * (1.0 + result_c.val);
///     result->err  = s * result_c.err;
///     result->err += 2.0 * GSL_DBL_EPSILON * (fabs(x) + 1.0) * fabs(result->val);
///     return GSL_SUCCESS;
///   }
///   else if(x <= -4.0) {
///     const double s = 1.0/x * ( scale ? 1.0 : exp(-x) );
///     gsl_sf_result result_c;
///     cheb_eval_e(&AE12_cs, (40.0/x+7.0)/3.0, &result_c);
///     result->val  = s * (1.0 + result_c.val);
///     result->err  = s * result_c.err;
///     result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
///     return GSL_SUCCESS;
///   }
///   else if(x <= -1.0) {
///     const double ln_term = -log(fabs(x));
///     const double scale_factor = ( scale ? exp(x) : 1.0 );
///     gsl_sf_result result_c;
///     cheb_eval_e(&E11_cs, (2.0*x+5.0)/3.0, &result_c);
///     result->val  = scale_factor * (ln_term + result_c.val);
///     result->err  = scale_factor * (result_c.err + GSL_DBL_EPSILON * fabs(ln_term));
///     result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
///     return GSL_SUCCESS;
///   }
///   else if(x == 0.0) {
///     DOMAIN_ERROR(result);
///   }
///   else if(x <= 1.0) {
///     const double ln_term = -log(fabs(x));
///     const double scale_factor = ( scale ? exp(x) : 1.0 );
///     gsl_sf_result result_c;
///     cheb_eval_e(&E12_cs, x, &result_c);
///     result->val  = scale_factor * (ln_term - 0.6875 + x + result_c.val);
///     result->err  = scale_factor * (result_c.err + GSL_DBL_EPSILON * fabs(ln_term));
///     result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
///     return GSL_SUCCESS;
///   }
///   else if(x <= 4.0) {
///     const double s = 1.0/x * ( scale ? 1.0 : exp(-x) );
///     gsl_sf_result result_c;
///     cheb_eval_e(&AE13_cs, (8.0/x-5.0)/3.0, &result_c);
///     result->val  = s * (1.0 + result_c.val);
///     result->err  = s * result_c.err;
///     result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
///     return GSL_SUCCESS;
///   }
///   else if(x <= xmax || scale) {
///     const double s = 1.0/x * ( scale ? 1.0 : exp(-x) );
///     gsl_sf_result result_c;
///     cheb_eval_e(&AE14_cs, 8.0/x-1.0, &result_c);
///     result->val  = s * (1.0 +  result_c.val);
///     result->err  = s * (GSL_DBL_EPSILON + result_c.err);
///     result->err += 2.0 * (x + 1.0) * GSL_DBL_EPSILON * fabs(result->val);
///     if(result->val == 0.0)
///       UNDERFLOW_ERROR(result);
///     else
///       return GSL_SUCCESS;
///   }
///   else {
///     UNDERFLOW_ERROR(result);
///   }
/// }
/// ```
///
/// # Errors
/// See `Error`.
#[inline]
#[expect(clippy::todo, reason = "TODO: REMOVE")]
#[expect(clippy::single_call_fn, reason = "to mirror the C implementation")]
pub(crate) fn E1(x: NonZero<Finite<f64>>) -> Result<Approx, Error> {
    match **x {
        ..constants::NXMAX => Err(Error::ArgumentTooNegative(x.get().get())),
        ..-10. => Ok(piecewise::lt_10(x)),
        _ => todo!(),
    }

    /*
    else if(x <= -4.0) {
      const double s = 1.0/x * ( 0 ? 1.0 : exp(-x) );
      gsl_sf_result result_c;
      cheb_eval_e(&AE12_cs, (40.0/x+7.0)/3.0, &result_c);
      result->val  = s * (1.0 + result_c.val);
      result->err  = s * result_c.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x <= -1.0) {
      const double ln_term = -log(fabs(x));
      const double scale_factor = ( 0 ? exp(x) : 1.0 );
      gsl_sf_result result_c;
      cheb_eval_e(&E11_cs, (2.0*x+5.0)/3.0, &result_c);
      result->val  = scale_factor * (ln_term + result_c.val);
      result->err  = scale_factor * (result_c.err + GSL_DBL_EPSILON * fabs(ln_term));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x == 0.0) {
      DOMAIN_ERROR(result);
    }
    else if(x <= 1.0) {
      const double ln_term = -log(fabs(x));
      const double scale_factor = ( 0 ? exp(x) : 1.0 );
      gsl_sf_result result_c;
      cheb_eval_e(&E12_cs, x, &result_c);
      result->val  = scale_factor * (ln_term - 0.6875 + x + result_c.val);
      result->err  = scale_factor * (result_c.err + GSL_DBL_EPSILON * fabs(ln_term));
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x <= 4.0) {
      const double s = 1.0/x * ( 0 ? 1.0 : exp(-x) );
      gsl_sf_result result_c;
      cheb_eval_e(&AE13_cs, (8.0/x-5.0)/3.0, &result_c);
      result->val  = s * (1.0 + result_c.val);
      result->err  = s * result_c.err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
      return GSL_SUCCESS;
    }
    else if(x <= XMAX || 0) {
      const double s = 1.0/x * ( 0 ? 1.0 : exp(-x) );
      gsl_sf_result result_c;
      cheb_eval_e(&AE14_cs, 8.0/x-1.0, &result_c);
      result->val  = s * (1.0 +  result_c.val);
      result->err  = s * (GSL_DBL_EPSILON + result_c.err);
      result->err += 2.0 * (x + 1.0) * GSL_DBL_EPSILON * fabs(result->val);
      if(result->val == 0.0)
        UNDERFLOW_ERROR(result);
      else
        return GSL_SUCCESS;
    }
    else {
      UNDERFLOW_ERROR(result);
    }
    */
}
