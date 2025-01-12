//! Behind the curtain: actual implementations. May change (but almost surely won't).

/// Specialized approximations to be used on disjoint intervals of the domain,
/// instead of a one-size-fits-all approach.
pub(crate) mod piecewise {
    #![cfg_attr(
        not(test),
        expect(clippy::single_call_fn, reason = "disjoint, so that's kinda the point")
    )]

    use {
        crate::{Approx, chebyshev, constants},
        sigma_types::{Finite, NonNegative, NonZero, One, Sorted, less_than::usize::LessThan},
    };

    /// Between -4 and -1.
    /// # Original C code
    /// ```c
    /// const double ln_term = -log(fabs(x));
    /// gsl_sf_result result_c;
    /// cheb_eval_e(&E11_cs, (2.0*x+5.0)/3.0, &result_c);
    /// result->val  = (ln_term + result_c.val);
    /// result->err  = (result_c.err + GSL_DBL_EPSILON * fabs(ln_term));
    /// result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    /// return GSL_SUCCESS;
    /// ```
    #[inline]
    pub(crate) fn le_neg_1(x: NonZero<Finite<f64>>) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        let abs = Finite::new(x.abs());
        let ln = Finite::new(abs.ln());

        let cheb = chebyshev::eval(
            Sorted::new([Finite::new(-1_f64), Finite::new(1_f64)]),
            Finite::all(&constants::E11),
            LessThan::new(18),
            ((Finite::new(2_f64) * *x) + Finite::new(5_f64)) / Finite::new(3_f64),
        );

        let value = ln + cheb.value;
        let epsilon = NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON));
        let init_err = cheb.error + epsilon * NonNegative::new(Finite::new(ln.abs()));
        let addl_err = NonNegative::new(Finite::new(2_f64))
            * epsilon
            * NonNegative::new(Finite::new(value.abs()));

        Approx {
            value,
            error: init_err + addl_err,
        }
    }

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
    pub(crate) fn le_neg_10(x: NonZero<Finite<f64>>) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        let s: Finite<f64> = (Finite::<f64>::ONE / *x) * (-*x).map(libm::exp);

        let cheb = chebyshev::eval(
            Sorted::new([Finite::new(-1_f64), Finite::new(1_f64)]),
            Finite::all(&constants::AE11),
            LessThan::new(38),
            (Finite::new(20_f64) / *x) + One::ONE,
        );

        let value = s * (Finite::ONE + cheb.value);
        let abs_x: NonNegative<Finite<f64>> = x.map(|f| f.map(f64::abs));
        let abs_value: NonNegative<Finite<f64>> = NonNegative::new(value.map(f64::abs));
        let epsilon = NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON));
        let two = NonNegative::new(Finite::new(2_f64));
        let addl_error = two * epsilon * (abs_x + One::ONE) * abs_value;
        Approx {
            value,
            error: NonNegative::new(s * cheb.error.get() + addl_error.get()),
        }
    }

    /// Between -10 and -4.
    /// # Original C code
    /// ```c
    /// const double s = 1.0/x * ( 0 ? 1.0 : exp(-x) );
    /// gsl_sf_result result_c;
    /// cheb_eval_e(&AE12_cs, (40.0/x+7.0)/3.0, &result_c);
    /// result->val  = s * (1.0 + result_c.val);
    /// result->err  = s * result_c.err;
    /// result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    /// return GSL_SUCCESS;
    /// ```
    #[inline]
    pub(crate) fn le_neg_4(x: NonZero<Finite<f64>>) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        let s: Finite<f64> = (Finite::<f64>::ONE / *x) * (-*x).map(libm::exp);

        let cheb = chebyshev::eval(
            Sorted::new([Finite::new(-1_f64), Finite::new(1_f64)]),
            Finite::all(&constants::AE12),
            LessThan::new(24),
            ((Finite::new(40_f64) / *x) + Finite::new(7_f64)) / Finite::new(3_f64),
        );

        let value = s * (Finite::ONE + cheb.value);
        let abs_value: NonNegative<Finite<f64>> = NonNegative::new(value.map(f64::abs));
        let epsilon = NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON));
        let two = NonNegative::new(Finite::new(2_f64));
        let addl_error = two * epsilon * abs_value;
        Approx {
            value,
            error: NonNegative::new(s * cheb.error.get() + addl_error.get()),
        }
    }
}

use {
    crate::{Approx, Error, constants},
    core::{cmp::Ordering, hint::unreachable_unchecked},
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
#[cfg_attr(
    not(test),
    expect(clippy::single_call_fn, reason = "to mirror the C implementation")
)]
pub(crate) fn E1(x: NonZero<Finite<f64>>) -> Result<Approx, Error> {
    match (**x).partial_cmp(&0_f64) {
        // (-\infty, 0)
        Some(Ordering::Less) => match (**x).partial_cmp(&-10_f64) {
            // (-\infty, -10]
            Some(Ordering::Less | Ordering::Equal) => match (**x).partial_cmp(&constants::NXMAX) {
                // (-XMAX, -10]
                Some(Ordering::Greater) => Ok(piecewise::le_neg_10(x)),
                // (-\infty, -XMAX]
                Some(Ordering::Less | Ordering::Equal) => {
                    Err(Error::ArgumentTooNegative(x.get().get()))
                }
                // SAFETY:
                // absurd case, since `x` is finite
                None => unsafe { unreachable_unchecked() },
            },
            // (-10, 0)
            Some(Ordering::Greater) => match (**x).partial_cmp(&-4_f64) {
                // (-10, -4]
                Some(Ordering::Less | Ordering::Equal) => Ok(piecewise::le_neg_4(x)),
                // (-4, 0)
                Some(Ordering::Greater) => match (**x).partial_cmp(&-1_f64) {
                    // (-4, -1]
                    Some(Ordering::Less | Ordering::Equal) => Ok(piecewise::le_neg_1(x)),
                    // (-1, 0)
                    Some(Ordering::Greater) => todo!(),
                    // SAFETY:
                    // absurd case, since `x` is finite
                    None => unsafe { unreachable_unchecked() },
                },
                // SAFETY:
                // absurd case, since `x` is finite
                None => unsafe { unreachable_unchecked() },
            },
            // SAFETY:
            // absurd case, since `x` is finite
            None => unsafe { unreachable_unchecked() },
        },
        // (0, +\infty)
        Some(Ordering::Greater) => todo!(),
        // SAFETY:
        // absurd case, since `x` is finite and nonzero
        Some(Ordering::Equal) | None => unsafe { unreachable_unchecked() },
    }

    /*
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
