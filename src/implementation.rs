//! Behind the curtain: actual implementations. May change (but almost surely won't).

pub(crate) mod neg {
    //! E1 for inputs less than 0.

    use {
        crate::{Approx, constants, implementation::piecewise, neg::HugeArgument},
        core::{cmp::Ordering, hint::unreachable_unchecked},
        sigma_types::{Finite, Negative},
    };

    /// See `implementation::E1` for the original C code,
    /// since the original code isn't partitioned by sign.
    /// # Errors
    /// If `x` is so large that floating-point operations will fail down the line (absolute value of just over 710).
    #[inline]
    pub(crate) fn E1(
        x: Negative<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Result<Approx, HugeArgument> {
        match (**x).partial_cmp(&-10_f64) {
            // = -10
            Some(Ordering::Equal) => Ok(piecewise::le_neg_10(
                x,
                #[cfg(feature = "precision")]
                max_precision,
            )),
            // (-\infty, -10)
            Some(Ordering::Less) => match (**x).partial_cmp(&constants::NXMAX) {
                // (-XMAX, -10]
                Some(Ordering::Greater) => Ok(piecewise::le_neg_10(
                    x,
                    #[cfg(feature = "precision")]
                    max_precision,
                )),
                // (-\infty, -XMAX]
                Some(Ordering::Less | Ordering::Equal) => Err(HugeArgument(x)),
                // SAFETY:
                // absurd case: `x` is finite
                None => unsafe { unreachable_unchecked() },
            },
            // (-10, 0)
            Some(Ordering::Greater) => Ok(match (**x).partial_cmp(&-4_f64) {
                // (-10, -4]
                Some(Ordering::Less | Ordering::Equal) => piecewise::le_neg_4(
                    x,
                    #[cfg(feature = "precision")]
                    max_precision,
                ),
                // (-4, 0)
                Some(Ordering::Greater) => match (**x).partial_cmp(&-1_f64) {
                    // (-4, -1]
                    Some(Ordering::Less | Ordering::Equal) => piecewise::le_neg_1(
                        x,
                        #[cfg(feature = "precision")]
                        max_precision,
                    ),
                    // (-1, 0)
                    Some(Ordering::Greater) => piecewise::le_pos_1(
                        x.also(),
                        #[cfg(feature = "precision")]
                        max_precision,
                    ),
                    // SAFETY:
                    // absurd case: `x` is finite
                    None => unsafe { unreachable_unchecked() },
                },
                // SAFETY:
                // absurd case: `x` is finite
                None => unsafe { unreachable_unchecked() },
            }),
            // SAFETY:
            // absurd case: `x` is finite
            None => unsafe { unreachable_unchecked() },
        }
    }
}

/// Specialized approximations to be used on disjoint intervals of the domain,
/// instead of a one-size-fits-all approach.
pub(crate) mod piecewise {
    #![cfg_attr(
        not(test),
        expect(clippy::single_call_fn, reason = "disjoint, so that's kinda the point")
    )]

    use {
        crate::{Approx, chebyshev, constants},
        sigma_types::{Finite, Negative, NonZero, One as _, Positive},
    };

    #[cfg(feature = "error")]
    use sigma_types::NonNegative;

    #[cfg(feature = "precision")]
    use sigma_types::usize::LessThan;

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
    pub(crate) fn le_neg_1(
        x: Negative<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        let abs = Finite::new(x.abs());
        let ln = Finite::new(abs.ln());
        let nln = -ln;

        let cheb = chebyshev::eval(
            Finite::all(&constants::E11),
            ((Finite::new(2_f64) * *x) + Finite::new(5_f64)) / Finite::new(3_f64),
            #[cfg(feature = "precision")]
            LessThan::new(max_precision.min(const { constants::size::E11 - 1 })),
        );

        let value = nln + cheb.value;
        #[cfg(feature = "error")]
        let epsilon = NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON));
        #[cfg(feature = "error")]
        let init_err = cheb.error + epsilon * NonNegative::new(Finite::new(nln.abs()));
        #[cfg(feature = "error")]
        let addl_err = NonNegative::new(Finite::new(2_f64))
            * epsilon
            * NonNegative::new(Finite::new(value.abs()));

        Approx {
            value,
            #[cfg(feature = "error")]
            error: init_err + addl_err,
        }
    }

    /// Between the minimum input (around -710) and -10.
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
    pub(crate) fn le_neg_10(
        x: Negative<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        let s: Finite<f64> = (Finite::<f64>::ONE / *x) * (-*x).map(libm::exp);

        let cheb = chebyshev::eval(
            Finite::all(&constants::AE11),
            (Finite::new(20_f64) / *x) + Finite::<f64>::ONE,
            #[cfg(feature = "precision")]
            LessThan::new(max_precision.min(const { constants::size::AE11 - 1 })),
        );

        let value = s * (Finite::<f64>::ONE + cheb.value);
        #[cfg(feature = "error")]
        let init_err = s * *cheb.error;
        #[cfg(feature = "error")]
        let addl_err = {
            let abs_x: NonNegative<Finite<f64>> = x.map(|f| f.map(f64::abs));
            let abs_value: NonNegative<Finite<f64>> = NonNegative::new(value.map(f64::abs));
            let epsilon = NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON));
            let two = NonNegative::new(Finite::new(2_f64));
            two * epsilon * (abs_x + NonNegative::<Finite<f64>>::ONE) * abs_value
        };
        Approx {
            value,
            #[cfg(feature = "error")]
            error: NonNegative::new(init_err + addl_err.get()),
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
    pub(crate) fn le_neg_4(
        x: Negative<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        let s: Finite<f64> = (Finite::<f64>::ONE / *x) * (-*x).map(libm::exp);

        let cheb = chebyshev::eval(
            Finite::all(&constants::AE12),
            ((Finite::new(40_f64) / *x) + Finite::new(7_f64)) / Finite::new(3_f64),
            #[cfg(feature = "precision")]
            LessThan::new(max_precision.min(const { constants::size::AE12 - 1 })),
        );

        let value = s * (Finite::<f64>::ONE + cheb.value);
        #[cfg(feature = "error")]
        let init_err = s * *cheb.error;
        #[cfg(feature = "error")]
        let addl_err = {
            let abs_value: NonNegative<Finite<f64>> = NonNegative::new(value.map(f64::abs));
            let two = NonNegative::new(Finite::new(2_f64));
            let epsilon = NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON));
            two * epsilon * abs_value
        };
        Approx {
            value,
            #[cfg(feature = "error")]
            error: NonNegative::new(init_err + addl_err.get()),
        }
    }

    /// Between -1 and +1.
    /// # Original C code
    /// ```c
    /// const double ln_term = -log(fabs(x));
    /// gsl_sf_result result_c;
    /// cheb_eval_e(&E12_cs, x, &result_c);
    /// result->val  = (ln_term - 0.6875 + x + result_c.val);
    /// result->err  = (result_c.err + GSL_DBL_EPSILON * fabs(ln_term));
    /// result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    /// return GSL_SUCCESS;
    /// ```
    #[inline]
    pub(crate) fn le_pos_1(
        x: NonZero<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        let abs = Finite::new(x.abs());
        let ln = Finite::new(abs.ln());
        let nln = -ln;

        let cheb = chebyshev::eval(
            Finite::all(&constants::E12),
            *x,
            #[cfg(feature = "precision")]
            LessThan::new(max_precision.min(const { constants::size::E12 - 1 })),
        );

        let value = nln - Finite::new(0.6875_f64) + *x + cheb.value;
        #[cfg(feature = "error")]
        let epsilon = NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON));
        #[cfg(feature = "error")]
        let init_err = cheb.error + epsilon * NonNegative::new(Finite::new(nln.abs()));
        #[cfg(feature = "error")]
        let addl_err = NonNegative::new(Finite::new(2_f64))
            * epsilon
            * NonNegative::new(Finite::new(value.abs()));

        Approx {
            value,
            #[cfg(feature = "error")]
            error: init_err + addl_err,
        }
    }

    /// Between +1 and +4.
    /// # Original C code
    /// ```c
    /// const double s = 1.0/x * exp(-x);
    /// gsl_sf_result result_c;
    /// cheb_eval_e(&AE13_cs, (8.0/x-5.0)/3.0, &result_c);
    /// result->val  = s * (1.0 + result_c.val);
    /// result->err  = s * result_c.err;
    /// result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    /// return GSL_SUCCESS;
    /// ```
    #[inline]
    pub(crate) fn le_pos_4(
        x: Positive<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        let s = (Finite::<f64>::ONE / *x) * (-*x).map(f64::exp);

        let cheb = chebyshev::eval(
            Finite::all(&constants::AE13),
            (Finite::new(8_f64) / *x - Finite::new(5_f64)) / Finite::new(3_f64),
            #[cfg(feature = "precision")]
            LessThan::new(max_precision.min(const { constants::size::AE13 - 1 })),
        );

        let value = s * (Finite::<f64>::ONE + cheb.value);
        #[cfg(feature = "error")]
        let init_err = s * *cheb.error;
        #[cfg(feature = "error")]
        let addl_err = {
            let epsilon = NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON));
            NonNegative::new(Finite::new(2_f64))
                * epsilon
                * NonNegative::new(Finite::new(value.abs()))
        };

        Approx {
            value,
            #[cfg(feature = "error")]
            error: NonNegative::new(init_err + *addl_err),
        }
    }

    /// Between +4 and the maximum input (around 710).
    /// # Original C code
    /// ```c
    /// const double s = 1.0/x * exp(-x);
    /// gsl_sf_result result_c;
    /// cheb_eval_e(&AE14_cs, 8.0/x-1.0, &result_c);
    /// result->val  = s * (1.0 +  result_c.val);
    /// result->err  = s * (GSL_DBL_EPSILON + result_c.err);
    /// result->err += 2.0 * (x + 1.0) * GSL_DBL_EPSILON * fabs(result->val);
    /// if(result->val == 0.0)
    ///   UNDERFLOW_ERROR(result);
    /// else
    ///   return GSL_SUCCESS;
    /// ```
    #[inline]
    pub(crate) fn le_pos_max(
        x: Positive<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Approx {
        #![expect(
            clippy::arithmetic_side_effects,
            reason = "property-based testing ensures this never happens"
        )]

        let s = (Finite::<f64>::ONE / *x) * (-*x).map(f64::exp);

        let cheb = chebyshev::eval(
            Finite::all(&constants::AE14),
            (Finite::new(8_f64) / *x) - Finite::new(1_f64),
            #[cfg(feature = "precision")]
            LessThan::new(max_precision.min(const { constants::size::AE14 - 1 })),
        );

        let value = s * (Finite::<f64>::ONE + cheb.value);
        #[cfg(feature = "error")]
        let epsilon = NonNegative::new(Finite::new(constants::GSL_DBL_EPSILON));
        #[cfg(feature = "error")]
        let init_err = s * *(epsilon + cheb.error);
        #[cfg(feature = "error")]
        let addl_err = {
            let also_x: NonNegative<Finite<f64>> = x.also();
            NonNegative::new(Finite::new(2_f64))
                * (also_x + NonNegative::new(Finite::new(1_f64)))
                * epsilon
                * NonNegative::new(Finite::new(value.abs()))
        };

        Approx {
            value,
            #[cfg(feature = "error")]
            error: NonNegative::new(init_err + *addl_err),
        }
    }
}

pub(crate) mod pos {
    //! E1 for inputs greater than 0.

    use {
        crate::{Approx, constants, implementation::piecewise, pos::HugeArgument},
        core::{cmp::Ordering, hint::unreachable_unchecked},
        sigma_types::{Finite, Positive},
    };

    /// See `implementation::E1` for the original C code,
    /// since the original code isn't partitioned by sign.
    /// # Errors
    /// If `x` is so large that floating-point operations will fail down the line (absolute value of just over 710).
    #[inline]
    pub(crate) fn E1(
        x: Positive<Finite<f64>>,
        #[cfg(feature = "precision")] max_precision: usize,
    ) -> Result<Approx, HugeArgument> {
        match (**x).partial_cmp(&4_f64) {
            // = 4
            Some(Ordering::Equal) => Ok(piecewise::le_pos_4(
                x,
                #[cfg(feature = "precision")]
                max_precision,
            )),
            // (0, +4)
            Some(Ordering::Less) => Ok(match (**x).partial_cmp(&1_f64) {
                // (0, +1]
                Some(Ordering::Less | Ordering::Equal) => piecewise::le_pos_1(
                    x.also(),
                    #[cfg(feature = "precision")]
                    max_precision,
                ),
                // (+1, +\infty]
                Some(Ordering::Greater) => piecewise::le_pos_4(
                    x,
                    #[cfg(feature = "precision")]
                    max_precision,
                ),
                // SAFETY:
                // absurd case: `x` is finite
                None => unsafe { unreachable_unchecked() },
            }),
            // (+4, +\infty)
            Some(Ordering::Greater) => match (**x).partial_cmp(&constants::XMAX) {
                Some(Ordering::Less) => Ok(piecewise::le_pos_max(
                    x,
                    #[cfg(feature = "precision")]
                    max_precision,
                )),
                Some(Ordering::Equal | Ordering::Greater) => Err(HugeArgument(x)),
                // SAFETY:
                // absurd case: `x` is finite
                None => unsafe { unreachable_unchecked() },
            },
            // SAFETY:
            // absurd case: `x` is finite
            None => unsafe { unreachable_unchecked() },
        }
    }
}

use {
    crate::{Approx, Error},
    core::{cmp::Ordering, hint::unreachable_unchecked},
    sigma_types::{Finite, NonZero},
};

/// # Errors
/// If `x` is so large that floating-point operations will fail down the line (absolute value of just over 710).
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
#[cfg_attr(
    not(test),
    expect(clippy::single_call_fn, reason = "to mirror the C implementation")
)]
#[expect(clippy::absolute_paths, reason = "always a collision except full path")]
pub(crate) fn E1(
    x: NonZero<Finite<f64>>,
    #[cfg(feature = "precision")] max_precision: usize,
) -> Result<Approx, Error> {
    match (**x).partial_cmp(&0_f64) {
        // (-\infty, 0)
        Some(Ordering::Less) => neg::E1(
            x.also(),
            #[cfg(feature = "precision")]
            max_precision,
        )
        .map_err(|crate::neg::HugeArgument(arg)| Error::ArgumentTooNegative(arg)),
        // (0, +\infty)
        Some(Ordering::Greater) => pos::E1(
            x.also(),
            #[cfg(feature = "precision")]
            max_precision,
        )
        .map_err(|crate::pos::HugeArgument(arg)| Error::ArgumentTooPositive(arg)),
        // SAFETY:
        // absurd case: `x` is finite and nonzero
        Some(Ordering::Equal) | None => unsafe { unreachable_unchecked() },
    }
}
