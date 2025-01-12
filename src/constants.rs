//! Relevant internal constants. Not user-facing.

#![expect(
    clippy::excessive_precision,
    clippy::unreadable_literal,
    reason = "copy & paste"
)]

/// Known sizes of constant arrays.
pub(crate) mod size {
    /// AE11
    pub(crate) const AE11: usize = 39;
}

/// AE11
pub(crate) const AE11: [f64; size::AE11] = [
    0.121503239716065790,
    -0.065088778513550150,
    0.004897651357459670,
    -0.000649237843027216,
    0.000093840434587471,
    0.000000420236380882,
    -0.000008113374735904,
    0.000002804247688663,
    0.000000056487164441,
    -0.000000344809174450,
    0.000000058209273578,
    0.000000038711426349,
    -0.000000012453235014,
    -0.000000005118504888,
    0.000000002148771527,
    0.000000000868459898,
    -0.000000000343650105,
    -0.000000000179796603,
    0.000000000047442060,
    0.000000000040423282,
    -0.000000000003543928,
    -0.000000000008853444,
    -0.000000000000960151,
    0.000000000001692921,
    0.000000000000607990,
    -0.000000000000224338,
    -0.000000000000200327,
    -0.000000000000006246,
    0.000000000000045571,
    0.000000000000016383,
    -0.000000000000005561,
    -0.000000000000006074,
    -0.000000000000000862,
    0.000000000000001223,
    0.000000000000000716,
    -0.000000000000000024,
    -0.000000000000000201,
    -0.000000000000000082,
    0.000000000000000017,
];

/*
pub(crate) const AE11_F: &[Finite<f64>; size::AE11] = {
    let ptr: *const [f64; size::AE11] = &AE11;
    let cast: *const [Finite<f64>; size::AE11] = ptr.cast();
    // SAFETY:
    // `Finite<..>` is `repr(transparent)`
    unsafe { &*cast }
};
*/

/// I'd guess that this is the maximum (average?) error between adjacent `f64` values.
pub(crate) const GSL_DBL_EPSILON: f64 = 2.220_446_049_250_313_1e-16;

// pub(crate) const XMAXT: f64 = 708.396_418_532_264_08;

/// No original C code: equal to `-XMAX`.
/// See `XMAX` for its original C code.
pub(crate) const NXMAX: f64 = -XMAX;

/// # Original C code
/// ```c
/// const double XMAX = XMAXT - f64::ln(XMAXT);
/// ```
pub(crate) const XMAX: f64 = 701.833_414_682_1; // XMAXT - f64::ln(XMAXT);
