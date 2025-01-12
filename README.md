The exponential integral, often written $\text{Ei}$,
equal to the the integral of an exponentiated input over the input itself:
$\text{Ei}(t) = \int_{-\infty}^{t} \frac{ e^{u} }{ u } \text{d}u$

Inspired by [GSL's implementation](https://github.com/ampl/gsl/blob/ff49e28bdffb893a1c0f6e3eff151296e0e71f82/specfunc/expint.c#L8).
