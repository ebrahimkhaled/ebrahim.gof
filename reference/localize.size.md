# How large is each part of a misfit?

Gives the size of the four parts of
[`localize.update`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.update.md)
on the risk scale, for frozen predictions on validation data. Four
nested updates are fitted: the intercept update, the logistic
recalibration, the flexible recalibration `ns(eta, 4)`, and a revision
that adds the main effects and the `COV` bases as one logistic model.
The size of a part is the mean absolute change in predicted risk that
its update makes beyond the simpler one, in percentage points.

## Usage

``` r
localize.size(y, p, X, R = 200L, seed = NULL)
```

## Arguments

- y, p, X:

  as in
  [`localize.update`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.update.md).

- R:

  number of chance draws for each part.

- seed:

  optional integer, passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the chance
  draws. The random number stream of the session is restored on exit.

## Value

An object of class `"gof_localize_size"`: a list with `parts` (a data
frame with `part`, `size`, `chance`, `excess` and the noise interval
`lo`, `hi`, in percentage points), `calibration` (the calibration
intercept and slope with Wald 95% intervals) and `ici` (the integrated
calibration index of Austin and Steyerberg, in percentage points).

## Details

An update fitted to noise also moves the predictions, so each size comes
with its chance level: the mean size over `R` outcome vectors drawn from
the simpler update, under which the part is absent. `excess` is the size
minus the chance level. The noise interval is the size minus the 97.5%
and 2.5% points of the chance draws; under no misfit it contains 0 with
probability of about 0.95. It is not a confidence interval for the true
size: it only shows whether a size can be told from chance.

`excess` is a cautious lower measure, not an estimate of the true size.
Every size is fitted with noise, so at moderate sample sizes the excess
understates the true size, most of all for `COV`. A size that cannot be
told from its chance level does not show that a part is absent;
[`localize.update`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.update.md)
carries that evidence.

## References

Ebrahim, E. K., Hussein, O. A. E.-A. and El-Kotory, A. (2026). Which
update does a clinical risk model need? Localization of misfit in
logistic risk models. Submitted. Reproduction materials:
[doi:10.5281/zenodo.24025418](https://doi.org/10.5281/zenodo.24025418)

Austin, P. C. and Steyerberg, E. W. (2019). The Integrated Calibration
Index (ICI) and related metrics for quantifying the calibration of
logistic regression models. *Statistics in Medicine*, 38(21), 4051–4065.
[doi:10.1002/sim.8281](https://doi.org/10.1002/sim.8281)

## See also

[`localize.update`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.update.md)
for the verdict.

## Examples

``` r
set.seed(1)
n <- 500
X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)
y <- rbinom(n, 1, plogis(-0.7 + qlogis(p)))   # overall risk too high
localize.size(y, p, X, R = 20, seed = 1)      # R = 20 to keep the example fast; use 200
#> 
#> Size of each part of the misfit, percentage points of risk
#> 
#>            size chance excess noise.lo noise.hi
#> INTERCEPT 12.77   1.65  11.11     9.59    12.37
#> SLOPE      1.21   1.15   0.07    -1.74     1.13
#> LINK       1.93   2.04  -0.12    -1.85     1.03
#> COV        4.48   3.58   0.90    -0.62     2.41
#> 
#> calibration intercept -0.697 (-0.913 to -0.482), slope 1.107 (0.849 to 1.366)
#> integrated calibration index 12.94
#> 
```
