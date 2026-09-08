# DeepGOF-1: a pretrained goodness-of-fit test for logistic regression

Tests whether a fitted binomial `glm` is correctly specified, using a
convolutional network that was trained once, offline, on simulated
departures and is shipped frozen with this package. The analyst never
trains anything: the network reads the model's residual map and the
p-value is the rank of the observed score within the analyst's own
parametric bootstrap, so the level does not depend on what the network
learned.

## Usage

``` r
deepgof1(fit, B = 199L, K = 6L)
```

## Arguments

- fit:

  a fitted `glm` with `family = binomial()` and at least two covariates.

- B:

  number of parametric-bootstrap replicates. The p-value lies on a grid
  of `1/(B+1)`, so `B = 199` makes the nominal .05 attainable exactly.

- K:

  grid resolution. Leave at 6: the shipped weights were trained at
  `K = 6` and are not valid at any other resolution.

## Value

An object of class `"deepgof1"`: a list with `statistic` (the observed
score), `p.value`, `B`, `K`, `axes` (the two covariates the grid was
built on), `boot` (the `B` bootstrap scores) and `method`.

## Details

The residual map is a `K` by `K` grid over the empirical ranks of the
two covariates with the largest \\\|\hat\beta_j\| \hat\sigma_j\\; each
cell holds a standardized residual sum, approximately standard normal
under a correct model. Misfit therefore has a location on the map – an
omitted quadratic paints a stripe, an omitted interaction a saddle –
which is what the convolutional statistic reads.

The test is a small-sample instrument. Against the classical partition
tests it gains most at \\n\\ of 50 to 200 and the gain decays as \\n\\
grows; because the grid uses only two covariates, misfit that lives off
those axes is harder for it to see than for covariate-space or smoothing
tests. It is not one of the tests
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
selects: call it directly on the same fitted model and read its p-value
beside the panel. See
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
for the classical battery.

## Reproducibility

A bootstrap refit that fails to converge is scored `+Inf`, so it counts
against rejection – the conservative direction. Set a seed before
calling for a reproducible p-value.

## References

Ebrahim EK (2026). "DeepGOF-1: A Pretrained Convolutional
Goodness-of-Fit Test for Logistic Regression with a Computable
Consistency Certificate." Manuscript under review. Reproduction
materials and frozen weights:
[doi:10.5281/zenodo.22113220](https://doi.org/10.5281/zenodo.22113220)

Besag, J. and Clifford, P. (1989). Generalized Monte Carlo significance
tests. *Biometrika* **76**, 633–642.
[doi:10.1093/biomet/76.4.633](https://doi.org/10.1093/biomet/76.4.633)

## See also

[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md),
[`ef.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md),
[`legoft`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/legoft.md)

## Examples

``` r
set.seed(1)
n  <- 150
x1 <- runif(n, -3, 3); x2 <- rnorm(n)
# a model with an omitted quadratic term
y  <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.9 * (x1^2 - mean(x1^2))))
fit <- glm(y ~ x1 + x2, family = binomial())
deepgof1(fit, B = 49)   # B = 49 to keep the example fast; use the default in practice
#> 
#>  DeepGOF-1: pretrained residual-map goodness-of-fit test
#> 
#> data:  fit
#> S = 6.8465, B = 49, p-value = 0.0200
#> grid: 6 x 6 over ranks of x1 and x2
#> alternative hypothesis: the logistic model is misspecified
#> 
```
