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
deepgof1(
  fit,
  B = 199L,
  K = 6L,
  reading = c("axes", "allpairs", "combined", "columns"),
  covariates = NULL
)
```

## Arguments

- fit:

  a fitted `glm` with `family = binomial()`, a 0/1 outcome and at least
  one covariate.

- B:

  number of parametric-bootstrap replicates. The p-value lies on a grid
  of `1/(B+1)`, so `B = 199` makes the nominal .05 attainable exactly.

- K:

  grid resolution. Leave at 6: the shipped weights were trained at
  `K = 6` and are not valid at any other resolution.

- reading:

  `"axes"` (the default), `"allpairs"`, `"combined"` or `"columns"` (the
  rule of versions 2.7.0 and 2.8.0). See Details for when to use which.

- covariates:

  optional character vector naming the covariates to form the pairs
  from, for the all-pairs and combined readings. By default every
  covariate of the model that is not constant is used.

## Value

An object of class `"deepgof1"`: a list with `statistic` (the observed
score), `p.value`, `B`, `K`, `reading`, `axes` (the covariates of the
map that gave the statistic), `map` (that map, a `K` by `K` matrix of
standardized residual sums whose rows follow the first axis), `pairs`
(for the all-pairs and combined readings, the observed score of every
pair), `components` (for `"combined"`, the p-values of the axis rule and
the all-pairs reading), `boot` (the `B` bootstrap statistics of the
reported statistic) and `method`.

## Details

A residual map is a `K` by `K` grid over the empirical ranks of two
covariates; each cell holds a standardized residual sum, approximately
standard normal under a correct model. Misfit therefore has a location
on the map – an omitted quadratic paints a stripe, an omitted
interaction a saddle – which is what the convolutional statistic reads.

Four readings of the map are offered. The default, `reading = "axes"`,
draws one map, over the two covariates whose terms contribute most to
the linear predictor (the standard deviation of each covariate's total
contribution; for a covariate that enters as one untransformed column
this is \\\|\hat\beta_j\| \hat\sigma_j\\), and chooses them again in
every bootstrap replicate. A covariate that enters as `ns(x, 3)`,
`poly(x, 2)`, `I(x^2)` or a factor is one covariate, and the map is
drawn over the ranks of `x` itself. When the misfit lies on the
covariates with the strongest effects, among others that carry little
signal, this rule finds them almost every time and has the most power.
It reads a covariate by its effect in the fitted model, so it can pass
over a covariate whose effect is a pure U-shape with no linear slope.

`reading = "allpairs"` scores the map of every pair of covariates and
takes the largest score as the statistic. The same maximum is taken in
every bootstrap replicate, so the p-value needs no correction for the
choice of pair. It does not depend on the fitted effects, so it finds a
U-shaped covariate, and with few covariates (three or so) it has more
power than the default; with many covariates the maximum over
\\p(p-1)/2\\ pairs costs power when one pair carries the misfit. The
pair that reaches the maximum is returned as `axes`, and its map as
`map`, so the result also says where the misfit lies; `covariates`
restricts the pairs to a chosen set.

`reading = "combined"` runs both from the same bootstrap refits and
reports the smaller of their two p-values, calibrated exactly: the
observed data and the `B` replicates are exchangeable under the null, so
the observed minimum is ranked among the `B + 1` minima. It costs no
more than the all-pairs reading. `reading = "columns"` is the rule of
versions 2.7.0 and 2.8.0, over two model-matrix columns, kept to
reproduce earlier results; for models whose covariates all enter as one
untransformed column it gives the same p-value as the default.

Transformed covariates are read from the data the model was fitted to;
factors enter by their level codes.

With one covariate there is no pair to choose, and every reading is the
same map of 36 quantile cells along its ranks. The shipped network was
trained on two-covariate maps only: on such models the bootstrap still
gives it its level, but it has less power than a network trained on this
map would.

Ties among covariate values, as with binary, categorical or rounded
covariates, are broken at random. The random order is drawn once per
call and used for the observed map and for every bootstrap map, so the
p-value does not depend on the order of the rows in the data. With
heavily tied axes the p-value can vary noticeably from one seed to the
next; report the seed.

The test is a small-sample instrument. Against the classical partition
tests it gains most at \\n\\ of 50 to 200 and the gain decays as \\n\\
grows; because each map uses two covariates at a time, misfit that
depends on three or more covariates jointly is harder for it to see than
for covariate-space or smoothing tests. It is not one of the tests
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
selects: call it directly on the same fitted model and read its p-value
beside the panel. See
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
for the classical battery.

## Reproducibility

A bootstrap refit that fails to converge is scored `+Inf`, so it counts
against rejection – the conservative direction. Set a seed before
calling for a reproducible p-value. When some covariate column has tied
values, the random tie-breaking uses the same seed, so from version
2.8.0 such calls give a different p-value for a given seed than earlier
versions did; calls without ties give the same p-value as before.
Version 2.9.0 refits each bootstrap sample on the fitted model's design
matrix. Earlier versions refitted the formula on the model frame, which
fails for every term that transforms a covariate (`log(x)`, `ns(x, 3)`,
`poly(x, 2)`): each replicate was then scored `+Inf` and the p-value
was 1. For models without such terms the default reading gives the same
p-value as in 2.8.0.

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
# every pair of covariates, with the pair where the misfit is largest
deepgof1(fit, B = 49, reading = "allpairs")$axes
#> [1] "x1" "x2"
```
