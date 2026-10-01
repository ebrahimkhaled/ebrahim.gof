# EDGE: Directed Goodness-of-Fit Test for Binary Logistic Regression

`edge.gof()` is the primary interface to the EDGE test (Efficient
Directed Grouped Examination): a grouped, *directed* goodness-of-fit
test for binary logistic regression under sparse data. EDGE projects the
grouped standardized residuals onto a small pre-specified basis of
calibration shapes (cubic `"poly3"` by default) and refers the resulting
quadratic form to its closed-form weighted chi-squared null distribution
– no refit, no resampling, no tuning.

`edge.gof()` computes exactly the same statistic as the legacy name
[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md)
(retained for backward compatibility); the returned `Test` label is
`"EDGE"`.

## Usage

``` r
edge.gof(
  object,
  predicted_probs = NULL,
  X = NULL,
  G = c("auto", 10),
  basis = "poly3",
  method = "satterthwaite",
  weights = "unit",
  external = FALSE,
  y = NULL
)
```

## Arguments

- object:

  A fitted binary logistic [`glm`](https://rdrr.io/r/stats/glm.html), or
  a binary (0/1) response vector `y` (then supply `predicted_probs`).

- predicted_probs:

  Numeric predicted probabilities; required when `object` is a `y`
  vector, ignored when it is a glm.

- X:

  Optional design matrix, used only with the `y`/`predicted_probs` form:
  it enables the exact estimation-adjusted (\\\Omega\\) calibration
  (logit working weights assumed). Without it the conservative
  \\\chi^2_k\\ reference is used and a warning is issued. Ignored when
  `object` is a glm, and ignored (with a warning) when
  `external = TRUE`.

- G:

  Number of groups: `"auto"`, a number, or a vector of them; the default
  `c("auto", 10)` reports both partitions.

- basis:

  One of `"poly3"` (default), `"poly2"`, `"stukel"`, `"sym"`, or
  `"ensemble"`. `"sym"` is one column, \\\eta\|\eta\|\\ at the logit
  \\\eta\\ of each group's mean fitted risk: Stukel's (1988) symmetric
  direction, aimed at tails that are too heavy or too light on both
  sides, for example a probit or cauchit truth fitted by a logit.

- method:

  One of `"satterthwaite"` (default) or `"imhof"`. Ignored when
  `weights = "score"`.

- weights:

  `"unit"` (default) is the statistic as published, \\S = r'P_Z r\\
  referred to a weighted chi-squared law. `"score"` multiplies each
  column by the square root of its group's variance, which for a logit
  fit makes the statistic the score test for adding the grouped shape to
  the model (a score-type test for other links). It is referred to
  chi-squared on the rank of its information matrix, which is the number
  of columns unless one is redundant (see Details).

- external:

  Logical, default `FALSE`. `TRUE` treats the predicted probabilities as
  frozen (external validation of a fixed model): \\\Omega = I\\, a
  constant column joins the basis, and the statistic is referred to
  \\\chi^2\\ on the number of basis columns (see Details of
  [`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md)).
  Supply `y` and `predicted_probs`, or a glm whose fitted probabilities
  are then taken as frozen. Requires `weights = "unit"` and a basis
  other than `"ensemble"`; `method` is not used.

- y:

  Optional alias for `object` when frozen predictions are tested:
  `edge.gof(y = y, predicted_probs = p, external = TRUE)`.

## Value

A `data.frame` with one row per partition and columns `Test` (`"EDGE"`),
`Partition`, `Role` (`"verdict"` for the first row, `"check"` for the
second), `Basis`, `Test_Statistic`, `df`, `Method` and `p_value`, as
documented in
[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md).

## Details

**Two partitions, side by side.** By default the test is reported at two
partitions, as the EDGE paper recommends. The *default partition*,
`G = "auto"` (\\\max(10, \lceil n/25 \rceil)\\ groups), refines with the
sample and has more power, above all against misfit at the extremes of
risk. *Ten groups*, `G = 10`, keep many records in the extreme groups
and so tolerate more corrupted records, at some cost in power. Choose
which one decides in the analysis plan, from what is known about how the
data were collected, not after seeing the results. The first row is
marked `Role = "verdict"` and the second `Role = "check"`: with the
default `G` the default partition decides and ten groups are the
robustness check; to let ten groups decide, give `G = c(10, "auto")`. If
only the default partition rejects, the signal sits in the extreme
groups: check those records. Give a single `G` to get one row.

## References

Ebrahim EK, El-Kotory A (2026). "EDGE: A Closed-Form Directed
Goodness-of-Fit Test for Sparse Logistic Regression." arXiv:2608.20511
\[stat.ME\].
[doi:10.48550/arXiv.2608.20511](https://doi.org/10.48550/arXiv.2608.20511)

Ebrahim EK, El-Kotory A (2026). "A Directional Hosmer-Lemeshow
Goodness-of-Fit Test for Sparse Logistic Regression." arXiv:2607.15454
\[stat.ME\].
[doi:10.48550/arXiv.2607.15454](https://doi.org/10.48550/arXiv.2607.15454)

## See also

[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md)
(legacy name),
[`ef.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md),
[`def.ensemble.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md),
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md).

## Author

Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>

## Examples

``` r
set.seed(1)
x <- runif(500, -3, 3)
y <- rbinom(500, 1, plogis(0.6 * x))
fit <- glm(y ~ x, family = binomial())
edge.gof(fit)                      # cubic basis, at the default partition and at ten groups
#>   Test           Partition    Role Basis Test_Statistic       df        Method
#> 1 EDGE    default (G = 20) verdict poly3      0.9477818 2.006467 satterthwaite
#> 2 EDGE ten groups (G = 10)   check poly3      1.2911935 2.024052 satterthwaite
#>     p_value
#> 1 0.6221428
#> 2 0.5264813
edge.gof(fit, G = 10)              # one partition only
#>   Test           Partition    Role Basis Test_Statistic       df        Method
#> 1 EDGE ten groups (G = 10) verdict poly3       1.291194 2.024052 satterthwaite
#>     p_value
#> 1 0.5264813
edge.gof(fit, basis = "stukel")    # Stukel-shape basis
#>   Test           Partition    Role  Basis Test_Statistic       df        Method
#> 1 EDGE    default (G = 20) verdict stukel      0.8736174 1.843343 satterthwaite
#> 2 EDGE ten groups (G = 10)   check stukel      1.2225715 1.858582 satterthwaite
#>     p_value
#> 1 0.5537401
#> 2 0.4449660
edge.gof(fit, basis = "sym")       # Stukel's symmetric direction, one column
#>   Test           Partition    Role Basis Test_Statistic df        Method
#> 1 EDGE    default (G = 20) verdict   sym     0.08741675  1 satterthwaite
#> 2 EDGE ten groups (G = 10)   check   sym     0.17947447  1 satterthwaite
#>     p_value
#> 1 0.3833456
#> 2 0.2175474
edge.gof(fit, basis = "sym", weights = "score")   # its score form
#>   Test           Partition    Role Basis Test_Statistic df Method   p_value
#> 1 EDGE    default (G = 20) verdict   sym      0.8681059  1  score 0.3514802
#> 2 EDGE ten groups (G = 10)   check   sym      1.6959547  1  score 0.1928179
edge.gof(fit, G = "auto")          # the default partition alone: max(10, ceiling(n / 25)) = 20 here
#>   Test        Partition    Role Basis Test_Statistic       df        Method
#> 1 EDGE default (G = 20) verdict poly3      0.9477818 2.006467 satterthwaite
#>     p_value
#> 1 0.6221428
```
