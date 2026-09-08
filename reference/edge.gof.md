# EDGE: Directed Goodness-of-Fit Test for Binary Logistic Regression

`edge.gof()` is the primary interface to the EDGE test (Ebrahim Directed
Goodness-of-fit Evaluation): a grouped, *directed* goodness-of-fit test
for binary logistic regression under sparse data. EDGE projects the
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
  G = 10,
  basis = "poly3",
  method = "satterthwaite"
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
  `object` is a glm.

- G:

  Integer number of equal-frequency groups (default 10; must be \>= 3).

- basis:

  One of `"poly3"` (default), `"poly2"`, `"stukel"`, or `"ensemble"`.

- method:

  One of `"satterthwaite"` (default) or `"imhof"`.

## Value

A one-row `data.frame` with columns `Test` (`"EDGE"`), `Basis`,
`Test_Statistic`, `df`, `Method`, and `p_value`, as documented in
[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md).

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
edge.gof(fit)                      # default cubic basis, G = 10
#>   Test Basis Test_Statistic       df        Method   p_value
#> 1 EDGE poly3       1.291194 2.024052 satterthwaite 0.5264813
edge.gof(fit, basis = "stukel")    # Stukel-shape basis
#>   Test  Basis Test_Statistic       df        Method  p_value
#> 1 EDGE stukel       1.222571 1.858582 satterthwaite 0.444966
```
