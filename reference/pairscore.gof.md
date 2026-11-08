# Pair score test: squares and product of the DeepGOF-1 pair

The pair score test, an extra check that comes with
[`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md):
it takes the pair of covariates that DeepGOF-1's default reading selects
and tests whether the fitted logistic model needs the squares of the two
covariates and their product. The statistic is the score (Lagrange
multiplier) statistic for adding \\x_a^2\\, \\x_b^2\\ and \\x_a x_b\\ to
the model, computed at the fit, so nothing is refitted to obtain it.

## Usage

``` r
pairscore.gof(fit, B = 199L, ...)
```

## Arguments

- fit:

  a fitted `glm` with `family = binomial()`, a 0/1 outcome and at least
  two covariates.

- B:

  number of parametric-bootstrap replicates; the p-value lies on a grid
  of `1/(B+1)`.

- ...:

  further arguments passed to
  [`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md):
  `K` and `reading` (which do not change the score test) or
  `covariates`.

## Value

An object of class `"pairscore.gof"`: a list with `statistic` (the score
statistic), `df`, `p.value` (the bootstrap p-value), `p.asymptotic` (the
chi-squared p-value), `B`, `axes` (the selected pair), `boot` (the `B`
bootstrap statistics), `method` and `data.name`.

## Details

The pair is chosen from the data by the axis rule of
[`deepgof1()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md):
the two covariates whose terms contribute most to the linear predictor.
The pair is therefore chosen again in every bootstrap sample, and the
p-value is the rank of the observed statistic among the `B` bootstrap
statistics. This keeps the level although the pair was picked by looking
at the data. The chi-squared p-value on `df` degrees of freedom is also
returned (`p.asymptotic`); it ignores the selection of the pair and is
given for comparison only.

The square of a binary covariate is the covariate itself, so it is left
out and the test has two degrees of freedom; a column that the model
already contains is dropped in the same way. Factors enter by their
level codes, as in
[`deepgof1()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md),
so the test is meant for numeric and binary covariates.

The test uses the same pair selection, the same parametric bootstrap and
the same random numbers as `deepgof1(fit, B)`: with the same seed,
`pairscore.gof(fit, B)` and
`deepgof1(fit, B, score.test = TRUE)$score.test` give the same result,
and the latter also returns DeepGOF-1 itself from the same refits.

Use it when the misfit is expected to be a smooth curvature or an
interaction of the two leading covariates; it then has more power than
DeepGOF-1. Use DeepGOF-1 when the shape of the misfit is unknown (for
example a local bump or an oscillation) and to see where on the map the
model fails.

## References

Ebrahim EK, Hussein OAE-A, El-Kotory A (2026). "Where Does a Logistic
Regression Model Fail? A Pretrained Neural Goodness-of-Fit Test with
Bootstrap Validity." Submitted.

Pregibon, D. (1980). Goodness of link tests for generalized linear
models. *Journal of the Royal Statistical Society, Series C* **29**,
15–24. [doi:10.2307/2346405](https://doi.org/10.2307/2346405)

## See also

[`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md),
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)

## Examples

``` r
set.seed(1)
n  <- 150
x1 <- runif(n, -3, 3); x2 <- rnorm(n)
y  <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.9 * (x1^2 - mean(x1^2))))
fit <- glm(y ~ x1 + x2, family = binomial())
pairscore.gof(fit, B = 49)   # B = 49 to keep the example fast; use the default in practice
#> 
#>  Pair score test: squares and product of the DeepGOF-1 pair
#> 
#> data:  fit
#> T = 55.8153, df = 3, B = 49, p-value = 0.0200 (chi-squared p = 4.6e-12)
#> pair: x1 and x2
#> alternative hypothesis: the model needs the squares or the product of the pair
#> 
```
