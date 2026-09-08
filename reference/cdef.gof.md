# Covariate-Space Directed Ebrahim-Farrington (CDEF) Goodness-of-Fit Test

A directed goodness-of-fit test for binary logistic regression whose
direction lives in *covariate space* (functions of the predictors)
rather than in fitted-probability space like
[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md).
It projects the standardized residuals onto a covariate-space basis
(polynomials and pairwise products, natural splines, or a combination
that also includes fitted-probability bends) and calibrates the
quadratic form with the Farrington estimation-adjusted projection,
exactly as in `def.gof`. This makes it sensitive to omitted interactions
and to local / oscillatory departures that fitted-probability grouping
can miss.

## Usage

``` r
cdef.gof(
  object,
  predicted_probs = NULL,
  X = NULL,
  basis = c("poly", "spline", "combined"),
  method = c("satterthwaite", "imhof")
)
```

## Arguments

- object:

  A fitted binary logistic [`glm`](https://rdrr.io/r/stats/glm.html), or
  a binary (0/1) response vector `y` (then supply `predicted_probs` and
  `X`).

- predicted_probs:

  Numeric predicted probabilities; required when `object` is a `y`
  vector.

- X:

  Design/covariate matrix (with or without an intercept column);
  required when `object` is a `y` vector. Ignored when `object` is a
  glm.

- basis:

  One of `"poly"` (squares, cubes, pairwise products), `"spline"`
  (natural cubic splines per covariate plus a pairwise term; needs
  splines), or `"combined"` (covariate polynomials plus
  fitted-probability bends).

- method:

  One of `"satterthwaite"` (default) or `"imhof"`.

## Value

A one-row `data.frame` with `Test`, `Basis`, `Test_Statistic`, `df`,
`Method`, and `p_value`.

## Details

Let \\\tilde r_i=(y_i-\hat p_i)/\sqrt{\hat p_i(1-\hat p_i)}\\ be the
standardized residuals and \\Z\\ a covariate-space basis matrix. The
statistic is \\S=(Z'\tilde r)'(Z'Z)^{-1}(Z'\tilde r)\\, whose null
distribution is a weighted sum of \\\chi^2_1\\ variables with weights
the eigenvalues of \\(Z'Z)^{-1}Z'\Omega Z\\, where
\\\Omega=I-V^{1/2}X(X'VX)^{-1}X'V^{1/2}\\ adjusts for estimating
\\\hat\beta\\. The p-value uses a Satterthwaite scaled-\\\chi^2\\
approximation (default) or Imhof's method (CompQuadForm). Rank-deficient
bases are reduced automatically.

## References

Farrington, C. P. (1996). On Assessing Goodness of Fit of Generalized
Linear Models to Sparse Data. *JRSS-B* 58(2), 349-360.

## See also

[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md),
[`ef.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md),
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md).

## Examples

``` r
set.seed(1)
n <- 600; x1 <- runif(n, -3, 3); x2 <- rnorm(n)
# truth has an omitted interaction; fit the additive model
y <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.4 * x1 * x2))
fit <- glm(y ~ x1 + x2, family = binomial())
cdef.gof(fit)                    # covariate-space directed test (poly basis)
#>                          Test Basis Test_Statistic       df        Method
#> 1 Covariate-space Directed EF  poly       86.40403 4.271854 satterthwaite
#>        p_value
#> 1 2.805972e-20
cdef.gof(fit, basis = "spline")  # for local / oscillatory misfit
#>                          Test  Basis Test_Statistic       df        Method
#> 1 Covariate-space Directed EF spline       89.49195 7.099734 satterthwaite
#>        p_value
#> 1 1.112775e-16
```
