# Directed Ebrahim-Farrington (DEF) Goodness-of-Fit Test

Performs the Directed Ebrahim-Farrington (DEF) goodness-of-fit test for
a fitted binary logistic regression model. DEF concentrates its power on
a small set of calibration-curve "shape" directions by projecting the
grouped standardized residuals onto a low-dimensional basis and testing
the squared length of that projection.

**Naming note:** this test is published under the name **EDGE** (Ebrahim
Directed Goodness-of-fit Evaluation), and
[`edge.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.gof.md)
is the primary interface going forward. `def.gof()` is retained,
unchanged, as a fully supported legacy name.

## Usage

``` r
def.gof(
  object,
  predicted_probs = NULL,
  X = NULL,
  G = 10,
  basis = c("poly3", "poly2", "stukel", "ensemble"),
  method = c("satterthwaite", "imhof")
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

A one-row `data.frame` with columns `Test`, `Basis`, `Test_Statistic`
(the statistic \\S\\), `df`, `Method`, and `p_value`. When
`basis = "ensemble"`, the return is that of
[`def.ensemble.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md).

## Details

The observations are sorted by predicted probability and split into `G`
equal-frequency groups; the standardized grouped residual vector \\r\\
is projected onto a basis matrix \\Z\\ of smooth shapes, giving \\S =
(Z'r)'(Z'Z)^{-1}(Z'r)\\. Its null distribution is a weighted sum of
\\\chi^2_1\\ variables with weights equal to the eigenvalues of
\\(Z'Z)^{-1}Z'\Omega Z\\, where \\\Omega = I - U(X'WX)^{-1}U'\\ is the
estimation-adjusted covariance of the grouped residuals. The p-value
uses a Satterthwaite scaled-\\\chi^2\\ approximation (default) or
Imhof's method (if the CompQuadForm package is installed). Bases:
`"poly2"`, `"poly3"` (default), `"stukel"`; `"ensemble"` runs all three
and combines them via
[`def.ensemble.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md).

## References

Ebrahim EK, El-Kotory A (2026). "A Directional Hosmer-Lemeshow
Goodness-of-Fit Test for Sparse Logistic Regression." arXiv:2607.15454
\[stat.ME\].
[doi:10.48550/arXiv.2607.15454](https://doi.org/10.48550/arXiv.2607.15454)

Ebrahim EK, El-Kotory A (2026). "EDGE: A Closed-Form Directed
Goodness-of-Fit Test for Sparse Logistic Regression." arXiv:2608.20511
\[stat.ME\].
[doi:10.48550/arXiv.2608.20511](https://doi.org/10.48550/arXiv.2608.20511)

## See also

[`ef.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md),
[`def.ensemble.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md).

## Author

Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>

## Examples

``` r
## gof_demo carries a documented smooth calibration misfit: the risk bends in age,
## and a model linear in age misses it. The point of a directed test is to see that.
data("gof_demo", package = "ebrahim.gof")
wrong <- glm(outcome ~ age + bmi + sex + treatment,
             data = gof_demo, family = binomial())
def.gof(wrong)                       # default poly3 basis
#>                          Test Basis Test_Statistic       df        Method
#> 1 Directed Ebrahim-Farrington poly3       8.329512 2.072814 satterthwaite
#>      p_value
#> 1 0.01496412
def.gof(wrong, basis = "stukel")     # tail-shape basis
#>                          Test  Basis Test_Statistic       df        Method
#> 1 Directed Ebrahim-Farrington stukel       6.499717 1.564118 satterthwaite
#>     p_value
#> 1 0.0118151
def.gof(wrong, basis = "ensemble")   # combine all three (CCT)
#>           Test Combiner         Components k     p_value
#> 1 DEF ensemble      cct poly2+poly3+stukel 3 0.006959379

## give the model the term it was missing, and the same test stands down
right <- glm(outcome ~ poly(age, 2) + bmi + sex + treatment,
             data = gof_demo, family = binomial())
def.gof(right)
#>                          Test Basis Test_Statistic       df        Method
#> 1 Directed Ebrahim-Farrington poly3      0.4259731 2.046178 satterthwaite
#>     p_value
#> 1 0.7966492
```
