# Directed Ebrahim-Farrington (DEF) Goodness-of-Fit Test

Performs the Directed Ebrahim-Farrington (DEF) goodness-of-fit test for
a fitted binary logistic regression model. DEF concentrates its power on
a small set of calibration-curve "shape" directions by projecting the
grouped standardized residuals onto a low-dimensional basis and testing
the squared length of that projection.

**Naming note:** this test is published under the name **EDGE**
(Efficient Directed Grouped Examination), and
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
  basis = c("poly3", "poly2", "stukel", "sym", "ensemble"),
  method = c("satterthwaite", "imhof"),
  weights = c("unit", "score"),
  external = FALSE
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

  Integer number of equal-frequency groups (default 10; must be \>= 3),
  or `"auto"` for `max(10, ceiling(n / 25))`, the partition rule of the
  EDGE paper.

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
  \\\chi^2\\ on the number of basis columns (see Details of `def.gof`).
  Supply `y` and `predicted_probs`, or a glm whose fitted probabilities
  are then taken as frozen. Requires `weights = "unit"` and a basis
  other than `"ensemble"`; `method` is not used.

## Value

A one-row `data.frame` with columns `Test`, `Basis`, `Test_Statistic`
(the statistic \\S\\), `df`, `Method`, and `p_value`. For
`weights = "score"`, `Method` is `"score"` and `df` is the integer rank
the statistic is referred to. For `external = TRUE`, `Method` is
`"external"` and `df` is the integer number of basis columns, constant
included. When `basis = "ensemble"`, the return is that of
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
`"poly2"`, `"poly3"` (default), `"stukel"`, `"sym"`; `"ensemble"` runs
`"poly2"`, `"poly3"` and `"stukel"` and combines them via
[`def.ensemble.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md).

Equal-frequency groups split tied fitted risks by row order. With many
ties, as with grouped data or a model on discrete covariates, the result
can therefore depend on the order of the rows, and randomising the row
order is advised.

With `weights = "score"` each column of \\Z\\ is multiplied by
\\\sqrt{V_g}\\, the square root of its group's variance, so that \\Z'r\\
becomes \\\sum_g z_g (O_g - E_g)\\: for a logit fit, the score for
adding the grouped shape to the model as a step covariate. Its
information after adjusting for the fitted coefficients is \\Z'\Omega
Z\\, and the statistic \\u'I^{-1}u\\ is referred to a \\\chi^2\\ law on
the rank of that information (the number of columns unless one is
redundant), read from it after scaling to a correlation matrix. For a
logit fit this is the Rao score test for adding the grouped columns, and
it agrees with `anova(..., test = "Rao")` up to glm's convergence
tolerance; for other links it is a score-type test. A column whose
information after the fit is below \\10^{-10}\\ times its information
before the fit (\\Z_s'Z_s\\, with \\Z_s\\ the weighted columns) is one
the model already spans, as when the fitted logit is constant. It is
left out, and when no column is left the p-value is `NA`, with a warning
of class `def_no_information`. The unit form is the statistic as
published; the score form keeps a shape on the logit scale from losing
its signal when the group variances differ strongly, as they do at high
discrimination.

**External mode.** With `external = TRUE` the predicted probabilities
are taken as frozen, as when a published model is checked on new data
with its coefficients fixed. Nothing is estimated from these data, so
\\\Omega = I\\ exactly, and no score equation absorbs the overall level,
so a column of ones is added to the basis before redundant columns are
dropped. The statistic \\S = r'P_Z r\\ is then referred to \\\chi^2\\ on
the number of columns kept: \\d + 1\\ for a \\d\\-column basis (4 for
`"poly3"` and `"stukel"`, 3 for `"poly2"`, 2 for `"sym"`; one fewer for
`"stukel"` when every group lies on one side of 0.5). The groups are the
same equal-frequency groups as in the default mode. No "conservative"
warning is given, since \\\Omega = I\\ is exact here, and `Method` is
`"external"`. Only `weights = "unit"` is available: the score form is a
different projection, and it has not been validated for frozen
predictions. When `object` is a glm, its response and fitted
probabilities are used as the frozen predictions; the fit is not
otherwise used, so this is a test of those predictions and not the
estimation-adjusted test of the model.

With fewer events (or fewer non-events) than groups, the grouped
reference distribution is unreliable. The p-value is still returned,
with a warning; a smaller `G` avoids it. With no event, or no non-event,
the model has no maximum-likelihood fit, and the p-value is `NA`, with a
warning of class `def_degenerate`.

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
def.gof(wrong, basis = "sym")        # symmetric tail direction, one column
#>                          Test Basis Test_Statistic df        Method     p_value
#> 1 Directed Ebrahim-Farrington   sym       1.518808  1 satterthwaite 0.009579613
def.gof(wrong, weights = "score")    # score form of the poly3 basis
#>                          Test Basis Test_Statistic df Method    p_value
#> 1 Directed Ebrahim-Farrington poly3       10.47366  3  score 0.01494063
def.gof(wrong, basis = "ensemble")   # combine poly2, poly3 and stukel (CCT)
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

## external validation: freeze the model fitted on one half, test it on the other
dev <- gof_demo[1:250, ]; val <- gof_demo[-(1:250), ]
frozen <- glm(outcome ~ age + bmi + sex + treatment, data = dev, family = binomial())
p_val <- predict(frozen, newdata = val, type = "response")
def.gof(val$outcome, predicted_probs = p_val, external = TRUE)
#>                          Test Basis Test_Statistic df   Method   p_value
#> 1 Directed Ebrahim-Farrington poly3       5.008141  4 external 0.2864632
```
