# Closed-form goodness-of-fit test for penalized logistic regression (CALM)

Refers the shrinkage-corrected grouped goodness-of-fit statistics to an
analytic reference distribution, so that no bootstrap is needed. It is
the companion of
[`shrink.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/shrink.gof.md),
which refers the same statistics to a prepivoting bootstrap: the
statistics are identical, only the reference differs. CALM stands for
Calibration Assessment under Lambda-shrunk Models.

## Usage

``` r
calm.gof(
  X,
  y,
  lambda,
  G = 10,
  basis = c("decile", "adaptive", "edge"),
  lambda_scale = c("theory", "glmnet"),
  tau = 0.2,
  inflate = c("kappa", "none")
)
```

## Arguments

- X:

  numeric matrix or data frame of predictors, without an intercept
  column. Standardize the columns as you would before any ridge fit.

- y:

  numeric or integer vector of 0/1 responses.

- lambda:

  the ridge penalty. On the theory scale by default, that is on the
  scale of the log-likelihood; pass `lambda_scale = "glmnet"` to give
  `glmnet`'s value instead, which is `lambda / n`.

- G:

  number of equal-frequency groups. Default 10.

- basis:

  which statistics to compute: any of `"decile"`, `"adaptive"` and
  `"edge"`. Default: all three.

- lambda_scale:

  `"theory"` (default) or `"glmnet"`.

- tau:

  threshold of the degree rule, in (0, 1). Default 0.2. Values below 0.2
  keep the cubic direction too often in the proportional regime.

- inflate:

  whether to inflate the reference for the selection effect, `"kappa"`
  (default) or `"none"`. On the adaptive basis the inflation changes the
  size by at most 0.008 in the designs examined.

## Value

An object of class `"calm.gof"`: a list with components `SC.HL`,
`SC.EDGE` and `SC.EDGE.adaptive` (each a list with `statistic` and
`p.value`, and for the adaptive basis also the chosen `degree` and the
observable `rho_hat`), together with the penalty on both scales and the
observables used to build the reference.

## Details

Under a ridge penalty the grouped standardized residuals are displaced
by shrinkage. Subtracting an estimate of that displacement restores the
maximum likelihood covariance to first order, but the resulting
reference is conservative, because it standardizes by the Bernoulli
variance at the *fitted* probabilities, which shrinkage inflates towards
one quarter. CALM replaces it by a de-noised estimate of the null
Bernoulli variance, obtained from the fit alone through the observable
adjustments of Bellec (2025), and inflates the result for the effect of
grouping on an index that depends on the response.

Two bases are returned. The EDGE basis projects the corrected residual
onto orthogonal polynomials in the group-mean fitted probability, and
`SC.EDGE.adaptive` chooses the degree from the data: it keeps the cubic
direction when the aspect ratio is small, or when the observable index
correlation \\\hat\rho\\ satisfies \\\hat\rho^{6} \ge \tau\\, and uses
degree two otherwise. That basis is the one to prefer. The decile basis
`SC.HL` is the shrinkage-corrected Hosmer-Lemeshow statistic and is
reported for continuity with that tradition; because it spends every
group direction it cannot avoid the direction the selection effect
occupies, and its reference relies on a constant calibrated by
simulation, which does not transfer to every design. See the reference
for the designs in which it fails.

## Scope

The reference is validated for aspect ratios \\p/n\\ up to 0.4 and for
penalties that shrink towards zero; the lasso is not covered, because
the displacement requires a differentiable penalty. The test asks
whether the logistic form is correct along the fitted index. It does not
ask whether the shrunk probabilities are calibrated, which under a
penalty they are not, by an amount the analyst chose when selecting
`lambda`.

The reference is validated for outcomes that are not strongly
unbalanced. Once \\p/n\\ is an appreciable fraction the level is lost as
the prevalence falls: at \\p/n = 0.25\\ the smooth basis rejects 0.140,
0.574 and 0.884 of correctly specified models at prevalences 0.30, 0.15
and 0.08, and the decile basis 0.060, 0.204 and 0.492. At fixed
dimension the decile basis is unaffected and the smooth one degrades far
more slowly, to 0.130 at prevalence 0.08. What fails there is the
response-dependent grouping inherited from the Hosmer-Lemeshow
construction rather than the reference itself. With few events and
\\p/n\\ an appreciable fraction, neither statistic is validated and
[`shrink.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/shrink.gof.md)
is the less badly behaved of the two.

## References

Bellec, P. C. (2025). Observable adjustments in single-index models for
regularized M-estimators with bounded p/n. *The Annals of Statistics*,
**53**(2), 531–560.

Davies, R. B. (1980). Algorithm AS 155: the distribution of a linear
combination of chi-squared random variables. *Journal of the Royal
Statistical Society, Series C*, **29**(3), 323–333.

Ebrahim, E. K. (2026). A closed-form reference distribution for
goodness-of-fit testing under penalized logistic regression in the
proportional regime.

## See also

[`shrink.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/shrink.gof.md)
for the bootstrap reference for the same statistics, and
[`edge.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.gof.md)
for the EDGE test on an unpenalized fit.

## Examples

``` r
set.seed(1)
n <- 300; p <- 20
X <- matrix(rnorm(n * p), n)
y <- rbinom(n, 1, 1 / (1 + exp(-(X[, 1] - 0.5 * X[, 2]))))
calm.gof(X, y, lambda = 100)
#> 
#> CALM: closed-form goodness of fit for penalized logistic regression
#> lambda = 100 (theory scale) = 0.3333 on glmnet's scale | G = 10 | p/n = 0.067
#>   SC.HL             statistic    7.133   p = 0.3733
#>   SC.EDGE           statistic    0.531   p = 0.7304
#>   SC.EDGE.adaptive  statistic    0.531   p = 0.7304   [degree 3, rho_hat 0.79]
#> 
#> SC.EDGE.adaptive is the statistic to prefer; see ?calm.gof for the scope of the reference.
#> 
```
