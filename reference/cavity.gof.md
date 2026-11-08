# Leave-one-out (cavity) Hosmer-Lemeshow test for penalised logistic regression

Tests the fit of a ridge, lasso or elastic-net logistic regression in
high dimensions, including \\p \> n\\. The textbook Hosmer-Lemeshow test
groups on fitted probabilities that are shrunk and have seen their own
responses, and after a penalised fit it rejects most correct models.
This test groups, recalibrates and standardises on the leave-one-out
linear predictor instead, computed from the one fit by a single Newton
step, and refers the statistic to the familiar \\\chi^2\_{G-2}\\ law.

## Usage

``` r
cavity.gof(
  X,
  y,
  lambda,
  alpha = 0,
  G = 10,
  lambda_scale = c("theory", "glmnet"),
  recalibration = c("auto", "logistic", "logistic-normal"),
  smooth.degree = 0,
  split = FALSE,
  reference = c("closed", "simulate"),
  B = 199,
  check.twins = TRUE,
  twin.cut = 0.999
)
```

## Arguments

- X:

  numeric matrix of predictors (no intercept column). Standardise the
  predictors before the call if wanted, not inside a procedure that uses
  all \\n\\ observations to fit each one; the fit itself does not
  rescale them.

- y:

  binary response (0/1).

- lambda:

  penalty. On the theory scale (log-likelihood minus \\\lambda/2\\ times
  the squared norm of the slopes; \\\lambda = n \lambda\_{glmnet}\\)
  unless `lambda_scale = "glmnet"`.

- alpha:

  elastic-net mixing parameter: 0 (default) is the ridge, fitted
  internally; values in (0, 1\] use glmnet (lasso at 1), and the
  leave-one-out index is computed on the active set.

- G:

  number of groups (default 10).

- lambda_scale:

  `"theory"` (default) or `"glmnet"`.

- recalibration:

  `"auto"` (logistic for \\n \le 5000\\, logistic-normal above),
  `"logistic"` or `"logistic-normal"`.

- smooth.degree:

  degree of the smooth version (orthogonal polynomials in the group
  means of the index); 0 skips it.

- split:

  if `TRUE`, also compute the sample-split test (index from a fit on the
  other half, so no approximation and no dependence between indices) and
  its cross-fitted version, which combines the two halves by the Cauchy
  combination.

- reference:

  `"closed"` (default) or `"simulate"`: refer \\S\\ also to its
  distribution under the fitted model, by simulating responses from the
  fitted probabilities, refitting and recomputing \\S\\, `B` times. Use
  it when `twin.share` is not near zero.

- B:

  number of simulated data sets for `reference = "simulate"`.

- check.twins:

  if `TRUE` (default), compute the share of observations with a twin. It
  costs of order \\n^2\\ times the size of the active set.

- twin.cut:

  cosine above which two observations count as twins (default 0.999).

## Value

An object of class `"cavity.gof"`: a list with `statistic` (\\S\\),
`df`, `p.value` (\\\chi^2\_{G-2}\\), `p.exact` (the weighted chi-squared
law of the grouped residuals), `p.smooth`, `p.split`, `p.crossfit`,
`p.simulated` (when requested), `twin.share`, `table` (observed and
expected counts and standardised residuals by group), `recalibration`
(fitted intercept, slope and noise), `index` (the cavity index), `beta`,
`lambda`, `alpha`.

## Details

The reference law is derived for predictors that behave like independent
draws. When many rows of `X` repeat, observations with the same row
share the error of their index and the closed-form reference rejects too
often; in a data set where half of the rows repeated, it rejected 57 per
cent of correct models. The function then warns, and
`reference = "simulate"` gives a reference conditional on the design.

## References

Ebrahim, E. K. and Arashi, M. (2026). A leave-one-out goodness-of-fit
test for penalised logistic regression in high dimensions. Submitted.
Code archive:
[doi:10.5281/zenodo.23861454](https://doi.org/10.5281/zenodo.23861454) .

## Examples

``` r
set.seed(1)
n <- 300; p <- 60
X <- matrix(rnorm(n * p), n)
y <- rbinom(n, 1, 1 / (1 + exp(-drop(X %*% rnorm(p, sd = 0.3)))))
res <- cavity.gof(X, y, lambda = 50)
res
#> 
#> Leave-one-out (cavity) Hosmer-Lemeshow test
#>   fit: ridge, lambda = 50 (theory scale), n = 300, active predictors = 60
#>   recalibration: logistic, slope = 1.531
#>   S = 8.794 on 8 df:  P = 0.36  (exact weighted law: P = 0.3626)
#>   share of observations with a twin: 0.000
```
