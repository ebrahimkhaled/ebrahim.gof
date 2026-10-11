# DeepGOF-1 for frozen predictions: external validation of a risk model

Tests whether given probabilities are calibrated for given 0/1 outcomes
within subgroups of the covariates, for a model that was fitted
elsewhere and is not refitted: a published risk score checked on new
patients, or the predictions of any model, logistic or not, on a
validation set. The null hypothesis is that each \\y_i\\ is
Bernoulli(\\p_i\\) with \\p_i\\ as given.

## Usage

``` r
deepgof1.external(
  y,
  p,
  X,
  B = 199L,
  K = 6L,
  reading = c("combined", "axes", "allpairs"),
  axes = NULL
)
```

## Arguments

- y:

  0/1 outcomes.

- p:

  the predicted probabilities to be checked, strictly between 0 and 1,
  one per outcome.

- X:

  a numeric matrix or data frame of covariates, one row per outcome,
  with column names. At least one column.

- B:

  number of Monte Carlo draws; the p-value lies on a grid of
  `1 / (B + 1)`.

- K:

  grid resolution; leave at 6, the resolution the shipped weights were
  trained at.

- reading:

  `"combined"` (the default), `"axes"` or `"allpairs"`, as in
  [`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md).
  The combined reading is the exact minimum of the other two, calibrated
  over the same draws.

- axes:

  optional names of two columns of `X` for the axis-rule map, fixing it
  in advance.

## Value

An object of class `"deepgof1"`, as returned by
[`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md).

## Details

The residual map of
[`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md)
is drawn over the ranks of the covariates and standardized by the given
probabilities, and the shipped network scores it. The reference
distribution is built by drawing \\y^\* \sim\\ Bernoulli(\\p\\) and
scoring again; with nothing estimated, the law of the score is the same
under the whole null, so the rank p-value is exactly valid at every
sample size (Besag and Clifford 1989).

Because the map lays the residuals out over the covariates, the test
checks calibration within patient subgroups (strong calibration in the
hierarchy of Van Calster et al. 2016). A wrong overall rate or a wrong
calibration slope is the same for every patient and is found better by
tests along the predicted risk, such as the calibration belt
([`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
with GiViTI); a miscalibration that differs between patients with the
same predicted risk, a missed U-shape, threshold or interaction, is
found better by this test, which also shows where it lies.

The axis rule needs a ranking of the covariates by their effect. Without
coefficients, it regresses the logit of `p` on the covariates by least
squares and scores each covariate by \\\|b_j\| \hat\sigma_j\\; for a
published logistic model in these covariates this recovers its own
coefficients exactly. Name the axes with `axes` to fix them in advance.

## References

Ebrahim EK, Hussein OAE-A, El-Kotory A (2026). "Where Does a Logistic
Risk Model Fail? An Audited Neural Goodness-of-Fit Test for Model
Development and External Validation." arXiv:2609.29575 \[stat.ME\].
[doi:10.48550/arXiv.2609.29575](https://doi.org/10.48550/arXiv.2609.29575)

Besag, J. and Clifford, P. (1989). Generalized Monte Carlo significance
tests. *Biometrika*, 76(4), 633–642.

Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina, M. J.
and Steyerberg, E. W. (2016). A calibration hierarchy for risk models
was defined: from utopia to empirical data. *Journal of Clinical
Epidemiology*, 74, 167–176.

## See also

[`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md)

## Examples

``` r
set.seed(1)
n <- 400
X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)          # the published model
y <- rbinom(n, 1, plogis(qlogis(p) + 0.8 * X$x1 * X$x2)) # the new patients
deepgof1.external(y, p, X, B = 99)
#> 
#>  DeepGOF-1 for frozen predictions: exact Monte Carlo calibration test
#> 
#> data:  y and p
#> S = 4.2987, B = 99, p-value = 0.0200
#> grid: 6 x 6 over ranks of x1 and x2
#> combined reading: axis rule p = 0.0100, all-pairs p = 0.0200
#> alternative hypothesis: the probabilities are miscalibrated within covariate subgroups
#> 
```
