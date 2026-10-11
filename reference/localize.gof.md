# Localize the misfit of a fitted logistic model, with error control

The in-sample counterpart of
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md):
checks a logistic `glm` on the data it was fitted to and names which
part of the misfit is present, with the familywise error rate held at
`alpha`. After a fit with an intercept the `INTERCEPT` and `SLOPE` parts
are identically zero, so two groups remain:

- `LINK`:

  bends of the map from the linear predictor to risk.

- `COV`:

  misfit among observations that share a fitted risk: a missed
  interaction or curved covariate.

The bases are those of
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md),
also made orthogonal to the model matrix. The reference is the
parametric bootstrap: `B` outcome vectors are drawn from the fitted
model, the model is refitted to each on the same design matrix, and
every member p-value is recomputed. A refit that fails is dropped; the
number used is returned.

## Usage

``` r
localize.gof(
  fit,
  X = NULL,
  B = 199L,
  alpha = 0.05,
  dealias = FALSE,
  cov_df = 5,
  naming = c("closure", "holm", "bonferroni"),
  robust = FALSE,
  seed = NULL,
  plot = interactive()
)
```

## Arguments

- fit:

  a fitted `glm` with `family = binomial(link = "logit")`, a 0/1
  outcome, one row per observation, no prior weights and no offset.

- X:

  a numeric matrix or data frame of covariates, one row per observation.
  By default the numeric columns of the model frame other than the
  response; a term such as `log(x)` enters as its column `log(x)`, while
  matrix terms (`ns()`, [`poly()`](https://rdrr.io/r/stats/poly.html))
  and factors are left out. `COV` needs two columns.

- B:

  number of parametric-bootstrap replicates; each p-value lies on a grid
  of `1 / (B_used + 1)`.

- alpha:

  familywise level.

- dealias:

  if `TRUE`, de-alias the `LINK` group from model columns that are
  curved in the score; see Details.

- cov_df, naming:

  as in
  [`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md).

- robust:

  if `TRUE`, build the bases on normal scores; see Details.

- seed:

  optional integer, passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the
  bootstrap. The random number stream of the session is restored on
  exit.

- plot:

  if `TRUE`, draw the verdict with
  [`plot.gof_localize`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/plot.gof_localize.md)
  (a misfit compass beside the lattice of closed tests). The default
  draws it in an interactive session and not in scripts, simulations or
  examples.

## Value

An object of class `"gof_localize"`, as described in
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md);
`draws` is the number of bootstrap replicates used and `dealiased` names
the model columns removed from `LINK` (if any).

## Details

In-sample, a model column whose mean given the score is curved can make
`LINK` read covariate misfit. `dealias = TRUE` selects such columns
once, on the observed fit (a weighted F test of a 3-df spline in the
score against a line, level .05), and makes the `LINK` bases orthogonal
to their spline fits on the observed score, frozen for every bootstrap
draw, so that `LINK` reads link shape only. The price is power of `LINK`
against departures that resemble those columns.

`naming` is as in
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md).
`robust = TRUE` builds the bases on normal scores, so that a few extreme
covariate values cannot drive a verdict: the covariates and the model
columns are transformed once, from the data, and the score of every fit,
observed or refitted, is transformed too (de-aliasing then uses the
transformed score and columns). The level is unaffected, since the
bootstrap still refits the model on the original model matrix.

Only a logistic fit is served, since the bases and the refits assume the
logit link. A model with another link, or any model that is not a `glm`,
can be checked on validation data with
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md),
which needs only its predictions.

## Reading a verdict

Act on the highest group named. `COV`: revise the model; the remaining
parts are assessed again after the revision. `LINK` without `COV`:
adding a function of the score to the model, or changing the link, is
consistent with the data; it is not proof that no revision is needed.
Nothing named: no misfit was detected, a verdict limited by the sample
size.

## References

Ebrahim, E. K., El-Kotory, A. and Hussein, O. A. E.-A. (2026). One
goodness-of-fit test is not enough: error-controlled localization of
misfit in logistic risk models. Preprint. Reproduction materials:
[doi:10.5281/zenodo.23192535](https://doi.org/10.5281/zenodo.23192535)

Marcus, R., Peritz, E. and Gabriel, K. R. (1976). On closed testing
procedures with special reference to ordered analysis of variance.
*Biometrika*, 63(3), 655–660.
[doi:10.1093/biomet/63.3.655](https://doi.org/10.1093/biomet/63.3.655)

Liu, Y. and Xie, J. (2020). Cauchy combination test: a powerful test
with analytic p-value calculation under arbitrary dependency structures.
*Journal of the American Statistical Association*, 115(529), 393–402.
[doi:10.1080/01621459.2018.1554485](https://doi.org/10.1080/01621459.2018.1554485)

Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina, M. J.
and Steyerberg, E. W. (2016). A calibration hierarchy for risk models
was defined: from utopia to empirical data. *Journal of Clinical
Epidemiology*, 74, 167–176.
[doi:10.1016/j.jclinepi.2015.12.005](https://doi.org/10.1016/j.jclinepi.2015.12.005)

Steyerberg, E. W., Borsboom, G. J. J. M., van Houwelingen, H. C.,
Eijkemans, M. J. C. and Habbema, J. D. F. (2004). Validation and
updating of predictive logistic regression models: a study on sample
size and shrinkage. *Statistics in Medicine*, 23(16), 2567–2586.
[doi:10.1002/sim.1844](https://doi.org/10.1002/sim.1844)

## See also

[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md)
for frozen predictions on new data;
[`plot.gof_localize`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/plot.gof_localize.md)
for the display of a verdict.

## Examples

``` r
set.seed(2)
n  <- 500
x1 <- rnorm(n); x2 <- rnorm(n)
y  <- rbinom(n, 1, plogis(-0.3 + 0.8 * x1 + 0.6 * x2 + 0.8 * x1 * x2))
fit <- glm(y ~ x1 + x2, family = binomial())   # the interaction is missed
localize.gof(fit, B = 49, seed = 1)   # B = 49 to keep the example fast; use 199 or more
#> 
#> Localization of misfit: closed testing over orthogonal groups
#> 
#> data:     fit (n = 500)
#> setting:  in-sample checking of a fitted logistic model
#> reference: parametric bootstrap, 49 draws
#> 
#> p-values by group (alone, and adjusted by closed testing):
#>       alone adjusted named
#> LINK 0.0200   0.0200     *
#> COV  0.0200   0.0200     *
#> 
#> groups named (closed testing) at FWER 0.05: LINK + COV
#> recommended update: revise the model; recalibration cannot repair misfit among
#>   patients who share a score
#> 
```
