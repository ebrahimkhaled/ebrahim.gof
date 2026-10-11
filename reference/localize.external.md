# Localize the misfit of frozen predictions, with error control

Says *which part* of the misfit is present when given probabilities are
checked against given 0/1 outcomes: a published risk model on new
patients, or any model's predictions on a validation set. Nothing is
refitted. The misfit is split into four parts that follow the
calibration hierarchy (Van Calster et al. 2016), and each part is named
or not with the familywise error rate held at `alpha`:

- `INTERCEPT`:

  calibration in the large: overall risk too high or too low.

- `SLOPE`:

  the calibration slope: risks too extreme or too modest.

- `LINK`:

  bends of the map from the score \\\eta = \mathrm{logit}(p)\\ to risk.

- `COV`:

  misfit among patients who share a score: a missed interaction or
  curved covariate, which no recalibration can repair.

## Usage

``` r
localize.external(
  y,
  p,
  X,
  M = 999L,
  alpha = 0.05,
  calibration = c("montecarlo", "multiplier"),
  cov_df = 5,
  naming = c("closure", "holm", "bonferroni"),
  robust = FALSE,
  seed = NULL,
  plot = interactive()
)
```

## Arguments

- y:

  0/1 outcomes.

- p:

  the predicted probabilities to be checked, one per outcome, in \\\[0,
  1\]\\. Values of exactly 0 or 1 would give an infinite score and are
  moved to `1e-10` and `1 - 1e-10`. The predictions should be frozen: a
  model's own fitted values are calibrated in the large and in slope by
  construction, so `INTERCEPT` and `SLOPE` read \\p = 1\\, and a member
  p-value near 1 pulls every Cauchy intersection that contains it
  towards 1, which can hide a real `COV` misfit. The function warns in
  that case; use
  [`localize.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.gof.md)
  for a model checked on its own data.

- X:

  a numeric matrix or data frame of covariates, one row per outcome,
  without missing values or constant columns. At least one column; `COV`
  needs two.

- M:

  number of reference draws; each p-value lies on a grid of
  `1 / (M + 1)`.

- alpha:

  familywise level.

- calibration:

  `"montecarlo"` (the default, exact) or `"multiplier"`; see Details.

- cov_df:

  the degrees of freedom of the spline in \\\eta\\ that the `COV` bases
  are made orthogonal to; `"auto"` uses `max(5, ceiling(n^(1/3)))`.

- naming:

  how groups are named: `"closure"` (the default; closed testing over
  Cauchy intersections), `"holm"` or `"bonferroni"` (over the
  single-group p-values); see Details.

- robust:

  if `TRUE`, build the bases on normal scores of the score and the
  covariates; see Details.

- seed:

  optional integer, passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the
  reference draws. The random number stream of the session is restored
  on exit, so a call inside a simulation loop does not reset the loop's
  own draws.

- plot:

  if `TRUE`, draw the verdict with
  [`plot.gof_localize`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/plot.gof_localize.md)
  (a misfit compass beside the lattice of closed tests). The default
  draws it in an interactive session and not in scripts, simulations or
  examples.

## Value

An object of class `"gof_localize"`: a list with `named` (the groups
named, in hierarchy order), `highest` (the highest of them, or `NA`),
`action` (the update the decision table gives for it), `single` (the
p-value of each group tested alone), `adjusted` (each group's adjusted
p-value: under closure the largest p-value of an intersection containing
it, otherwise the Holm or Bonferroni adjustment of the single-group
p-values; a group is named when it is at most `alpha`), `intersection`
(the p-value of every intersection, named as in `"SLOPE+COV"`),
`members` (the chi-square p-value of each member score test), `groups`
(the members of each group), `alpha`, `naming`, `robust`, `setting`,
`calibration`, `draws` (the number of reference draws, `M`), `n`,
`method` and `data.name`.

## Details

Each group is a Cauchy combination (Liu and Xie 2020) of score tests
whose bases are confined to that group's part, weighted-orthogonal to
the parts before it: `1`; \\\eta\\; \\\eta^2, \eta^3\\, Stukel's two
terms and `ns(eta, 4)`; and, for `COV`, squares and cubes of the
covariates, their pairwise products and `ns(x_j, 3)`, each made
orthogonal to every function of \\\eta\\ through `ns(eta, cov_df)`. A
covariate with four or fewer distinct values enters the products only.
With a single covariate every function of it is a function of \\\eta\\,
so `COV` is dropped.

A group is named when every intersection of groups containing it
rejects: closed testing (Marcus et al. 1976), so the probability of
naming any group whose part is absent is at most `alpha`. Each
intersection is tested by the mean Cauchy coordinate of its members'
p-values, referred to its rank among `M` reference draws.

With `calibration = "montecarlo"` (the default) the draws are \\y^\*
\sim\\ Bernoulli(\\p\\). Nothing is estimated, so under no misfit the
p-values are exactly valid at every sample size, and for a part that is
absent the error control holds for local departures in the other parts.
`calibration = "multiplier"` uses a robust variance and sign-flipped
residuals, so that an absent part keeps its level for departures of any
size in the other parts; that guarantee is asymptotic, and in
simulations at \\n \le 4000\\ this reference was liberal. Use it for
study, not yet for decisions.

`naming = "holm"` or `"bonferroni"` names the groups by Holm's or
Bonferroni's correction of their single-group p-values instead of by
closure. Both are closed procedures with Bonferroni intersections, so
they keep the same familywise guarantee; every intersection is still
computed and returned.

`robust = TRUE` builds the bases on normal scores,
`qnorm(rank / (n + 1))`, of the score and of every covariate with more
than four distinct values, so that a few extreme covariate values cannot
drive a verdict. The level is unaffected: the Monte Carlo reference is
exact for any fixed choice of bases.

## Reading a verdict

Act on the highest group named, in the order
`INTERCEPT < SLOPE < LINK < COV`: it points to the lightest update that
repairs the model (Steyerberg et al. 2004).

|             |                                                      |
|-------------|------------------------------------------------------|
| `INTERCEPT` | update the intercept                                 |
| `SLOPE`     | logistic recalibration \\a + b\eta\\                 |
| `LINK`      | flexible recalibration \\f(\eta)\\ or another link   |
| `COV`       | revise the model, even if other groups are named too |

When nothing is named, no misfit was detected; that verdict is limited
by the sample size. When the miscalibration is gross, as after transport
to another population, a large error on the logit scale also bends the
probability curve, so the named rung can be too high. Then split the
validation data at random, reach the verdict and update on one half, and
test the updated (frozen) model on the other half with this function.

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

[`localize.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.gof.md)
for a fitted model checked on its own data;
[`deepgof1.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.external.md)
and
[`run.all.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.external.md)
for one overall test.

## Examples

``` r
set.seed(1)
n <- 500
X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)              # the published model
y <- rbinom(n, 1, plogis(qlogis(p) + 0.8 * X$x1 * X$x2))  # new patients: a missed interaction
localize.external(y, p, X, M = 199, seed = 1)   # M = 199 to keep the example fast
#> 
#> Localization of misfit: closed testing over orthogonal groups
#> 
#> data:     y and p (n = 500)
#> setting:  external validation of frozen predictions
#> reference: Monte Carlo, 199 draws
#> 
#> p-values by group (alone, and adjusted by closed testing):
#>            alone adjusted named
#> INTERCEPT 0.1950   0.5400      
#> SLOPE     0.8350   0.8350      
#> LINK      0.0050   0.0050     *
#> COV       0.0050   0.0050     *
#> 
#> groups named (closed testing) at FWER 0.05: LINK + COV
#> recommended update: revise the model; recalibration cannot repair misfit among
#>   patients who share a score
#> 
```
