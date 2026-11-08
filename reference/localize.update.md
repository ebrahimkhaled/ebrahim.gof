# Which update does a risk model need? Localization with the updated-model reference

Checks frozen predictions on validation data and names the parts of the
misfit that are present, so that the verdict points to the update the
model needs. The four parts and their bases are those of
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md);
what changes is the reference. Each part is tested against the model
already updated for the simpler parts:

|             |                                                  |
|-------------|--------------------------------------------------|
| `INTERCEPT` | against the frozen model \\p\\ (Monte Carlo)     |
| `SLOPE`     | against the intercept update \\a + \eta\\        |
| `LINK`      | against the logistic recalibration \\a + b\eta\\ |
| `COV`       | against the flexible recalibration `ns(eta, 4)`  |

with \\\eta = \mathrm{logit}(p)\\. Except for `INTERCEPT` this is a
parametric bootstrap: `B` outcome vectors are drawn from the fitted
update, the update is refitted to each draw, and each part's bases are
projected off that update's own directions with its own fitted weights.
The hypothesis of each test is therefore "this part is zero, whatever
the simpler parts are": a large error in overall risk or slope is
absorbed by the update and cannot inflate the tests of the more
demanding parts. This makes `localize.update()` the function to use
after transport to another population, where
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md)
can name a rung that is too high.

## Usage

``` r
localize.update(
  y,
  p,
  X,
  B = 199L,
  alpha = 0.05,
  naming = c("holm", "bonferroni"),
  cov_df = 5,
  cluster = NULL,
  tau = 0,
  main = TRUE,
  seed = NULL,
  plot = interactive()
)
```

## Arguments

- y:

  0/1 outcomes.

- p:

  the frozen predicted probabilities, one per outcome, in \\\[0, 1\]\\;
  values of exactly 0 or 1 are moved to `1e-10` and `1 - 1e-10`.

- X:

  a numeric matrix or data frame of covariates, one row per outcome,
  without missing values or constant columns. At least one column; `COV`
  needs two.

- B:

  number of reference draws for each part; each p-value lies on a grid
  of `1 / (B + 1)`.

- alpha:

  familywise level.

- naming:

  `"holm"` (the default) or `"bonferroni"`.

- cov_df:

  as in
  [`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md).

- cluster:

  optional vector of centre labels, one per outcome.

- tau:

  the standard deviation of the centre effect on the logit scale; used
  with `cluster`.

- main:

  if `TRUE`, `COV` includes the main-effects member; see Details.

- seed:

  optional integer, passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the
  reference draws. The random number stream of the session is restored
  on exit.

- plot:

  if `TRUE`, draw the verdict with
  [`plot.gof_localize`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/plot.gof_localize.md).
  The default draws it in an interactive session and not in scripts,
  simulations or examples.

## Value

An object of class `"gof_localize"`, as described in
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md).
Here `intersection` holds the Bonferroni test of each set of parts,
whose closure is Holm's procedure; `draws` is the smallest number of
usable draws over the parts (a refit that fails is dropped).

## Details

With this reference `COV` has a fourth member beside those of
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md):
the main effects of the covariates, so that a changed coefficient, of a
binary covariate too, has a test of its own. Risks pooled by
cross-validation are not fully frozen: their coefficients differ from
fold to fold, and the main-effects member responds to that difference
even when the form of the model is right. Use `main = FALSE` for such
risks, and check the error rate by repeating the refits in simulated
data.

The parts are named by Holm's correction of their p-values (the default)
or by Bonferroni's.

When the validation data come from several centres, outcomes within a
centre are correlated and the plain reference names `INTERCEPT` too
often. Give `cluster` and `tau`, the standard deviation of a centre
random effect on the logit scale, usually estimated by a
random-intercept model; every reference draw then carries a new effect
for each centre.

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

Ebrahim, E. K., Hussein, O. A. E.-A. and El-Kotory, A. (2026). Which
update does a clinical risk model need? Localization of misfit in
logistic risk models. Submitted. Reproduction materials:
[doi:10.5281/zenodo.24025418](https://doi.org/10.5281/zenodo.24025418)

Holm, S. (1979). A simple sequentially rejective multiple test
procedure. *Scandinavian Journal of Statistics*, 6(2), 65–70.

Vergouwe, Y., Nieboer, D., Oostenbrink, R., Debray, T. P. A., Murray, G.
D., Kattan, M. W., Koffijberg, H., Moons, K. G. M. and Steyerberg, E. W.
(2017). A closed testing procedure to select an appropriate method for
updating prediction models. *Statistics in Medicine*, 36(28), 4529–4539.
[doi:10.1002/sim.7179](https://doi.org/10.1002/sim.7179)

## See also

[`localize.size`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.size.md)
for the size of each part;
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md)
for the frozen-model reference;
[`localize.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.gof.md)
for a fitted model checked on its own data.

## Examples

``` r
set.seed(1)
n <- 500
X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)                     # the published model
y <- rbinom(n, 1, plogis(-0.7 + qlogis(p) + 0.8 * X$x1 * X$x2))  # new patients
localize.update(y, p, X, B = 49, seed = 1)   # B = 49 to keep the example fast; use 199 or more
#> 
#> Localization of misfit: updated-model reference
#> 
#> data:     y and p (n = 500)
#> setting:  external validation of frozen predictions
#> reference: updated model, 49 draws
#> 
#> p-values by group (alone, and adjusted by Holm over the groups):
#>            alone adjusted named
#> INTERCEPT 0.0200   0.0800      
#> SLOPE     0.5000   0.5000      
#> LINK      0.0200   0.0800      
#> COV       0.0200   0.0800      
#> 
#> groups named (Holm over the groups) at FWER 0.05: none
#> recommended update: keep the model; no misfit detected (the verdict is limited
#>   by the sample size)
#> 
```
