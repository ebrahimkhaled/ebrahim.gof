<!-- On the day of submission: set Date in DESCRIPTION to that day and rebuild the tarball
     (R CMD build) on that day, or the incoming check notes that the Date field and the build
     time stamp are over a month old. -->

# ebrahim.gof 2.9.0

CRAN has 2.8.0 (published 2026-09-26). This release follows it after a few days because it
corrects two bugs in `deepgof1()` that make the function return a wrong p-value.

## Why this release

`deepgof1()` computes its p-value from a parametric bootstrap. Up to 2.8.0 each bootstrap sample
was refitted by evaluating the model formula on the model frame. When the formula transforms a
covariate, for example `log(x)`, `splines::ns(x, 3)` or `poly(x, 2)`, the model frame holds the
transformed column and not `x`, so every refit failed. Every replicate was then scored `+Inf`, and
the p-value was 1 whatever the data, so the test could never reject such a model and gave no
warning. The refits now use `glm.fit()` on the fitted model's design matrix. For models without
such terms the p-value for a given seed is the same as in 2.8.0 (checked on the same seeds).

`deepgof1()` also refitted probit and complementary log-log models with the logit link; the refits
now use the model's own family and link (logit fits are unchanged). It now stops with a clear
message for a grouped (`cbind(successes, failures)`) or weighted binomial fit, which its Bernoulli
bootstrap never served correctly.

## New functions

* `deepgof1.external()`: the network test for frozen predictions (external validation of a fixed
  model), with an exact Monte Carlo p-value.
* `run.all.external()`: the external-validation tests on outcomes and frozen predicted
  probabilities, in the battery format of `run.all.gof()`. GiViTI runs only when the suggested
  'givitiR' and 'callr' are installed, as in `run.all.gof()`.
* `edge.stream()` with `update()`, `summary()` and `print()` methods: calibration monitoring of a
  frozen model on streaming data, in constant time per record.

None adds a dependency.

## New options

* `edge.gof(external = TRUE)` and `def.gof(external = TRUE)`: the directed test for frozen
  predictions. The default is `FALSE`.
* `run.all.gof()` gains `HL-largeN`, the large-sample Hosmer-Lemeshow test of Nattino, Pennell and
  Lemeshow (2020, Biometrics); it reproduces their published example.
* `deepgof1(reading = "allpairs")` and `reading = "combined"`; the default reading is unchanged in
  result for models whose covariates enter as single untransformed columns.

## Behaviour changes

* `edge.gof()` now reports two rows by default, the default partition (the verdict) and ten groups
  (a robustness check), with new `Partition` and `Role` columns. `edge.gof(fit, G = 10)` reproduces
  the earlier single-row result. `def.gof()` is unchanged.
* `G = "auto"` in the EDGE functions now uses `max(10, ceiling(n / 25))` groups, the published rule,
  instead of `round(n / 25)`; the number of groups moves by at most one.
* `deepgof1()` returns a different (correct) p-value for models with transformed covariates or a
  non-logit link, and an error instead of a p-value for grouped or weighted binomial fits.

## Test environments

* local: Windows 11 x64 (build 26200), R 4.4.1 (2024-06-14 ucrt), `R CMD check --as-cran` on the
  built tarball -- 0 errors, 0 warnings, 2 notes (see below)

## R CMD check results

0 errors | 0 warnings | 2 notes

    checking CRAN incoming feasibility ... NOTE
    Days since last update: 5

This update corrects the `deepgof1()` bugs described above; the first silently returns p = 1 for
any model with a transformed covariate.

    checking for future file timestamps ... NOTE
    unable to verify current time

This is the local machine being unable to reach the time server it uses, not a property of the
package.

## Reverse dependencies

None.
