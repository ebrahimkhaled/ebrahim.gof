<!-- On the day of submission: set Date in DESCRIPTION to that day and rebuild the tarball
     (R CMD build) on that day, or the incoming check notes that the Date field and the build
     time stamp are over a month old. -->

# ebrahim.gof 2.9.0

CRAN has 2.8.0 (published 2026-09-26). This release follows it after a few days because it
corrects a bug in `deepgof1()` that makes the function return a wrong p-value.

## Why this release

`deepgof1()` computes its p-value from a parametric bootstrap. Up to 2.8.0 each bootstrap sample
was refitted by evaluating the model formula on the model frame. When the formula transforms a
covariate, for example `log(x)`, `splines::ns(x, 3)` or `poly(x, 2)`, the model frame holds the
transformed column and not `x`, so every refit failed. Every replicate was then scored `+Inf`, and
the p-value was 1 whatever the data, so the test could never reject such a model and gave no
warning. The refits now use `glm.fit()` on the fitted model's design matrix. For models without
such terms the p-value for a given seed is the same as in 2.8.0 (checked on the same seeds).

`deepgof1()` also now stops with a clear message for a grouped (`cbind(successes, failures)`) or
weighted binomial fit, which its Bernoulli bootstrap never served correctly.

## New options

* `deepgof1(reading = "allpairs")`: the largest score over the residual maps of every pair of
  covariates, taken again in every bootstrap replicate. The default reading is unchanged.
* `deepgof1()` now accepts a model with one covariate (earlier versions stopped with an error), and
  its result also returns the map that gave the statistic.

## Behaviour changes

* `deepgof1()` returns a different (correct) p-value for models with transformed covariates, and an
  error instead of a p-value for grouped or weighted binomial fits. Nothing else changes by default.

## Test environments

* local: Windows 11 x64 (build 26200), R 4.4.1 (2024-06-14 ucrt), `R CMD check --as-cran` on the
  built tarball -- 0 errors, 0 warnings, 2 notes (see below)
* win-builder, R-devel and R-release -- to be run on the tarball rebuilt on the day of submission

## R CMD check results

0 errors | 0 warnings | 2 notes

    checking CRAN incoming feasibility ... NOTE
    Days since last update: 3

This update corrects the `deepgof1()` bug described above, which silently returns p = 1 for any
model with a transformed covariate.

    checking for future file timestamps ... NOTE
    unable to verify current time

This is the local machine being unable to reach the time server it uses, not a property of the
package.

## Reverse dependencies

None.
