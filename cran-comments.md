## Submission

This is a major update of an existing CRAN package (1.0.0 -> 2.0.0). It adds
the Directed Ebrahim-Farrington test (`def.gof()`), an ensemble combiner
(`def.ensemble.gof()`), and a one-call battery runner (`run.all.gof()`).

It contains one intentional **breaking change**: `ef.gof()` now defaults to
`method = "chisq"` (referring the statistic to a chi-square_{G-2} distribution)
instead of the standardized-normal reference used in 1.0.0. The previous
behaviour remains available via `method = "normal"`. This is documented in
NEWS.md and the function help.

## Test environments

* local: Windows 11, R 4.4.1 (2024-06-14 ucrt)
* (planned before submission) win-builder: R-devel and R-release
* (planned before submission) GitHub Actions R-CMD-check: Windows / macOS / Ubuntu

## R CMD check results

0 errors | 0 warnings | 1 note

* The only NOTE is "checking for future file timestamps ... unable to verify
  current time". This is a local environment issue (the check machine cannot
  reach a time server) and does not occur on the CRAN/win-builder test machines.

## Dependencies

All packages used only inside individual tests are in Suggests and are accessed
conditionally via `requireNamespace()` (CompQuadForm, statmod, mgcv, BAGofT,
ResourceSelection). The package imports only base/recommended packages (stats).

## Reverse dependencies

There are no reverse dependencies (this package has no downstream packages on
CRAN).
