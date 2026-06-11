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
* win-builder: R-devel (R Under development, 2026-06-09 r90126 ucrt) -- Status: OK
* R-hub: Linux (R-devel) and Windows (R-devel) -- Status: OK
* GitHub Actions (R-CMD-check): Windows-release, macOS-release, and
  Ubuntu (R-devel, R-release, R-oldrel-1) -- all passing

## R CMD check results

0 errors | 0 warnings | 0 notes

* win-builder (R-devel) and R-hub return Status: OK with no notes.
* The local check shows one NOTE only -- "checking for future file timestamps
  ... unable to verify current time" -- which is a local environment issue (the
  machine cannot reach a time server). It does not appear on win-builder, R-hub,
  or the CRAN test machines.

## Dependencies

All packages used only inside individual tests are in Suggests and are accessed
conditionally via `requireNamespace()` (CompQuadForm, statmod, mgcv, BAGofT,
ResourceSelection). The package imports only base/recommended packages (stats).

## Reverse dependencies

There are no reverse dependencies (this package has no downstream packages on
CRAN).
