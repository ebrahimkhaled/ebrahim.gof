## Submission

This is a feature update of an existing CRAN package (2.0.0 -> 2.1.0). It adds:

* New exported functions: `cdef.gof()` (covariate-space directed test),
  `gof.features()` (the GOF evidence vector), and `deploy.gof()` (a deployable
  learned-ensemble test).
* `run.all.gof()` gains the `McCullagh` test (exact conditional-moment
  standardization of the Pearson statistic) and the `GiViTI` polynomial
  calibration test, completing the recommended core. GiViTI wraps the suggested
  `givitiR` package run inside an isolated `callr` subprocess, so a failure in
  its compiled dependencies returns `NA` rather than aborting the session.
* `run.all.gof()` now returns a `gof_battery` object with `print` and `plot`
  methods (the plot draws the GiViTI calibration belt), `include_slow` defaults
  to `TRUE`, and a new `calibration_plot` argument is available.
* The `BAGofT` wrapper now also runs on single-predictor models, and all
  diagnostic `Note` messages were rewritten to clear phrases.

There are no breaking changes in this release.

## Test environments

* local: Windows 11, R 4.4.1 (2024-06-14 ucrt)
* win-builder: R-devel and R-release (to be confirmed on submission)
* GitHub Actions (R-CMD-check): Windows-release, macOS-release, and
  Ubuntu (R-devel, R-release, R-oldrel-1)

## R CMD check results

0 errors | 0 warnings | 1 note

* The only NOTE is "checking for future file timestamps ... unable to verify
  current time", a local environment issue (the machine cannot reach a time
  server). It does not appear on win-builder, R-hub, or the CRAN machines.

## Dependencies

All packages used only inside individual tests are in Suggests and are accessed
conditionally via `requireNamespace()` (CompQuadForm, statmod, mgcv, BAGofT,
ResourceSelection, givitiR, callr). The GiViTI test additionally runs `givitiR`
in a `callr` subprocess so that a crash in its compiled dependencies cannot
abort the user's session. The package imports only base/recommended packages
(stats, graphics).

## Reverse dependencies

There are no reverse dependencies (this package has no downstream packages on
CRAN).
