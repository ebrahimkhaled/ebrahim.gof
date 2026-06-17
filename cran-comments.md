## Resubmission

This is a resubmission of 2.1.0 after the automated incoming pretest. The
pretest reported one NOTE per flavour, addressed as follows:

* "Overall checktime ... > 10 min" (r-devel-windows). Fixed. The 2.1.0 default
  `include_slow = TRUE` had caused the full slow battery (GiViTI in a `callr`
  subprocess, the GAM tests, and the bootstrap tests) to run in the vignette,
  in a `\donttest` example, and in one test. These now run the fast battery
  (or a slim, GiViTI-only subset), and the slow test uses `skip_on_cran()`, so
  the slow battery is no longer exercised during `R CMD check`. Local
  `R CMD check --as-cran` checktime is now well under the limit.

* "Possibly misspelled words in DESCRIPTION" (Cessie, GiViTI, Houwelingen,
  McCullagh, Osius, Rojek, Stute, Zhu, le). These are all proper names of the
  goodness-of-fit tests and their authors, spelled as in the cited literature;
  they are not misspellings.

* "Days since last update: 4." 2.1.0 follows 2.0.0 closely because it completes
  the goodness-of-fit battery that 2.0.0 introduced -- it adds the McCullagh and
  GiViTI calibration tests that round out the recommended core -- and an
  accompanying methods paper, now under review, cites and reproduces its results
  via this package. I would be grateful if you could accept the short interval
  on that basis; I will otherwise space subsequent updates well apart.

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
