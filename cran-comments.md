## Submission

This is a small feature update of an existing CRAN package (2.1.0 -> 2.1.1). It adds:

* `gof_install_suggests()`, a helper that installs the optional ("Suggests")
  packages used by the slow tests in `run.all.gof()`.
* an `install` argument to `run.all.gof()` that, when a selected test needs an
  optional package that is not installed, offers to install it.

There are no breaking changes.

## Note on installing Suggests (use of utils::install.packages)

`gof_install_suggests()` calls `utils::install.packages()`, but only with the
user's explicit consent, and never during checks:

* It is reached either by the user calling `gof_install_suggests()` directly, or
  by `run.all.gof(..., install = "ask"/"yes")`.
* It only installs after an interactive `readline()` confirmation (the default),
  unless the user explicitly passes `ask = FALSE`.
* Every install path is guarded by `interactive()`, so in a non-interactive
  session (scripts, `R CMD check`, the CRAN machines) nothing is ever installed
  and the library is never modified. The examples that would install are wrapped
  in `\dontrun{}`.

The library is therefore never modified without the user's explicit action, in
keeping with the CRAN policy on writing to the user's filespace.

## Test environments

* local: Windows 11, R 4.6.0 (ucrt)
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
