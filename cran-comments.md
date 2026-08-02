## Submission

This is a feature update of an existing CRAN package (2.1.0 -> 2.4.0). It
reframes the package as a unified goodness-of-fit and calibration *toolbox* for
binary logistic regression and adds the following, without breaking changes:

* `run.all.gof()`, a single-call battery that runs a wide range of classical,
  modern, and author tests and returns a tidy `gof_battery` object (with
  `print` and `plot` methods). Optional-package tests are accessed
  conditionally and are simply skipped when their package is not installed.
* `edges.gof()`, a brand-name alias for the existing `def.ensemble.gof()`
  (the EDGES Cauchy-combination ensemble of the directed EDGE bases), and
  `edge.gof()`, the primary interface to the directed EDGE test (the legacy
  name `def.gof()` is retained unchanged). Both aliases compute exactly the
  same statistics as the legacy names and are documented.
* a `parallel` argument (with `ncores`) to `run.all.gof()`, which optionally
  runs the slow resampling loops on a PSOCK cluster via `parallel::parLapply()`.
  It is `FALSE` by default; the parallel path is reproducible via
  `parallel::clusterSetRNGStream()`.
* two small example datasets, `gof_demo` and `gof_demo_grouped`, used in the
  examples and vignette.

There are no breaking changes; all previously exported functions keep their
names and behaviour.

## Note on installing Suggests (use of utils::install.packages)

The optional helper `gof_install_suggests()` (and `run.all.gof(..., install =)`)
calls `utils::install.packages()`, but only with the user's explicit consent,
and never during checks:

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

## Startup message

The package emits a single, conditional startup message via
`packageStartupMessage()` in `.onAttach()`: only when one or more of the
optional `Suggests` packages is not installed, it names them and points to
`gof_install_suggests()`. When all optional packages are present it is silent,
and it is fully suppressible with `suppressPackageStartupMessages()`. It never
installs or writes anything.

## Test environments

* local: Windows 11, R 4.4.1, `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

* The only NOTE is "checking for future file timestamps ... unable to verify
  current time", which is a local/offline environment issue (the check machine
  cannot reach a time server). It is not a property of the package.

## Dependencies

All packages used only inside individual tests are in Suggests and are accessed
conditionally via `requireNamespace()` (CompQuadForm, statmod, mgcv, BAGofT,
ResourceSelection, givitiR, callr). The GiViTI test additionally runs `givitiR`
in a `callr` subprocess so that a crash in its compiled dependencies cannot
abort the user's session. The package imports only base/recommended packages
(parallel, stats).

## Reverse dependencies

There are no reverse dependencies (this package has no downstream packages on
CRAN).
