## Release summary

ebrahim.gof 3.1.0. CRAN currently has 3.0.0 (published 2026-11-01).

Sorry for another update so soon after 3.0.0. The main reason is a bug fix: in 3.0.0, `cavity.gof()` with
`recalibration = "logistic-normal"` could stop with "system is computationally singular" when the noise
parameter of the recalibration was weakly identified, which happens mostly under a strong misfit. The
covariance of the grouped residuals now uses a pseudo-inverse in that case, so the test returns both
p-values. Results are unchanged whenever 3.0.0 ran. Two manuscripts now being submitted describe the
functions of this release and cite version 3.1.0.

The release also adds exported functions; with their defaults, no existing function changes its results.

* New `localize.update()` and `localize.size()`: which update a risk model needs on validation data, and the
  size of each part of the misfit (Ebrahim, Hussein and El-Kotory, 2026, submitted).
* New ridge network functions, from the package ridgegraph 0.1.1 by Mohammad Arashi, who joins as an
  author: `ridge_precision()`, `ridge_pcor()`, `effective_penalty()`, `effective_penalty_floor()`,
  `penalty_for_kappa()`, `extrapolated_ridge()` and `mps_select()`, with a new vignette. No new dependency.

## Test environments

* local Windows 11, R 4.4.1

## R CMD check results

0 errors | 0 warnings | 2 notes

* "Days since last update: 6", "Number of updates in past 6 months: 9": see the release summary above.
  After this release we plan no further update for several months.
* "checking for future file timestamps ... unable to verify current time": the check machine could not
  reach the time server.

## Reverse dependencies

None.
