<!-- On the day of submission: set Date in DESCRIPTION to that day and rebuild the tarball
     (R CMD build) on that day, or the incoming check notes that the Date field and the build
     time stamp are over a month old. -->

# ebrahim.gof 2.8.0

## Why this release

It corrects a statistic. The `Stukel` row of `run.all.gof()` squared and summed two marginal score
statistics and referred the sum to chi-squared on 2 degrees of freedom. After fitting, the two
statistics are correlated, so the sum is not chi-squared on 2 degrees of freedom and the test was
liberal: in simulation it rejected up to about 7 per cent of correctly specified models at the 5 per
cent level. When every fitted risk lay on one side of 0.5 it returned `NaN` without a note. The row
now reports the joint score test for the same two directions, on 2 degrees of freedom, or on 1 with a
note when every fitted risk is on one side of 0.5. This is why 2.8.0 follows 2.7.0 closely.

The same release fixes `def.gof()` and `edge.gof()`, which stopped with a singular-system error when
one group's mean fitted risk lay just above 0.5 with `basis = "stukel"`.

## New options

* `run.all.gof(control = list(Stukel = list(form = ...)))`: `"lr"` gives the likelihood-ratio test
  for the same two directions, and `"marginal"` the old statistic, kept for reproducing earlier
  results.
* `def.gof()`, `edge.gof()` and `def.ensemble.gof()`: the basis `"sym"`, `weights = "score"` and
  `G = "auto"`. The defaults are unchanged.

## Behaviour changes

* The default `Stukel` statistic in `run.all.gof()` changes, so its p-values differ from 2.7.0.
* `def.gof()`, `edge.gof()` and `def.ensemble.gof()` now warn when there are fewer events, or fewer
  non-events, than groups. The p-value is still returned. Inside `run.all.gof()` the message goes to
  the `Note` column instead of a warning.
* The fast battery of `run.all.gof()` has one more row, `DEF.sym`.
* On a sample with no event, or no non-event, `def.gof()`, `edge.gof()` and `def.ensemble.gof()` return
  `NA` with a warning instead of stopping with an error, and the `Stukel` row returns `NA` with a note.
* The score form of `def.gof()` and the joint `Stukel` form leave out a column whose information after
  the fit is below 1e-10 of its information before the fit; with none left they return `NA`.

## Test environments

* local: Windows 11 x64 (build 26200), R 4.4.1 (2024-06-14 ucrt), `R CMD check --as-cran` on the
  built tarball -- 0 errors, 0 warnings, 2 notes (see below)
* win-builder, R-devel and R-release -- not yet run for 2.8.0; to be run on the tarball rebuilt on
  the day of submission
* GitHub Actions: Windows-release, macOS-release, Ubuntu devel/release/oldrel-1 -- not yet run for
  2.8.0

## R CMD check results

0 errors | 0 warnings | 2 notes

    checking CRAN incoming feasibility ... NOTE
    Maintainer: 'Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>'

    Days since last update: 5

2.7.0 was published on 2026-09-09. This update follows it closely because it corrects the liberal
Stukel statistic described above.

    checking for future file timestamps ... NOTE
    unable to verify current time

This is the local machine being unable to reach the time server it uses, not a property of the
package.

## Reverse dependencies

None.
