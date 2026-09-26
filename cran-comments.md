<!-- On the day of submission: set Date in DESCRIPTION to that day and rebuild the tarball
     (R CMD build) on that day, or the incoming check notes that the Date field and the build
     time stamp are over a month old. -->

# ebrahim.gof 2.8.0

CRAN has 2.7.0 (published 2026-09-09). No version was prepared in between, so this release carries
the changes listed below and nothing else.

## Why this release

It corrects two statistics.

* The `Stukel` row of `run.all.gof()` squared and summed two marginal score statistics and referred
  the sum to chi-squared on 2 degrees of freedom. After fitting, the two statistics are correlated,
  so the sum is not chi-squared on 2 degrees of freedom and the test was liberal: in simulation it
  rejected up to about 7 per cent of correctly specified models at the 5 per cent level. When every
  fitted risk lay on one side of 0.5 it returned `NaN` without a note. The row now reports the joint
  score test for the same two directions, on 2 degrees of freedom, or on 1 with a note when every
  fitted risk is on one side of 0.5.
* `deepgof1()` placed tied covariate values into the cells of its residual map by row order, so on
  data sorted by the outcome a correct model could be rejected (p = .005 in file order against a
  median of .54 over random row orders on a published data set). Ties are now broken at random,
  once per call, and the same order is used for the observed and every bootstrap map.

The same release fixes `def.gof()` and `edge.gof()`, which stopped with a singular-system error when
one group's mean fitted risk lay just above 0.5 with `basis = "stukel"`.

## New functions and options

* `projection.gof()`: the projection test of Escanciano (2006) for logistic regression, as defined
  by Liu et al. (2024, Statistics and Computing). Also a slow row, `Projection`, of the full battery.
* `bagoft.fast()`: the BAGofT test of Zhang, Ding and Yang (2023, JASA), returning the same p-values
  as the BAGofT package 1.0.0 for the same seed, at less than half the time. The `BAGofT` row of the
  battery now uses it when `randomForest` is installed.
* `run.all.gof(control = list(Stukel = list(form = ...)))`: `"lr"` gives the likelihood-ratio test
  for the same two directions, and `"marginal"` the old statistic, kept for reproducing earlier
  results.
* `def.gof()`, `edge.gof()` and `def.ensemble.gof()`: the basis `"sym"`, `weights = "score"` and
  `G = "auto"`. The defaults are unchanged.

## Behaviour changes

* The default `Stukel` statistic in `run.all.gof()` changes, so its p-values differ from 2.7.0.
* `deepgof1()` gives a different p-value for a given seed than 2.7.0 when some covariate has tied
  values; without ties the p-value is unchanged.
* `def.gof()`, `edge.gof()` and `def.ensemble.gof()` now warn when there are fewer events, or fewer
  non-events, than groups. The p-value is still returned. Inside `run.all.gof()` the message goes to
  the `Note` column instead of a warning.
* The fast battery of `run.all.gof()` has one more row, `DEF.sym`; the full battery also gains
  `Projection`.
* On a sample with no event, or no non-event, `def.gof()`, `edge.gof()` and `def.ensemble.gof()` return
  `NA` with a warning instead of stopping with an error, and the `Stukel` row returns `NA` with a note.
* The score form of `def.gof()` and the joint `Stukel` form leave out a column whose information after
  the fit is below 1e-10 of its information before the fit; with none left they return `NA`.

## Notes for the CRAN team

* `R/bagoft_fast.R` adapts code from the CRAN package BAGofT (GPL-3, the same licence as this
  package). Its three authors, Jiawei Zhang, Jie Ding and Yuhong Yang, are listed in `Authors@R` as
  contributors and copyright holders.
* Two packages were added to Suggests, `randomForest` and `dcov`, both on CRAN; they are used only
  by `bagoft.fast()` and the `BAGofT` row, and are checked with `requireNamespace()`. The tests that
  need them are skipped when they are absent.
* The example of `bagoft.fast()` and the full-size example of `projection.gof()` are in
  `\donttest{}` because the random forests and the bootstrap take longer than 5 seconds at useful
  settings; `projection.gof()` also has a small example that runs.

## Test environments

* local: Windows 11 x64 (build 26200), R 4.4.1 (2024-06-14 ucrt), `R CMD check --as-cran` on the
  built tarball -- 0 errors, 0 warnings, 1 note (see below)
* win-builder, R-devel and R-release -- not yet run for 2.8.0; to be run on the tarball rebuilt on
  the day of submission
* GitHub Actions: Windows-release, macOS-release, Ubuntu devel/release/oldrel-1 -- not yet run for
  2.8.0

## R CMD check results

0 errors | 0 warnings | 1 note

The incoming-feasibility check reported only the maintainer line. 2.7.0 was published on
2026-09-09; this update follows it within three weeks because it corrects the liberal Stukel
statistic and the row-order dependence of `deepgof1()` described above.

The tarball grows from 867 KB (2.7.0) to 905 KB: the two new functions, their help pages and tests.

    checking for future file timestamps ... NOTE
    unable to verify current time

This is the local machine being unable to reach the time server it uses, not a property of the
package.

## Reverse dependencies

None.
