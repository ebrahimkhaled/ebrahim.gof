## Release summary

This submission also fixes the test failures on CRAN's MKL check (issues reported 2026-10-29 for
2.10.0, to be corrected before 2026-11-21). They came from tests comparing two separately computed
floating-point results with identical(), which MKL does not reproduce bit for bit; these now use
all.equal() tolerance, and the tests that reproduce the BAGofT package exactly run off CRAN only.
No package code changed for this.

ebrahim.gof 3.0.0. CRAN currently has 2.10.0. This release adds two exported functions and
optional arguments to two others; with their defaults, no existing function changes its results.

* New `pairscore.gof()`: the pair score test, a score test of whether a fitted logistic model needs
  the squares and the product of the two covariates that `deepgof1()` selects, calibrated by the same
  parametric bootstrap (the pair is chosen again in every bootstrap sample).
* `deepgof1()` gains `score.test = FALSE`; with `TRUE` it also returns the pair score test from the
  same refits. With the default, every output is unchanged (tested).
* `legoft()` gains `panel` and `rule`: alternative panels of tests and combination rules, with the
  same bootstrap calibration. The default is the existing rule, unchanged (tested).
* New `cavity.gof()`: a leave-one-out Hosmer-Lemeshow test for ridge, lasso and elastic-net logistic
  regression in high dimensions (Ebrahim and Arashi, 2026, submitted). `glmnet` is added to Suggests
  and its tests are guarded by skip_if_not_installed().

## Test environments

* local Windows 11, R 4.4.1

## R CMD check results

0 errors | 0 warnings | 2 notes

* "Number of updates in past 6 months": sorry for the quick succession. This submission fixes the
  MKL test failures CRAN reported for 2.10.0 and adds the functions that accompany three manuscripts
  now under submission, in one release, so that readers can run them from CRAN. We do not plan
  another update for some months.
* "unable to verify current time" (the check machine was offline for that step). Some DOIs in the
  Rd files timed out from the same machine; each resolves at doi.org.

## Reverse dependencies

None.
