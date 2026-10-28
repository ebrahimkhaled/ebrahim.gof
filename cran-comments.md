## Release summary

ebrahim.gof 3.0.0. CRAN currently has 2.10.0. This release adds one exported function and one
optional argument; no existing function changes its results.

* New `pairscore.gof()`: the pair score test, a score test of whether a fitted logistic model needs
  the squares and the product of the two covariates that `deepgof1()` selects, calibrated by the same
  parametric bootstrap (the pair is chosen again in every bootstrap sample).
* `deepgof1()` gains `score.test = FALSE`; with `TRUE` it also returns the pair score test from the
  same refits. With the default, every output is unchanged (tested).

## Test environments

* local Windows 11, R 4.4.1

## R CMD check results

0 errors | 0 warnings | 1 note

* "unable to verify current time" (the check machine was offline for that step).
* "Found the following (possibly) invalid DOIs: 10.1002/sim.943": the DOI is valid (Crossref resolves it);
  the publisher's site refuses automated requests.
* "Number of updates in past 6 months": sorry for the quick succession. 2.10.0 was a resubmission that
  fixed a Linux test tolerance; this release adds the function that accompanies a manuscript now under
  submission, so that readers can run it from CRAN. We do not plan another update for some months.

## Reverse dependencies

None.
