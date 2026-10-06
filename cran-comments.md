<!-- On the day of submission: set Date in DESCRIPTION to that day and rebuild the tarball
     (R CMD build) on that day, or the incoming check notes that the Date field and the build
     time stamp are over a month old. -->

# ebrahim.gof 2.10.0

CRAN has 2.9.0 (published 2026-10-01).

## Why this release

* `ef.gof()` now computes the Ebrahim-Farrington test as its paper defines it (equal-size groups,
  tied risks in random order, the Theorem-3 normal reference). Its default p-value changes for some
  data; NEWS.md lists each change and gives the call that reproduces a 2.9.0 p-value.
* New functions for locating misfit: `localize.gof()` and `localize.external()`, with a `plot()`
  method.
* New functions `cubic.calib.gof()`, `edge.belt()` and `gof_certificate()`, and new methods
  `plot()` for `edge.gof()` and `deepgof1()` results and `summary()` for the battery.
  `edge.gof()` returns a data frame with the extra class `edge_gof`; its columns and values are
  unchanged.

None adds a dependency.

## Test environments

* local: Windows 11 x64 (build 26200), R 4.4.1 (2024-06-14 ucrt), `R CMD check --as-cran` on the
  built tarball -- 0 errors, 0 warnings, 2 notes (see below)

## R CMD check results

0 errors | 0 warnings | 2 notes

    checking CRAN incoming feasibility ... NOTE
    Days since last update: 5
    Number of updates in past 6 months: 7

The note will read differently on the day of submission.

    checking for future file timestamps ... NOTE
    unable to verify current time

This is the local machine being unable to reach the time server it uses, not a property of the
package.

## Reverse dependencies

None.
