# ebrahim.gof 2.10.0

CRAN has 2.9.0 (published 2026-10-01).

## Why this release, and why it follows 2.9.0 so soon

I am aware that this update comes only ten days after 2.9.0 and that the package has been updated
often this year. I would not send it so soon if it did not correct test results that users rely
on. The first item below is a correctness fix:

* `run.all.external()`: the `Stukel (offset)` row used the likelihood-ratio statistic, which
  rejects too often at low event rates when only a few records have predicted risks above 0.5
  (type I error 0.061 to 0.078 at a nominal 0.05 in 4,000-replicate null simulations). It is now
  the score test of the same two terms at the null, which holds its level (0.035 to 0.053) and
  needs no model fit. A new `Osius-Rojek (external)` row uses the exact moments of the Pearson
  statistic under frozen predictions.
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
    Number of updates in past 6 months: 7

The reason for this update so soon after 2.9.0 is given at the top: it corrects the type I error
of a test in run.all.external().

    checking for future file timestamps ... NOTE
    unable to verify current time

This is the local machine being unable to reach the time server it uses, not a property of the
package.

## Reverse dependencies

None.
