# ebrahim.gof 2.7.0

## What is new

`calm.gof()`, a closed-form reference distribution for the shrinkage-corrected goodness-of-fit
statistics. It is the companion of `shrink.gof()`, released in 2.5.0: the statistics are the same
and only the reference differs, so where `shrink.gof()` needs several hundred penalised refits
`calm.gof()` needs one fit. It accompanies a manuscript under review; the method and the designs in
which its decile variant is not reliable are documented in `?calm.gof`.

`CompQuadForm` moves from Suggests to Imports. `calm.gof()` computes the exact tail of a weighted
chi-squared law and cannot return a p-value without it, so it is no longer optional. No other
dependency changed.

## Test environments

* local: Windows 11, R 4.4.x, `R CMD check --as-cran`
* (fill in before submitting) win-builder devel and release
* (fill in before submitting) macOS builder

## R CMD check results

0 errors | 0 warnings | 1 note

The note is

    checking for future file timestamps ... NOTE
    unable to verify current time

which is the check machine being unable to reach the time server it uses, not a property of the
package. It does not appear on win-builder.

## Reverse dependencies

None.

## Notes for the CRAN team

The examples for `calm.gof()` run in about a second on a 300 by 20 design. The heavier comparisons
in the manuscript are not run by the examples or the vignettes.
