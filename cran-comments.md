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

* local: Windows 11, R 4.4.1, `R CMD check --as-cran` -- 1 note (see below)
* GitHub Actions: Windows-release, macOS-release, Ubuntu devel/release/oldrel-1 -- all passing
* win-builder, R-devel (2026-09-06 r90498 ucrt) -- 1 note (see below); install 5s, check 99s
* win-builder, R-release -- submitted, but the result mail has not been delivered to the
  maintainer address on three attempts; the R-release platform is covered by the GitHub
  Actions Windows-release job above, which passes

## R CMD check results

0 errors | 0 warnings | 1 note on each environment, and they are different notes.

On win-builder R-devel:

    checking CRAN incoming feasibility ... NOTE
    Possibly misspelled words in DESCRIPTION:
      prepivoting (29:26)

"Prepivoting" is the standard term for Beran's bootstrap transformation, from Beran (1987),
"Prepivoting to reduce level error of confidence sets", Biometrika 74(3), 457-468. It is spelled
correctly and is the name of the procedure `shrink.gof()` implements. It is listed in
`inst/WORDLIST`, which the `spelling` package honours but the incoming check's own aspell run does
not.

Locally:

    checking for future file timestamps ... NOTE
    unable to verify current time

which is this machine being unable to reach the time server it uses, not a property of the package.
It does not appear on win-builder.

## Reverse dependencies

None.

## Notes for the CRAN team

The examples for `calm.gof()` run in about a second on a 300 by 20 design. The heavier comparisons
in the manuscript are not run by the examples or the vignettes.
