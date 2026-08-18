# ebrahim.gof 2.5.0

## Test environments

* Windows 11, R 4.4.1 (local), `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

The single note is

```
* checking for future file timestamps ... NOTE
  unable to verify current time
```

which is raised because the check machine could not reach the remote time service. It
does not occur on machines with network access.

## What is in this release

This release ships a bug fix that has been sitting unreleased, and three new functions.

**Bug fix (important for users comparing results across versions).** `le-Cessie` formed
its moment reference as `(I - H) R (I - H)`; the residual expansion requires
`(I - H)' R (I - H)`. In a weighted fit `H = VX(X'VX)^{-1}X'` is idempotent but not
symmetric, so the transpose matters -- the inherited implementation came from the
unweighted linear-model setting, where the two forms coincide. The raw quadratic form was
always correct, but the *reported* statistic is that form rescaled by the matched moments,
so the printed statistic, degrees of freedom and p-value all change. NEWS.md states the
before and after values on the package's own regression example, and the regression test
was updated to the corrected values with a comment recording that disagreeing with the
inherited source is the intent of the fix.

**New functions.** `legoft()` and `legoft.localize()` accompany a manuscript on combining
and localizing goodness-of-fit evidence for binary logistic regression; `shrink.gof()`
accompanies a manuscript on goodness of fit under ridge penalization. NEWS.md documents,
for `legoft()`, both what it does and what it does not do -- it attains higher mean power
than each individual test studied but does not dominate them family by family, and it does
not improve on an equal-weight Cauchy combination of the same members.

**Documentation.** The help for `run.all.gof()` was rewritten to describe each of the
tests individually, grouped by mechanism, with the failure mode of each stated; ten
classical tests that were implemented but uncited now carry references. This is why the
generated `run.all.gof.Rd` is substantially larger than in 2.4.1.

## Reverse dependencies

None. This package has no reverse dependencies on CRAN.

## Other notes for the maintainers

* Examples that need bootstrap replication are wrapped in `\donttest{}`; the full battery
  with `include_slow = TRUE` is exercised in tests that carry `skip_on_cran()`.
* `shrink.gof()` is vendored byte-identical from the manuscript source, so the paper and
  the package compute the same numbers; its `stats` imports are declared through roxygen
  rather than by editing the vendored code.
