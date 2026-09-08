# Shrinkage-corrected goodness-of-fit test for penalized (ridge) logistic regression

The Hosmer–Lemeshow test is not valid when the coefficients are shrunk:
penalization biases the fitted probabilities, the grouped residuals
acquire a non-centrality, and the usual chi-squared reference is wrong.
`shrink.gof()` removes that non-centrality and refers the corrected
statistic to a bootstrap built from the debiased generator. For the same
correction referred to a closed-form reference, needing one fit rather
than `B` of them, see
[`calm.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md).

## Usage

``` r
shrink.gof(
  X,
  y,
  lambda,
  G = 10,
  basis = c("edge", "decile"),
  B = 999,
  seed = NULL,
  penalize = NULL,
  uncorrected = FALSE
)
```

## Arguments

- X:

  numeric design matrix, without an intercept column.

- y:

  binary response (0/1) of length `nrow(X)`.

- lambda:

  ridge penalty on the **theory** scale, `lambda = n * lambda_glmnet`.
  Passing a glmnet lambda directly is the commonest error; multiply by
  `n`.

- G:

  number of groups.

- basis:

  `"edge"`, `"decile"`, or both.

- B:

  bootstrap replicates.

- seed:

  optional integer for reproducibility.

- penalize:

  logical vector marking which columns of `X` are penalized; the
  intercept is never penalized.

- uncorrected:

  if `TRUE`, also return the uncorrected statistic, which is what a
  naive application of Hosmer–Lemeshow to a penalized fit computes.

## Value

An object of class `"shrink.gof"`: a list carrying the settings the test
ran under (`lambda` on the theory scale, `G`, `B`, `n`, `p`) and, for
each requested basis, a component named `SC.HL` or `SC.EDGE` holding its
`statistic` and bootstrap `p.value`, plus `p.uncorrected` when
`uncorrected = TRUE`. Printed by `print.shrink.gof`.

## Details

The correction subtracts the estimated shrinkage non-centrality and
prepivots against \\\pi(\tilde\beta)\\ rather than \\\pi(\hat\beta)\\
(Beran prepivoting), so the reference world matches the null being
tested.

Reference implementation for: Ebrahim, E. K. (2026), "Shrinkage
invalidates the Hosmer–Lemeshow test: goodness of fit for penalized
logistic regression, with an application to glaucoma diagnosis."

Depends only on base R and stats, so results do not move with package
versions.

## Choosing between this and calm.gof

The statistics are the same; only the reference differs.
[`calm.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md)
reads the exact tail of a weighted chi-squared law from a single fit, so
it returns in a fraction of a second and its p-value carries no Monte
Carlo error – which matters when the p-value is read as a magnitude
rather than compared with a threshold. `shrink.gof()` needs `B`
penalized refits and its p-value is granular to \\1/(B+1)\\.

Prefer
[`calm.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md)
unless one of these applies: the columns to shrink are chosen through
`penalize`, which
[`calm.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md)
does not accept; or \\p \ge n\\, which it refuses. Note also that
[`calm.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md)
takes a `lambda_scale` argument and `shrink.gof()` does not, so the same
number passed to both is a `glmnet` penalty in one and a theory-scale
penalty in the other – a factor of \\n\\ apart. This function expects
the theory scale.

## See also

[`calm.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md),
which refers the same corrected statistics to a closed-form reference
and needs a single fit, and the section above for when to prefer it;
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
for the unpenalized battery.

## Examples

``` r
# \donttest{
set.seed(1)
X <- matrix(rnorm(400 * 5), 400, 5)
y <- rbinom(400, 1, plogis(0.3 + X %*% c(0.8, -0.5, 0.3, 0, 0)))
shrink.gof(X, y, lambda = 40, basis = "decile", B = 99, seed = 1)
#> 
#> Shrinkage-corrected Hosmer-Lemeshow test
#> n = 400, p = 5 (p/n = 0.013), lambda = 40, G = 10, B = 99
#> 
#>   SC.HL    statistic =   8.5932   p = 0.3000
#> 
# }
```
