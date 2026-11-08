# Multi-penalty permutation selection of edges

Selects the pairs whose ridge partial correlation exceeds a permutation
threshold at some penalty in a grid. For each penalty the threshold is
an order statistic of the maximal absolute partial correlation over the
observed data and `B` data sets in which each column is permuted
independently, at level `alpha` divided by the number of penalties.

## Usage

``` r
mps_select(x, mult = c(0.3, 1, 3, 10), B = 99L, alpha = 0.1)
```

## Arguments

- x:

  numeric data matrix with observations in rows.

- mult:

  penalty grid in units of \\\mathrm{tr}(S) / p\\.

- B:

  number of permutations.

- alpha:

  level in \\(0, 1)\\.

## Value

A list with the logical `adjacency` matrix, the penalties `lambda` and
the thresholds `threshold`.

## Details

If the rows of `x` are independent and identically distributed and its
columns are mutually independent, the probability that any edge is
selected is at most `alpha`, for every sample size and dimension. No
edge can be selected unless `B + 1 >= length(mult) / alpha`.

## References

Arashi, M. and Ebrahim, E. K. (2026). What does a ridge network
estimate? Effective penalty and leakage in ridge precision matrix
estimation. Submitted.

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(60 * 8), 60, 8)
x[, 2] <- x[, 1] + rnorm(60, sd = 0.5)
sel <- mps_select(x, B = 49)
which(sel$adjacency & upper.tri(sel$adjacency), arr.ind = TRUE)
#>      row col
#> [1,]   1   2
```
