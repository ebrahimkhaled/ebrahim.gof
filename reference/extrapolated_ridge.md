# Extrapolation of ridge precision estimates across penalties

Combines ridge estimates at two effective penalties \\\kappa_1\\ and
\\\kappa_2 = r \kappa_1\\ with weights \\\kappa_2 / (\kappa_2 -
\kappa_1)\\ and \\-\kappa_1 / (\kappa_2 - \kappa_1)\\, after rescaling
each estimate by its estimated scale. The population version of this
combination differs from the precision matrix by a term of second order
in the penalties.

## Usage

``` r
extrapolated_ridge(x, lambda, ratio = 2)
```

## Arguments

- x:

  numeric data matrix with observations in rows.

- lambda:

  first (smaller) penalty.

- ratio:

  ratio \\\kappa_2 / \kappa_1 \> 1\\ of the effective penalties.

## Value

A list with the combined matrix `precision`, its partial correlations
`pcor`, the two penalties `lambda`, the two estimated effective
penalties `kappa` and the `weights`.

## Details

The sum of the absolute weights is \\(r + 1) / (r - 1)\\, which
multiplies the sampling variability.

## References

Arashi, M. and Ebrahim, E. K. (2026). What does a ridge network
estimate? Effective penalty and leakage in ridge precision matrix
estimation. Submitted.

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(50 * 20), 50, 20)
fit <- extrapolated_ridge(x, lambda = 0.5)
fit$kappa
#> [1] 0.6662174 1.3324349
```
