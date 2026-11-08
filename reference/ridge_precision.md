# Ridge precision matrix and ridge partial correlations

Computes the ridge precision estimator \\(S + \lambda I)^{-1}\\ from a
data matrix and the corresponding partial correlations \\-q\_{ij} /
(q\_{ii} q\_{jj})^{1/2}\\.

## Usage

``` r
ridge_precision(x, lambda)

ridge_pcor(x, lambda)
```

## Arguments

- x:

  numeric data matrix with observations in rows. Columns are centred
  internally and \\S = X^\top X / (n - 1)\\.

- lambda:

  a single positive penalty.

## Value

A symmetric \\p \times p\\ matrix. For `ridge_pcor` the diagonal is set
to one.

## Details

The partial correlations do not depend on the scale of the estimator.
For Gaussian data they estimate the partial correlations of \\(\Sigma +
\kappa I)^{-1}\\, where \\\kappa\\ is the effective penalty returned by
[`effective_penalty`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/effective_penalty.md).

## References

Arashi, M. and Ebrahim, E. K. (2026). What does a ridge network
estimate? Effective penalty and leakage in ridge precision matrix
estimation. Submitted.

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(40 * 10), 40, 10)
r <- ridge_pcor(x, lambda = 0.5)
range(r[upper.tri(r)])
#> [1] -0.2142514  0.2108523
```
