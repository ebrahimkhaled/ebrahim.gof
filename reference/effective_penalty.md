# Effective penalty of the ridge precision estimator

Estimates the effective penalty \\\kappa = \lambda (1 + \delta)\\ of the
ridge precision estimator through the degrees-of-freedom ratio \\\hat r
= \mathrm{tr}\\S (S + \lambda I)^{-1}\\ / (n - 1)\\, as \\\hat\kappa =
\lambda / (1 - \hat r)\\; returns the limit of \\\hat\kappa\\ as
\\\lambda \to 0\\; and finds the penalty that attains a prescribed
effective penalty. The columns of `x` are centred and \\S = X^\top X /
(n - 1)\\, so the effective sample size in all formulas is \\n - 1\\,
where \\n\\ is the number of rows.

## Usage

``` r
effective_penalty(x, lambda)

effective_penalty_floor(x, tol = 1e-10)

penalty_for_kappa(x, kappa)
```

## Arguments

- x:

  numeric data matrix with observations in rows.

- lambda:

  vector of positive penalties.

- tol:

  relative tolerance below which an eigenvalue of \\S\\ is treated as
  zero.

- kappa:

  target effective penalty, larger than `effective_penalty_floor(x)`.

## Value

`effective_penalty` returns a data frame with columns `lambda`,
`df_ratio` (\\\hat r\\), `scale` (\\1 - \hat r\\, an estimate of \\1 /
(1 + \delta)\\) and `kappa`. `effective_penalty_floor` returns the
harmonic mean of the positive eigenvalues of \\S\\ when their number
equals the effective sample size, and zero otherwise.
`penalty_for_kappa` returns a penalty.

## Details

The estimated effective penalty is continuous and strictly increasing in
\\\lambda\\, so `penalty_for_kappa` has a unique solution.

## References

Arashi, M. and Ebrahim, E. K. (2026). What does a ridge network
estimate? Effective penalty and leakage in ridge precision matrix
estimation. Submitted.

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(30 * 60), 30, 60)
effective_penalty(x, c(0.01, 0.1, 1))
#>   lambda  df_ratio       scale    kappa
#> 1   0.01 0.9914974 0.008502633 1.176106
#> 2   0.10 0.9246945 0.075305543 1.327923
#> 3   1.00 0.6074250 0.392575010 2.547284
effective_penalty_floor(x)
#> [1] 1.15844
```
