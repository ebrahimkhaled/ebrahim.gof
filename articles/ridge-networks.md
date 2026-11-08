# Ridge networks: effective penalty and edge selection

The ridge precision estimator $`(S + \lambda I)^{-1}`$ has no zero
entries. For Gaussian data its partial correlations estimate those of
$`(\Sigma + \kappa I)^{-1}`$, where the effective penalty $`\kappa`$
exceeds $`\lambda`$ (Arashi and Ebrahim, 2026, submitted). This vignette
shows how to estimate $`\kappa`$ and how to select edges. These
functions come from the package ridgegraph by Mohammad Arashi.

``` r
library(ebrahim.gof)
set.seed(1)
p <- 30; n <- 20
theta <- diag(p)
for (i in 1:(p - 1)) theta[i, i + 1] <- theta[i + 1, i] <- 0.4
x <- matrix(rnorm(n * p), n, p) %*% chol(solve(theta))
effective_penalty(x, c(0.05, 0.5, 5))
#>   lambda  df_ratio      scale    kappa
#> 1   0.05 0.9232675 0.07673254 0.651614
#> 2   0.50 0.6788698 0.32113025 1.557001
#> 3   5.00 0.2748789 0.72512113 6.895400
effective_penalty_floor(x)
#> [1] 0.4619051
```

With $`p > n`$ the effective penalty does not fall below the floor
however small $`\lambda`$ is. Edges can be selected with the permutation
threshold.

``` r
sel <- mps_select(x, B = 99)
sum(sel$adjacency[upper.tri(sel$adjacency)])
#> [1] 1
```
