# Projection Goodness-of-Fit Test for Binary Regression

`projection.gof()` computes the projection test of Escanciano (2006) as
defined for logistic regression by Liu et al. (2024, Sec. 2.1), and
refers it to their model-based bootstrap. The residual-marked empirical
process is taken along every direction of the covariate space and its
Cramer-von Mises norm is integrated over the unit sphere, so the test is
consistent against any departure of the mean function, including
departures that are invisible along the fitted linear predictor (where
the Stute-Zhu test looks).

## Usage

``` r
projection.gof(object, B = 1000, scale = FALSE, tol = 1e-10)
```

## Arguments

- object:

  A fitted binary [`glm`](https://rdrr.io/r/stats/glm.html)
  (`family = binomial`, any link) with unit prior weights and at least
  one covariate.

- B:

  Number of bootstrap replicates (default 1000, as in Liu et al.'s
  examples).

- scale:

  Logical; standardize each covariate column before computing the
  angles. Default `FALSE`, the statistic of Liu et al.

- tol:

  A difference vector shorter than `tol` times the largest covariate
  norm (at least 1) is treated as the zero vector.

## Value

An object of class `"htest"` with elements `statistic` (\\T\\),
`parameter` (`B`), `p.value`, `method`, `data.name`, and additionally
`boot` (the `B` bootstrap statistics), `n_failed` (refits that failed)
and `n_nonconverged` (refits that did not converge; their fitted values
are still used).

## Details

With \\e_i = y_i - \hat p_i\\ and \\X_i\\ the covariate vector of
observation \\i\\ (the model matrix without its intercept), the
statistic is \$\$T = n^{-2} \sum_i \sum_j \sum_l e_i e_j A\_{ijl},\$\$
\$\$A\_{ijl} = \int\_{S^p} I(X_i^T w \le X_l^T w) I(X_j^T w \le X_l^T
w)\\ dw,\$\$ with \\dw\\ the uniform probability measure on the unit
sphere. For \\u = X_i - X_l\\ and \\v = X_j - X_l\\ both nonzero,
\\A\_{ijl} = \\\pi - \angle(u, v)\\ / (2\pi)\\; it is \\1/2\\ when
exactly one of them is zero and \\1\\ when both are (Escanciano 2006,
Appendix). The constant does not affect the bootstrap p-value. Liu et
al. print the weight as a multiple of \\\arccos(\cdot)\\; the integral
fixes the orientation as \\\pi - \arccos(\cdot)\\, and a Monte Carlo
over random directions agrees with the form used here (see the package
tests).

The weight uses the covariates on their own scale, as Liu et al. do. The
angle is not invariant to rescaling one column, so `scale = TRUE`, which
standardizes each column first, gives a different test; it is offered
for covariates in unrelated units.

**Reference distribution.** The p-value comes from the model-based
bootstrap of Liu et al. (2024, Sec. 2.1; Dikta et al. 2006): `B`
responses are drawn from the fitted probabilities, the model is refitted
to each with the same design, link and offset, and \\p = (1 + \\\\T^\*
\ge T\\) / (B + 1)\\. A refit that fails scores \\T^\* = -\infty\\,
which can only make the test conservative; the count is returned. The
default `B = 1000` is the number of bootstrap samples Liu et al. use in
both of their data examples.

**Cost.** The weight matrix is computed once, in \\O(n^3 p)\\ time and
\\O(n^2)\\ memory, and each bootstrap replicate then costs one refit and
one quadratic form. In pure R the weight takes about 0.2 s at \\n =
200\\ and 3 to 5 s at \\n = 500\\, and the whole test with `B = 1000`
about 1 s and 6 s; at \\n = 2000\\ it takes minutes and several n-by-n
matrices of memory. With a single covariate the weight has an exact rank
form and is much faster.

## References

Escanciano JC (2006). "A consistent diagnostic test for regression
models using projections." *Econometric Theory*, 22(6), 1030–1051.
[doi:10.1017/S0266466606060506](https://doi.org/10.1017/S0266466606060506)

Liu H, Li X, Chen F, Haerdle W, Liang H (2024). "A comprehensive
comparison of goodness-of-fit tests for logistic regression models."
*Statistics and Computing*, 34, 175.
[doi:10.1007/s11222-024-10487-5](https://doi.org/10.1007/s11222-024-10487-5)

Dikta G, Kvesic M, Schmidt C (2006). "Bootstrap approximations in model
checks for binary data." *Journal of the American Statistical
Association*, 101(474), 521–530.
[doi:10.1198/016214505000001032](https://doi.org/10.1198/016214505000001032)

## See also

[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md),
where the test is the row `"Projection"`.

## Author

Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>

## Examples

``` r
set.seed(1)
n  <- 150
x1 <- rnorm(n); x2 <- rnorm(n)
y  <- rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2))
fit <- glm(y ~ x1 + x2, family = binomial())
projection.gof(fit, B = 99)
#> 
#>  Projection goodness-of-fit test (Escanciano 2006; Liu et al. 2024),
#>  model-based bootstrap
#> 
#> data:  y ~ x1 + x2
#> T = 0.012664, B = 99, p-value = 0.73
#> 

# \donttest{
## an omitted interaction: invisible to a test that looks only along the
## fitted linear predictor, visible along other directions
y2  <- rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2 + 1.5 * x1 * x2))
bad <- glm(y2 ~ x1 + x2, family = binomial())
projection.gof(bad)          # B = 1000
#> 
#>  Projection goodness-of-fit test (Escanciano 2006; Liu et al. 2024),
#>  model-based bootstrap
#> 
#> data:  y2 ~ x1 + x2
#> T = 0.027328, B = 1000, p-value = 0.06294
#> 
# }
```
