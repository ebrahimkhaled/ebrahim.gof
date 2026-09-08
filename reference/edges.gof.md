# EDGES: Cauchy-Combination Ensemble of Directed GOF Tests (Alias)

Alias for
[`def.ensemble.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md);
see the EDGES paper. `edges.gof()` is the brand name (EDGES = the
Cauchy-combination ensemble of the EDGE directed bases) used in the
manuscript. It takes exactly the same arguments as
[`def.ensemble.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md)
and returns exactly the same value; the legacy name
[`def.ensemble.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md)
is retained unchanged for back-compatibility.

## Usage

``` r
edges.gof(...)
```

## Arguments

- ...:

  Arguments passed on to
  [`def.ensemble.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md)
  (e.g. `object`, `predicted_probs`, `X`, `components`, `add_ef`,
  `combine`, `G`, `extra_pvalues`).

## Value

A one-row `data.frame` with columns `Test`, `Combiner`, `Components`,
`k`, and `p_value`.

## References

Liu, Y. and Xie, J. (2020). Cauchy combination test. *JASA*, 115(529),
393-402.

## See also

[`def.ensemble.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md)
(legacy name),
[`edge.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.gof.md),
[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md).

## Author

Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>

## Examples

``` r
set.seed(1)
x <- runif(500, -3, 3)
y <- rbinom(500, 1, plogis(0.6 * x))
fit <- glm(y ~ x, family = binomial())
edges.gof(fit)                 # identical to def.ensemble.gof(fit)
#>           Test Combiner         Components k   p_value
#> 1 DEF ensemble      cct poly2+poly3+stukel 3 0.7654395
```
