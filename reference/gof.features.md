# Goodness-of-fit evidence features for a fitted model

Builds the evidence vector used by the learned-ensemble goodness-of-fit
test: one-sided z-scores \\\Phi^{-1}(1-p)\\ from a panel of GOF tests
plus the covariate-space directed tests. Larger values mean stronger
evidence of misfit.

## Usage

``` r
gof.features(
  object,
  tests = c("HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF", "DEF.poly2",
    "DEF.poly3", "DEF.stukel")
)
```

## Arguments

- object:

  A fitted binary logistic [`glm`](https://rdrr.io/r/stats/glm.html).

- tests:

  Character vector of
  [`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
  test names to use as panel features (default: a fast partition +
  DEF-family panel).

## Value

A named numeric vector of evidence features.

## See also

[`deploy.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deploy.gof.md),
[`cdef.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/cdef.gof.md),
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md).
