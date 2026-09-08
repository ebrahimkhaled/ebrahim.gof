# Deployable learned-ensemble GOF test via parametric bootstrap

Turns a pre-trained ensemble `meta` into a deployable goodness-of-fit
test for *any* fitted model: it scores the model, then calibrates the
p-value by a per-dataset parametric bootstrap from the fitted model (so
no knowledge of the truth or the data-generating design is required).
Validity comes from the bootstrap, independent of how `meta` was
trained.

## Usage

``` r
deploy.gof(object, meta, B = 99, feature_fn = gof.features)
```

## Arguments

- object:

  A fitted binary logistic [`glm`](https://rdrr.io/r/stats/glm.html).

- meta:

  A pre-trained scorer: either a function `f(features)` returning a
  scalar misfit score, or an object with a `predict` method consuming a
  one-row feature `matrix`.

- B:

  Number of parametric-bootstrap resamples (default 99).

- feature_fn:

  Function mapping a fitted glm to its feature vector (default
  [`gof.features`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof.features.md)).

## Value

A one-row `data.frame` with the score, `B`, and the bootstrap `p_value`.

## See also

[`gof.features`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof.features.md),
[`cdef.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/cdef.gof.md).
