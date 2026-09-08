# Combine Directed GOF Tests into One Decision (Ensemble)

Combines the three Directed Ebrahim-Farrington (DEF) basis tests
(`"poly2"`, `"poly3"`, `"stukel"`) into a single goodness-of-fit
decision, so the user does not have to choose a basis. By default the
p-values are combined with the Cauchy Combination Test (CCT), which
controls the error rate under the strong dependence between tests
computed on the same fitted model. The omnibus EF test can optionally be
added to the vote.

## Usage

``` r
def.ensemble.gof(
  object,
  predicted_probs = NULL,
  X = NULL,
  components = c("poly2", "poly3", "stukel"),
  add_ef = FALSE,
  combine = c("cct", "minp", "fisher"),
  G = 10,
  extra_pvalues = NULL
)
```

## Arguments

- object:

  A fitted binary logistic [`glm`](https://rdrr.io/r/stats/glm.html), or
  a binary (0/1) vector `y` (then supply `predicted_probs`).

- predicted_probs:

  Numeric predicted probabilities; required when `object` is a `y`
  vector.

- X:

  Optional design matrix, threaded to
  [`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md)
  for the exact calibration (only used with the `y`/`predicted_probs`
  form).

- components:

  Character vector, a subset of `c("poly2","poly3","stukel")`. Default
  is all three.

- add_ef:

  Logical; if `TRUE`, the omnibus EF p-value
  ([`ef.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md))
  is appended to the components. Default `FALSE`.

- combine:

  One of `"cct"` (default), `"minp"`, `"fisher"`.

- G:

  Integer number of groups passed to `def.gof`/`ef.gof` (default 10).

- extra_pvalues:

  Optional named numeric vector of additional p-values to include (e.g.
  a Tsiatis test computed elsewhere). Default `NULL`.

## Value

A one-row `data.frame` with columns `Test`, `Combiner`, `Components`,
`k`, and `p_value`.

## Details

Because the component tests are computed on the same fit, their p-values
are strongly dependent. The CCT (`combine = "cct"`) has an asymptotic
standard-Cauchy null whose tail is robust to this dependence, so it
needs no calibration. The `"minp"` (Sidak) and `"fisher"` rules assume
independence and are offered for comparison only; under positive
dependence `"minp"` is conservative and `"fisher"` is anti-conservative,
so they should be calibrated by simulation before use (not done here).

## References

Liu, Y. and Xie, J. (2020). Cauchy combination test. *JASA*, 115(529),
393-402.

## See also

[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md),
[`ef.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md).

## Author

Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>

## Examples

``` r
data("gof_demo", package = "ebrahim.gof")
wrong <- glm(outcome ~ age + bmi + sex + treatment,
             data = gof_demo, family = binomial())
def.ensemble.gof(wrong)                 # CCT of the three DEF bases
#>           Test Combiner         Components k     p_value
#> 1 DEF ensemble      cct poly2+poly3+stukel 3 0.006959379
def.ensemble.gof(wrong, add_ef = TRUE)  # add the omnibus EF
#>           Test Combiner            Components k    p_value
#> 1 DEF ensemble      cct poly2+poly3+stukel+EF 4 0.00880698

## the corrected model, for contrast
right <- glm(outcome ~ poly(age, 2) + bmi + sex + treatment,
             data = gof_demo, family = binomial())
def.ensemble.gof(right)
#>           Test Combiner         Components k   p_value
#> 1 DEF ensemble      cct poly2+poly3+stukel 3 0.8090437
```
