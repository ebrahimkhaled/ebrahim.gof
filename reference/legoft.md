# LEGofT: frozen-weight combination goodness-of-fit test for binary logistic regression

Combines eleven classical and directed goodness-of-fit statistics with
weights that were fixed offline and ship frozen, and calibrates the
combination by a parametric bootstrap at the fitted parameters. Nothing
is retrained when you call it.

## Usage

``` r
legoft(
  object,
  B = 199,
  seed = NULL,
  weights = NULL,
  panel = c("legoft", "equal", "full", "classical", "one.per.group", "protected"),
  rule = c("cauchy", "hmp", "fisher", "minp")
)
```

## Arguments

- object:

  a fitted binomial `glm`.

- B:

  number of parametric-bootstrap replicates for the reference
  distribution. 199 gives a smallest attainable p-value of 0.005; raise
  it for smaller p-values.

- seed:

  optional integer for reproducibility.

- weights:

  optional named vector of member weights; defaults to the frozen rule.
  Supplying your own makes the result no longer the shipped procedure –
  say so if you report it. Members may be any of the eleven frozen-rule
  tests plus `"Stukel"`, `"Osius-Rojek"` and `"le-Cessie"`.

- panel:

  which tests to combine. `"legoft"` (default) is the frozen rule. The
  others combine a panel of the panel paper with equal weights:
  `"equal"` (the eleven frozen-rule tests), `"full"` (all fourteen
  tests), `"classical"` (Hosmer-Lemeshow, Osius-Rojek, le Cessie-van
  Houwelingen, Stukel, Tsiatis, Xie), `"one.per.group"` (one test from
  each of the five groups of correlated tests) and `"protected"`
  (Hosmer-Lemeshow, EF, Pigeon-Heyse, le Cessie-van Houwelingen; robust
  to a few recording errors in a covariate). Ignored when `weights` is
  given.

- rule:

  combination statistic: `"cauchy"` (default), `"hmp"` (harmonic mean
  p-value), `"fisher"` or `"minp"` (minimum p-value). The last three use
  equal weights, so they need a panel other than `"legoft"` or equal
  `weights`. Every rule is calibrated by the same parametric bootstrap.

## Value

an object of class `"legoft"`: the statistic, its bootstrap p-value, the
member p-values, the weights, the panel, the rule and `B_used`.

## Details

The p-value is exact in finite samples when the null parameters are
known. With the parameters estimated – the case here – the calibration
is asymptotic; simulation at \\n = 500\\ put the empirical size at 0.037
against a nominal 0.05.

The panel paper recommends `panel = "equal"` or `"full"` with the
default Cauchy rule for a single decision from many tests, and
`panel = "protected"` when a few covariate values may be recorded
wrongly. In its simulations the bootstrap-calibrated combination kept
the 5% level in three null designs. Panels that contain le Cessie-van
Houwelingen are slower, because that test builds an \\n \times n\\
matrix for the data and for every bootstrap replicate.

## See also

[`legoft.localize`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/legoft.localize.md)
for which domain of evidence carries the misfit.

## Examples

``` r
# \donttest{
set.seed(1)
x1 <- runif(300, -3, 3); x2 <- rnorm(300)
y  <- rbinom(300, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.25 * x1^2))
fit <- glm(y ~ x1 + x2, family = binomial())
legoft(fit, B = 99, seed = 1)                       # frozen LEGofT weights
#> 
#> LEGofT (frozen weights, parametric-bootstrap calibration)
#> 
#>   statistic = 9520.2744    p-value = 0.0200    (B = 99)
#> 
legoft(fit, B = 99, seed = 1, panel = "equal")      # equal weights, same tests
#> 
#> Panel combination (panel = equal, rule = cauchy, 11 tests, parametric-bootstrap calibration)
#> 
#>   statistic = 3840.2194    p-value = 0.0200    (B = 99)
#> 
legoft(fit, B = 99, seed = 1, panel = "protected")  # robust to recording errors
#> 
#> Panel combination (panel = protected, rule = cauchy, 4 tests, parametric-bootstrap calibration)
#> 
#>   statistic = 17.6830    p-value = 0.0300    (B = 99)
#> 
# }
```
