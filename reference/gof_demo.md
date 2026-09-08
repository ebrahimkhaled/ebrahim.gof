# Synthetic binary outcome data with a smooth calibration misfit

A small, fully synthetic dataset for demonstrating the goodness-of-fit
and calibration battery. It was generated reproducibly (see
`data-raw/make_gof_demo.R`) from a logistic data-generating process
whose true linear predictor includes a *quadratic* term in
(standardized) age. A model that regresses `outcome` on `age` linearly
(together with `bmi`, `sex` and `treatment`) is therefore mildly
misspecified, through a smooth, low-dimensional calibration distortion.
This is the regime in which the directed Ebrahim–Farrington / EDGE test
([`edge.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.gof.md),
[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md))
is designed to have more power than classical omnibus tests such as
Hosmer–Lemeshow.

## Usage

``` r
gof_demo
```

## Format

A data frame with 800 rows and 5 variables:

- outcome:

  binary response, 0/1 (event rate about 0.27).

- age:

  continuous covariate, years (range about 20–70). The true model
  depends on age quadratically.

- bmi:

  continuous covariate, body mass index in kg/m^2.

- sex:

  binary covariate, 0 = female, 1 = male.

- treatment:

  binary covariate, 0 = control, 1 = treated.

## Source

Simulated; see `data-raw/make_gof_demo.R` in the package sources.

## Details

The true data-generating linear predictor is \$\$\eta = -0.6 + 0.8 z_a -
0.7 z_a^2 + 0.5 z_b + 0.4\\\mathrm{sex} - 0.3\\\mathrm{treatment},\$\$
where \\z_a = (\mathrm{age} - 45)/14\\ and \\z_b = (\mathrm{bmi} -
27)/4\\, and \\\Pr(\mathrm{outcome} = 1) = \mathrm{plogis}(\eta)\\.

## Examples

``` r
data("gof_demo", package = "ebrahim.gof")
fit <- glm(outcome ~ age + bmi + sex + treatment,
           data = gof_demo, family = binomial)
edge.gof(fit)
#>   Test Basis Test_Statistic       df        Method    p_value
#> 1 EDGE poly3       8.329512 2.072814 satterthwaite 0.01496412
```
