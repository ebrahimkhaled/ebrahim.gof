# Grouped-covariate companion to `gof_demo` (replicated covariate patterns)

A companion dataset to
[`gof_demo`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof_demo.md)
built from the *same* data-generating process and seed discipline (see
`data-raw/make_gof_demo.R`), except that the covariates are coarsened
*before* the linear predictor is computed: `age` is rounded to 10-year
bins (20, 30, ..., 70) and `bmi` to whole integers. The recorded
covariates are therefore exactly the covariates the outcome was
generated from, and many observations share a covariate pattern (328
distinct patterns among 800 observations, versus one pattern per
observation in `gof_demo`).

## Usage

``` r
gof_demo_grouped
```

## Format

A data frame with 800 rows and 5 variables:

- outcome:

  binary response, 0/1 (event rate about 0.27).

- age:

  age in years, rounded to 10-year bins (20, 30, ..., 70). The true
  model depends on (binned) age quadratically.

- bmi:

  body mass index in kg/m^2, rounded to whole integers.

- sex:

  binary covariate, 0 = female, 1 = male.

- treatment:

  binary covariate, 0 = control, 1 = treated.

## Source

Simulated; see `data-raw/make_gof_demo.R` in the package sources.

## Details

Its purpose is to demonstrate the sparse-versus-grouped distinction that
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
surfaces: the battery reports the per-observation ("sparse") and
per-covariate-pattern ("grouped") forms of the pattern-sensitive tests
(Pearson, deviance, McCullagh) side by side, and on replicated-pattern
data such as this the two forms can disagree on the same fitted model.
On the all-continuous `gof_demo` (every observation its own pattern) the
two forms coincide – the degenerate case.

The true data-generating linear predictor has the same form as for
[`gof_demo`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof_demo.md),
\$\$\eta = -0.6 + 0.8 z_a - 0.7 z_a^2 + 0.5 z_b + 0.4\\\mathrm{sex} -
0.3\\\mathrm{treatment},\$\$ with \\z_a = (\mathrm{age} - 45)/14\\ and
\\z_b = (\mathrm{bmi} - 27)/4\\ computed from the *binned* `age` and
`bmi`, and \\\Pr(\mathrm{outcome} = 1) = \mathrm{plogis}(\eta)\\.

## See also

[`gof_demo`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof_demo.md),
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)

## Examples

``` r
data("gof_demo_grouped", package = "ebrahim.gof")
fit <- glm(outcome ~ age + bmi + sex + treatment,
           data = gof_demo_grouped, family = binomial)
# sparse and grouped forms reported side by side:
run.all.gof(fit, include_slow = FALSE, install = "no")
#> 
#> Goodness-of-fit battery: 24 tests  (7 reject at 0.05)
#> ========================================================
#>  Test                        Statistic   df  p-value    
#>  --- Global --------------------------------------------
#>  Pearson                         769.2  795   0.7377    
#>  Pearson (grouped) [a]           346.8  323   0.1738    
#>  Deviance                        833.7  795   0.1655    
#>  Deviance (grouped) [a]          369.5  323   0.0380 *  
#>  Information-Matrix              20.94    5   0.0008 ***
#>  --- Standardized --------------------------------------
#>  Osius-Rojek [a]                0.9196        0.3578    
#>  McCullagh                      -1.864        0.9689    
#>  McCullagh (grouped) [a]         0.847        0.1985    
#>  Copas-RSS                       2.189    1   0.0286 *  
#>  EF                              1.038    8   0.1445    
#>  EF-normal [b]                   1.038    8   0.1496    
#>  --- Partition -----------------------------------------
#>  HL                              11.79    8   0.1610    
#>  HL-equalwidth                   10.33    7   0.1706    
#>  Pigeon-Heyse                    11.57    9   0.2387    
#>  F-test [c]                      1.535    9   0.1312    
#>  --- Covariate-space -----------------------------------
#>  Tsiatis                         52.87    8  1.1e-08 ***
#>  Xie                             48.98    7  2.3e-08 ***
#>  Pulkstenis-Robinson [d]         84.93   43   0.0001 ***
#>  --- Directed ------------------------------------------
#>  DEF.poly2                       3.406 1.08   0.0645 .  
#>  DEF.poly3                       4.657 2.08   0.0974 .  
#>  DEF.stukel                        3.3 1.61   0.0936 .  
#>  Stukel                          6.057    2   0.0484 *  
#>  --- Ensemble ------------------------------------------
#>  Ensemble.Vote(3DEF) [e]                      0.0824 .  
#>  Ensemble.Univ(3DEF+EF) [e]                   0.0924 .  
#> --------------------------------------------------------
#>  Signif.:  *** <.001   ** <.01   * <.05   . <.1
#>  Notes:
#>    [a] grouped to 328 covariate patterns
#>    [b] normal reference (thesis)
#>    [c] deviance residuals ~ groups (ANOVA F)
#>    [d] split on categorical: age, sex, treatment
#>    [e] Cauchy combination of the directed tests
```
