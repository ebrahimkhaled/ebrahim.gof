# Print a goodness-of-fit battery

Formats the
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
result as a compact, readable table: rows grouped by test family,
p-values shown to four decimals (or scientific for very small values,
`"-"` when not available), and a significance flag. The object is still
a plain `data.frame` underneath, so all the raw columns remain available
for programmatic use.

## Usage

``` r
# S3 method for class 'gof_battery'
print(x, ...)
```

## Arguments

- x:

  A `gof_battery` object returned by
  [`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md).

- ...:

  Ignored.

## Value

`x`, invisibly.

## Details

The header counts the rows that reject at the 0.05 level. The rows
printed for comparison only (the sparse `Pearson` and `Deviance`
statistics and the `F-test`; see
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md))
carry no significance flag and are not counted. Notes are wrapped to
`getOption("width")`.

## Examples

``` r
fit <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial)
print(run.all.gof(fit, include_slow = FALSE, install = "no"))
#> 
#> Goodness-of-fit battery: 25 tests  (18 reject at 0.05)
#> ========================================================
#>  Test                        Statistic   df  p-value    
#>  --- Global --------------------------------------------
#>  Pearson [a]                       760  795   0.8087    
#>  Deviance [a]                    839.1  795   0.1352    
#>  Information-Matrix              26.39    5  7.5e-05 ***
#>  --- Standardized --------------------------------------
#>  Osius-Rojek [b]                -2.388        0.0169 *  
#>  McCullagh [c]                  -2.725        0.9968    
#>  Copas-RSS [b]                   3.352        0.0008 ***
#>  --- Partition -----------------------------------------
#>  HL                              15.52    8   0.0499 *  
#>  HL-largeN [d]                   15.52    8   0.0501 .  
#>  HL-equalwidth                   22.11    7   0.0024 ** 
#>  Pigeon-Heyse                    15.58    9   0.0761 .  
#>  F-test [e]                      2.088    9   0.0283    
#>  EF                              15.95    8   0.0431 *  
#>  EF-normal [f]                   1.767        0.0386 *  
#>  --- Covariate-space -----------------------------------
#>  Tsiatis                         63.53    9  2.8e-10 ***
#>  Xie                              51.9    7  6.1e-09 ***
#>  Pulkstenis-Robinson [g]         7.282    4   0.1217    
#>  --- Directed ------------------------------------------
#>  EDGE [h]                        11.21 2.01   0.0035 ** 
#>  EDGE.G10 [i]                     8.33 2.07   0.0150 *  
#>  DEF.poly2 [i]                   8.299 1.08   0.0036 ** 
#>  DEF.stukel [i]                    6.5 1.56   0.0118 *  
#>  DEF.sym [i]                     1.519    1   0.0096 ** 
#>  Stukel                          11.58    2   0.0031 ** 
#>  --- Calibration ---------------------------------------
#>  Cubic-LR [j]                    12.44    2   0.0020 ** 
#>  --- Ensemble ------------------------------------------
#>  Ensemble.Vote(3DEF) [k]                      0.0070 ** 
#>  Ensemble.Univ(3DEF+EF) [k]                   0.0088 ** 
#> --------------------------------------------------------
#>  Signif.:  *** <.001   ** <.01   * <.05   . <.1
#>  Notes:
#>    [a] sparse form; comparison only, not counted
#>    [b] normal deviate, two-sided
#>    [c] normal deviate, upper tail
#>    [d] H0: eps <= 2.74e-03 (n0 = 1e+06); eps_hat = 9.69e-02
#>    [e] ANOVA F on deviance residuals; comparison only, not counted
#>    [f] normal reference, many groups (G = 32)
#>    [g] split on categorical: sex, treatment
#>    [h] S/c ~ chi-square(df), Satterthwaite; G = 32 (auto)
#>    [i] S/c ~ chi-square(df), Satterthwaite
#>    [j] eta^2 and eta^3 added to a logit fit on eta
#>    [k] Cauchy combination of the directed tests
```
