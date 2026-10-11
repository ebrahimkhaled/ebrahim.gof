# Summarise a goodness-of-fit battery by family

Condenses a
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
or
[`run.all.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.external.md)
result to one line per test family: how many tests ran, how many
returned a p-value, how many reject at `alpha`, and the smallest p-value
with the test that gave it. It answers "which kind of departure do the
data point to" before the full table is read.

## Usage

``` r
# S3 method for class 'gof_battery'
summary(object, alpha = 0.05, ...)

# S3 method for class 'summary.gof_battery'
print(x, ...)
```

## Arguments

- object:

  A `gof_battery` object.

- alpha:

  Level at which a p-value counts as a rejection. Default `0.05`.

- ...:

  Not used.

- x:

  A `summary.gof_battery` object.

## Value

An object of class `"summary.gof_battery"`: a data frame with columns
`Family`, `Tests` (rows in the battery), `With_p` (rows with a p-value),
`Rejected` (p-values at most `alpha`, comparison rows left out),
`Smallest_p` and `Smallest_test`, in the family order of the printed
battery, with the attribute `alpha`.

## Details

The counts are descriptive. The tests within a family are strongly
correlated and no multiplicity adjustment is made, so a count of
rejections is not itself a test; for an error-controlled statement of
which part of the model is wrong, use
[`localize.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.gof.md)
or
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md).

The rows that
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
prints for comparison only (the sparse `Pearson` and `Deviance`
statistics and the `F-test`, which do not hold their level on ungrouped
binary data) are counted in `Tests` and `With_p` but not in `Rejected`,
and never give `Smallest_p`.

## See also

[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md),
[`print.gof_battery`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/print.gof_battery.md).

## Examples

``` r
fit <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial)
summary(run.all.gof(fit, include_slow = FALSE, install = "no"))
#> 
#> Goodness-of-fit battery by family: 25 tests, 25 with a p-value, 18 reject at
#> 0.05
#> 
#>  Family          Tests With p Reject Min p   Given by           
#>  Global          3     3      1      7.5e-05 Information-Matrix 
#>  Standardized    3     3      2      0.0008  Copas-RSS          
#>  Partition       7     7      4      0.0024  HL-equalwidth      
#>  Covariate-space 3     3      2      2.8e-10 Tsiatis            
#>  Directed        6     6      6      0.0031  Stukel             
#>  Calibration     1     1      1      0.0020  Cubic-LR           
#>  Ensemble        2     2      2      0.0070  Ensemble.Vote(3DEF)
#> 
#> Counts are descriptive: the tests in a family are correlated and not adjusted
#> for multiplicity. Rows printed for comparison only are not counted as
#> rejections.
```
