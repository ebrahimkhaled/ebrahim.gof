# Run the External-Validation Tests at Once

Checks the calibration of predictions that were made *without* the data
at hand – a published model, or any model, applied to new patients – and
returns one tidy `data.frame`, one row per test, in the format of
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md).
Only the outcomes `y` and the predicted probabilities `p` are needed, so
the predictions may come from a logistic regression, a random forest, a
neural network or a clinical score: nothing is refitted.

## Usage

``` r
run.all.external(y, p, G = 10, X = NULL, include_slow = FALSE)
```

## Arguments

- y:

  Binary (0/1) outcomes of the validation sample.

- p:

  Predicted probabilities for the same records, produced without them.

- G:

  Number of risk groups for the directed and Hosmer–Lemeshow tests
  (default `10`).

- X:

  Optional covariate matrix or data frame, for le Cessie's test only.

- include_slow:

  Logical; run le Cessie's test when `X` is given (default `FALSE`).

## Value

A `data.frame` of class `gof_battery` with columns `Test`, `Family`,
`Statistic`, `df`, `p_value` and `Note`, printed by the same method as
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md).

## Details

**Why a separate battery.** `run.all.gof(y, predicted_probs = p)` treats
the predictions as fitted to `y`, and each reference distribution then
allows for the parameters the fit spent. When the predictions are frozen
nothing was spent: the Hosmer–Lemeshow statistic is referred to
\\\chi^2_G\\, not \\\chi^2\_{G-2}\\; Stukel's terms are added to the
frozen linear predictor as an offset; and the directed test takes
\\\Omega = I\\ with a constant column, because no score equation absorbs
the overall level. Using the internal references on frozen predictions
makes every one of these tests conservative.

**The tests** (`Family` in brackets).

- `EDGE (G=...)` \[Directed\] – the directed grouped test in external
  mode, `def.gof(y, predicted_probs = p, external = TRUE)`: the cubic
  basis plus a constant, four degrees of freedom. Run at `G` (ten by
  default, the setting that keeps its level when a few records carry
  corrupted predictors) and at `G = "auto"`, \\\max(10, \lceil n/25
  \rceil)\\, which has more power on clean data.

- `Cox recalibration` \[Calibration\] – the likelihood-ratio test that
  the calibration intercept is 0 and the slope 1 in `glm(y ~ logit(p))`
  (Cox 1958; Miller et al. 1991). The Note gives both estimates.

- `Calibration-in-the-large` \[Calibration\] – the score test that the
  intercept is 0 with the slope held at 1: \\(O - E)/\sqrt{\sum
  p(1-p)}\\.

- `Spiegelhalter-z` \[Calibration\] – Spiegelhalter's (1986) \\z\\ test,
  two-sided.

- `GiViTI` \[Calibration\] – the GiViTI calibration test in its external
  mode (Nattino et al. 2014), run in an isolated callr process; needs
  givitiR and callr. When it selects a polynomial of degree one it
  returns the same \\p\\-value as the Cox test.

- `Hosmer-Lemeshow (external)` \[Partition\] – \\G\\ equal-frequency
  risk groups, referred to \\\chi^2_G\\.

- `Stukel (offset)` \[Directed\] – Stukel's two terms added to the
  frozen linear predictor, likelihood-ratio \\\chi^2_2\\.

- `le Cessie (external)` \[Smoothing\] – le Cessie and van Houwelingen's
  kernel statistic over covariate space with \\\Omega = I\\; only when
  `X` is given and `include_slow = TRUE`, because it builds an \\n
  \times n\\ kernel.

Three descriptive rows carry no \\p\\-value: the ratio of observed to
expected events, the calibration slope, and the c-statistic (area under
the ROC curve).

**Reading the panel.** In a simulation of external validation (EDGE
paper, Supporting Information) the directed test at ten groups had the
power of the Cox test and the GiViTI belt on average, led them on curved
departures and trailed them on a shift, and kept its level with ten
reversed predictions in 1000 records, where the Cox test and the belt
did not. The Cox test says whether the predictions need recalibrating;
the directed test says whether a recalibration would be enough.

## See also

[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md)
(its `external` argument),
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md).

## Author

Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>

## Examples

``` r
set.seed(1)
n <- 1000
x <- rnorm(n)
p <- plogis(-1 + 0.8 * x)                 # a published model, frozen
y <- rbinom(n, 1, plogis(-1 + 0.6 * x))    # new patients: the model is overfitted
run.all.external(y, p)
#> 
#> Goodness-of-fit battery: 11 tests  (5 reject at 0.05)
#> ============================================================
#>  Test                            Statistic   df  p-value    
#>  --- Partition ---------------------------------------------
#>  Hosmer-Lemeshow (external) [a]      24.66   10   0.0060 ** 
#>  --- Directed ----------------------------------------------
#>  EDGE (G=10) [b]                     14.86    4   0.0050 ** 
#>  EDGE (G=auto, 40) [b]               17.42    4   0.0016 ** 
#>  Stukel (offset) [c]                 5.342    2   0.0692 .  
#>  --- Calibration -------------------------------------------
#>  Cox recalibration [d]                13.2    2   0.0014 ** 
#>  Calibration-in-the-large [e]       0.8792    1   0.3484    
#>  Spiegelhalter-z [f]                 1.565        0.1176    
#>  GiViTI [g]                                       0.0014 ** 
#>  --- Descriptive -------------------------------------------
#>  O/E ratio [h]                      0.9566             -    
#>  Calibration slope [i]              0.6643             -    
#>  c-statistic (AUC) [j]              0.6474             -    
#> ------------------------------------------------------------
#>  Signif.:  *** <.001   ** <.01   * <.05   . <.1
#>  Notes:
#>    [a] 10 groups, chi-square(10)
#>    [b] external mode, cubic basis + constant
#>    [c] terms added to logit(p)
#>    [d] intercept -0.337, slope 0.664
#>    [e] O/E = 0.957
#>    [f] two-sided z
#>    [g] calibration belt; devel=external
#>    [h] 1 = calibrated in the large
#>    [i] 1 = calibrated; < 1 = overfitted
#>    [j] discrimination, not calibration
```
