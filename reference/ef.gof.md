# Ebrahim-Farrington Goodness-of-Fit Test for Logistic Regression

Performs the Ebrahim-Farrington (EF) goodness-of-fit test for logistic
regression models: Farrington's (1996) correction applied to the
Hosmer-Lemeshow risk groups. It keeps the protection of the risk-ordered
table against gross recording errors and is the score test for records
of similar fitted risk sharing a miscalibration.

## Usage

``` r
ef.gof(
  y,
  predicted_probs = NULL,
  model = NULL,
  m = NULL,
  G = 10,
  method = NULL,
  reference = c("chisq", "normal"),
  X = NULL,
  groups = NULL
)
```

## Arguments

- y:

  A fitted binary logistic `glm` (then `predicted_probs` and `X` are
  taken from it), or a numeric vector of binary responses (0/1) for
  binary data / counts of successes for grouped data.

- predicted_probs:

  Numeric vector of predicted probabilities from the logistic regression
  model. Must be same length as `y`.

- model:

  Optional `glm` object. Required only for the original Farrington test
  with grouped data (when `m` is provided and `G` is NULL).

- m:

  Optional numeric vector of trial counts for each observation (for
  grouped data). If NULL, data is assumed to be binary.

- G:

  Number of risk groups: an integer (default 10); `"auto"` for groups of
  about 25 records, `max(10, ceiling(n / 25))`; or `"paul"` for the rule
  of Paul, Pennell and Lemeshow (2013), \\\max\\10, \min(n_1/2,
  (n-n_1)/2, 2 + 8(n/1000)^2)\\\\ rounded down, \\n_1\\ the number of
  events (designed for \\1000 \< n \le 25000\\). Ignored when `groups`
  is given. If NULL (and no `groups`), no grouping is performed and `m`
  must be provided.

- method:

  Deprecated. `"chisq"` is `reference = "chisq"`; `"normal"` is the
  standardization \\(EF - (G-2))/\sqrt{2(G-2)}\\ of package versions up
  to 2.9.0, kept, with a warning, to reproduce earlier results. It is
  not the normal reference of the EF paper, which is
  `reference = "normal"`.

- reference:

  `"chisq"` (default) or `"normal"`; see Details.

- X:

  Design matrix of the fitted model (with the intercept column), needed
  for `reference = "normal"` when `y` is not a fitted `glm`.

- groups:

  `NULL` (default) for equal-size groups by rank of the fitted risk,
  ties at random; `"pattern"` for one group per distinct fitted risk; or
  a vector of group labels, one per record.

## Value

For the grouped test, a one-row data frame with columns

- Test:

  `"Ebrahim-Farrington"`

- Test_Statistic:

  \\EF\\ for `reference = "chisq"`; the standardized \\z\\ for
  `reference = "normal"` (and for the deprecated `method = "normal"`)

- df:

  \\G - 2\\ for the chi-squared reference, `NA` otherwise

- Reference:

  `"chisq"`, `"normal"` or `"thesis"` (the deprecated
  `method = "normal"`)

- p_value:

  the p-value

- G:

  the number of groups used

- HL:

  the Hosmer-Lemeshow statistic on the same groups

- C:

  the removed directional part, \\C = HL - EF\\

and an attribute `"groups"`: one row per group with `group`, `n`,
`observed`, `expected`, `pbar`, `r` (\\r_g\\), `b` (\\b_g\\) and `c`
(\\c_g = b_g r_g\\, summing to `C`). For the original Farrington test
(`G = NULL`), a data frame with `Test`, `Test_Statistic` (standardized)
and `p_value`.

## Details

The records are sorted by fitted risk and cut into `G` groups of sizes
that differ by at most one. For group \\g\\ with \\n_g\\ records,
\\o_g\\ observed and \\e_g = n_g \bar\pi_g\\ expected events and \\V_g =
n_g\bar\pi_g(1-\bar\pi_g)\\, the standardized residual is \\r_g = (o_g -
e_g)/\sqrt{V_g}\\ and \\HL = \sum_g r_g^2\\. With \\b_g = (1 -
2\bar\pi_g)/\sqrt{V_g}\\, \$\$EF = HL - C, \qquad C = \sum_g b_g r_g
.\$\$ \\EF\\ keeps only the within-group products of the residuals: \\HL
= G + C + Q\\ and \\EF = G + Q\\ (Ebrahim, Khattab and El-Kotory, 2026).

**Reference distributions.** `reference = "chisq"` (the default) refers
\\EF\\ to \\\chi^2\_{G-2}\\, the reference of the Hosmer-Lemeshow test,
which the fixed-\\G\\ theory justifies. `reference = "normal"` is the
plug-in normal reference for many groups (Theorem 3 of the paper): \\z =
(EF - \hat\mu)/\hat\sigma\\ with model-based moments evaluated at the
fit, referred to the upper tail of the standard normal. It needs the
model's design matrix: pass the fitted `glm`, or `X`.

**Which reference, how many groups.**

- With `G = 10` and \\\chi^2\_{G-2}\\, each group needs about one
  expected event and one expected non-event. When some group expects
  fewer than one half, the test is conservative and the Hosmer-Lemeshow
  reference is unreliable too: use fewer groups. A warning of class
  `ef_sparse_groups` says so.

- With many groups, as under the rule of Paul, Pennell and
  Lemeshow (2013) (`G = "paul"`), use `reference = "normal"`, not
  \\\chi^2\_{G-2}\\, and prefer groups of about 25 records
  (`G = "auto"`, that is `max(10, ceiling(n / 25))`): they keep the
  protection against gross errors and do not dilute a smooth misfit over
  very many small groups.

**Ties.** Records with tied fitted risks must not be ordered by the
response: a group boundary inside a covariate pattern then manufactures
misfit. Ties are broken at random, so a result with tied fitted risks
depends on the random number stream; call
[`set.seed()`](https://rdrr.io/r/base/Random.html) first to reproduce
it. Nothing is drawn when no fitted risks are tied. Alternatively,
`groups = "pattern"` keeps each pattern (records with the same fitted
risk) in one group of its own, as for binomial data expanded to binary
records, or `groups` takes any vector of group labels.

**Directional read-out.** The returned `C` and the per-group terms \\c_g
= b_g r_g\\ in `attr(, "groups")` describe the direction of the misfit.
Since \\b_g\\ changes sign at \\\bar\pi_g = 1/2\\, \\C \< 0\\ says that,
on balance, the model underestimates how fast the risk approaches one in
the high-risk groups or approaches zero in the low-risk groups. The
read-out is descriptive.

For grouped binomial data (`G = NULL` with `m` and `model`), the
original Farrington test is applied with its full variance calculations.

## Note

- For grouped data (`m` provided *and* `G` set to `NULL`): Use the
  original Farrington test, which requires the fitted model object.
  Leaving `G` at its default keeps the automatic grouping, and `m` and
  `model` are then ignored.

- For binary data with `m=1` for all observations and no grouping, the
  test is not applicable and will return a p-value of 1.

## References

Ebrahim EK, Khattab IG, El-Kotory A (2026). "A Modified Hosmer-Lemeshow
Goodness-of-Fit Test for Asymmetric Links: Second-Order Power and
Robustness." arXiv:2607.15454 \[stat.ME\].
[doi:10.48550/arXiv.2607.15454](https://doi.org/10.48550/arXiv.2607.15454)

Farrington CP (1996). "On Assessing Goodness of Fit of Generalized
Linear Models to Sparse Data." *Journal of the Royal Statistical
Society, Series B*, **58**(2), 349-360.
[doi:10.1111/j.2517-6161.1996.tb02086.x](https://doi.org/10.1111/j.2517-6161.1996.tb02086.x)

Hosmer DW, Lemeshow S (1980). "A goodness-of-fit test for the multiple
logistic regression model." *Communications in Statistics - Theory and
Methods*, **9**(10), 1043-1069.
[doi:10.1080/03610928008827941](https://doi.org/10.1080/03610928008827941)

Paul P, Pennell ML, Lemeshow S (2013). "Standardizing the power of the
Hosmer-Lemeshow goodness of fit test in large data sets." *Statistics in
Medicine*, **32**(1), 67-80.
[doi:10.1002/sim.5525](https://doi.org/10.1002/sim.5525)

Ebrahim EK (2026). "Goodness-of-Fit Tests and Calibration
Machine-Learning Algorithms for Logistic Regression with Sparse Data."
M.Sc. thesis, Alexandria University. arXiv:2608.11140 \[stat.ME\].
[doi:10.48550/arXiv.2608.11140](https://doi.org/10.48550/arXiv.2608.11140)

## See also

[`hoslem.test`](https://rdrr.io/pkg/ResourceSelection/man/hoslem.test.html)
for the Hosmer-Lemeshow test

## Author

Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>

## Examples

``` r
# Example 1: Binary data with automatic grouping (Ebrahim-Farrington test)
set.seed(123)
n <- 500
x <- rnorm(n)
linpred <- 0.5 + 1.2 * x
prob <- 1 / (1 + exp(-linpred))
y <- rbinom(n, 1, prob)

# Fit logistic regression
model <- glm(y ~ x, family = binomial())

# Ten groups, chi-square(G - 2) reference
result <- ef.gof(model)
print(result)
#>                 Test Test_Statistic df Reference   p_value  G       HL
#> 1 Ebrahim-Farrington       2.997733  8     chisq 0.9344997 10 2.855296
#>            C
#> 1 -0.1424368
attr(result, "groups")              # r_g, b_g and the read-out c_g
#>    group  n observed  expected      pbar           r           b             c
#> 1      1 50        8  6.862475 0.1372495  0.46749624  0.29816388  0.1393904939
#> 2      2 50       10 13.349544 0.2669909 -1.07077398  0.14897555 -0.1595191431
#> 3      3 50       20 18.988563 0.3797713  0.29472498  0.07006745  0.0206506279
#> 4      4 50       25 22.796111 0.4559222  0.62579037  0.02503161  0.0156645437
#> 5      5 50       27 26.875711 0.5375142  0.03525373 -0.02128123 -0.0007502427
#> 6      6 50       31 30.799159 0.6159832  0.05839942 -0.06744976 -0.0039390273
#> 7      7 50       34 34.584493 0.6916899 -0.17899647 -0.11740714  0.0210154635
#> 8      8 50       38 38.624532 0.7724906 -0.21068021 -0.18384438  0.0387323716
#> 9      9 50       41 42.477714 0.8495543 -0.58454747 -0.27655029  0.1616567711
#> 10    10 50       48 46.641700 0.9328340  0.76742131 -0.48909071 -0.3753386298

# Groups of about 25 records with the normal reference
ef.gof(model, G = "auto", reference = "normal")
#>                 Test Test_Statistic df Reference   p_value  G       HL
#> 1 Ebrahim-Farrington      0.2564992 NA    normal 0.3987827 20 19.74056
#>            C
#> 1 0.09711975

# Frozen predictions: y and the predicted probabilities
ef.gof(y, fitted(model), G = 10)
#>                 Test Test_Statistic df Reference   p_value  G       HL
#> 1 Ebrahim-Farrington       2.997733  8     chisq 0.9344997 10 2.855296
#>            C
#> 1 -0.1424368

# Example 2: Grouped data (original Farrington test)
set.seed(456)
n_groups <- 50
m_trials <- sample(5:20, n_groups, replace = TRUE)
x_grouped <- rnorm(n_groups)
linpred_grouped <- -0.5 + 1.0 * x_grouped
prob_grouped <- 1 / (1 + exp(-linpred_grouped))
y_grouped <- rbinom(n_groups, m_trials, prob_grouped)

# Fit model for grouped data
data_grouped <- data.frame(successes = y_grouped, trials = m_trials, x = x_grouped)
model_grouped <- glm(cbind(successes, trials - successes) ~ x,
                     data = data_grouped, family = binomial())
predicted_probs_grouped <- fitted(model_grouped)

# Original Farrington test. G = NULL is required: left at its default of 10 the
# call takes the automatic-grouping branch instead, which ignores 'model' and 'm'
# and refers binomial counts to the binary statistic.
result_grouped <- ef.gof(y_grouped, predicted_probs_grouped,
                         model = model_grouped, m = m_trials,
                         G = NULL)
print(result_grouped)
#>                  Test Test_Statistic   p_value
#> 1 Farrington-Original     -0.2010413 0.5796669
```
