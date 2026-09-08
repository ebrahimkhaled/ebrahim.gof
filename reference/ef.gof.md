# Ebrahim-Farrington Goodness-of-Fit Test for Logistic Regression

Performs the Ebrahim-Farrington goodness-of-fit test for logistic
regression models. This test is particularly effective for binary data
and sparse datasets, providing an improved alternative to the
traditional Hosmer-Lemeshow test.

## Usage

``` r
ef.gof(
  y,
  predicted_probs = NULL,
  model = NULL,
  m = NULL,
  G = 10,
  method = c("chisq", "normal")
)
```

## Arguments

- y:

  A fitted binary logistic `glm` (then `predicted_probs` is taken from
  it automatically), or a numeric vector of binary responses (0/1) for
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

  Optional integer specifying the number of groups for binary data
  grouping. Default is 10. If NULL, no grouping is performed and `m`
  must be provided.

- method:

  Reference distribution for the grouped EF statistic: `"chisq"`
  (default) refers \\T\_{EF}\\ to a \\\chi^2\_{G-2}\\ distribution;
  `"normal"` uses the standardized \\Z\_{EF}\\ (the behaviour of package
  versions \<= 1.0.0).

## Value

A data frame with the following columns:

- Test:

  Character string identifying the test performed

- Test_Statistic:

  Numeric value of the standardized test statistic

- p_value:

  Numeric p-value for the test

## Details

The Ebrahim-Farrington test is based on Farrington's (1996) theoretical
framework but simplified for practical implementation with binary data.
The test uses a modified Pearson chi-square statistic with
data-dependent grouping, where observations are grouped by their
predicted probabilities.

For binary data (when `G` is specified), the test automatically groups
observations into `G` groups based on predicted probabilities and
applies the simplified Ebrahim-Farrington statistic:

\$\$Z\_{EF} = \frac{T\_{EF} - (G - 2)}{\sqrt{2(G-2)}}\$\$

where \\T\_{EF}\\ is the modified Pearson chi-square statistic, and
\\G\\ is the number of groups.

For grouped data (when `m` is provided), the test applies the original
Farrington test with full variance calculations.

## Note

- For binary data with automatic grouping (`G` specified): Use the
  Ebrahim-Farrington test which is computationally efficient and doesn't
  require the model specification.

- For grouped data (`m` provided *and* `G` set to `NULL`): Use the
  original Farrington test, which requires the fitted model object.
  Leaving `G` at its default keeps the automatic grouping, and `m` and
  `model` are then ignored.

- The test statistic follows a standard normal distribution under the
  null hypothesis of adequate model fit.

- For binary data with `m=1` for all observations and no grouping, the
  test is not applicable and will return a p-value of 1.

## References

Farrington CP (1996). "On Assessing Goodness of Fit of Generalized
Linear Models to Sparse Data." *Journal of the Royal Statistical
Society, Series B*, **58**(2), 349-360.
[doi:10.1111/j.2517-6161.1996.tb02086.x](https://doi.org/10.1111/j.2517-6161.1996.tb02086.x)

Hosmer DW, Lemeshow S (1980). "A goodness-of-fit test for the multiple
logistic regression model." *Communications in Statistics - Theory and
Methods*, **9**(10), 1043-1069.
[doi:10.1080/03610928008827941](https://doi.org/10.1080/03610928008827941)

Ebrahim EK, El-Kotory A (2026). "A Directional Hosmer-Lemeshow
Goodness-of-Fit Test for Sparse Logistic Regression." arXiv:2607.15454
\[stat.ME\].
[doi:10.48550/arXiv.2607.15454](https://doi.org/10.48550/arXiv.2607.15454)

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
predicted_probs <- fitted(model)

# Perform Ebrahim-Farrington test with 10 groups
result <- ef.gof(y, predicted_probs, G = 10)
print(result)
#>                 Test Test_Statistic   p_value
#> 1 Ebrahim-Farrington      -1.250567 0.9344997

# Example 2: Compare with different number of groups
result_4 <- ef.gof(y, predicted_probs, G = 4)
result_20 <- ef.gof(y, predicted_probs, G = 20)

# Example 3: Grouped data (original Farrington test)
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
