# Introduction to the Ebrahim-Farrington Goodness-of-Fit Test

## Introduction

The **ebrahim.gof** package implements the Ebrahim-Farrington
goodness-of-fit test for logistic regression models. This test is
particularly effective for binary data and sparse datasets, providing an
improved alternative to the traditional Hosmer-Lemeshow test.

## Background and Motivation

Goodness-of-fit testing is crucial in logistic regression to assess
whether the fitted model adequately describes the data. The most
commonly used test is the Hosmer-Lemeshow test, but it has several
limitations:

1.  **Limited power** for detecting certain types of model
    misspecification
2.  **Dependency on grouping strategy** which can affect results
3.  **Poor performance** with sparse data or continuous covariates

The Ebrahim-Farrington test addresses these limitations by using a
modified Pearson chi-square statistic based on Farrington’s (1996)
theoretical framework, but simplified for practical implementation with
binary data.

## Installation and Loading

``` r

# Install from GitHub
devtools::install_github("ebrahimkhaled/ebrahim.gof")

# Load the package
library(ebrahim.gof)
```

``` r

library(ebrahim.gof)
```

## Basic Usage

The main function
[`ef.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md)
performs the goodness-of-fit test:

``` r

# Simulate binary data
set.seed(123)
n <- 500
x <- rnorm(n)
linpred <- 0.5 + 1.2 * x
prob <- plogis(linpred)  # Convert to probabilities
y <- rbinom(n, 1, prob)

# Fit logistic regression
model <- glm(y ~ x, family = binomial())
predicted_probs <- fitted(model)

# Perform Ebrahim-Farrington test
result <- ef.gof(y, predicted_probs, G = 10)
print(result)
#>                 Test Test_Statistic   p_value
#> 1 Ebrahim-Farrington      -1.250567 0.9344997
```

## Understanding the Test Statistic

For binary data with automatic grouping, the Ebrahim-Farrington test
statistic is:

``` math
Z_{EF} = \frac{T_{EF} - (G - 2)}{\sqrt{2(G-2)}}
```

Where: - $`T_{EF}`$ is the modified Pearson chi-square statistic - $`G`$
is the number of groups - The test statistic follows a standard normal
distribution under $`H_0`$

The null hypothesis is that the model fits the data adequately.

## Comparing with Different Group Numbers

The number of groups $`G`$ can affect the test’s performance:

``` r

# Test with different numbers of groups
group_sizes <- c(4, 8, 10, 15, 20)
results <- data.frame(
  Groups = group_sizes,
  P_value = sapply(group_sizes, function(g) {
    ef.gof(y, predicted_probs, G = g)$p_value
  })
)
print(results)
#>   Groups   P_value
#> 1      4 0.7449740
#> 2      8 0.3317745
#> 3     10 0.9344997
#> 4     15 0.7347885
#> 5     20 0.3532473
```

## Comparison with Hosmer-Lemeshow Test

Let’s compare the Ebrahim-Farrington test with the traditional
Hosmer-Lemeshow test:

``` r

# Hosmer-Lemeshow test (requires ResourceSelection package)
if (requireNamespace("ResourceSelection", quietly = TRUE)) {
  library(ResourceSelection)
  
  # Perform both tests
  ef_result <- ef.gof(y, predicted_probs, G = 10)
  hl_result <- hoslem.test(y, predicted_probs, g = 10)
  
  # Compare results
  comparison <- data.frame(
    Test = c("Ebrahim-Farrington", "Hosmer-Lemeshow"),
    P_value = c(ef_result$p_value, hl_result$p.value),
    Test_Statistic = c(ef_result$Test_Statistic, hl_result$statistic)
  )
  print(comparison)
} else {
  cat("ResourceSelection package not available for comparison\n")
}
#> ResourceSelection 0.3-6   2023-06-27
#>                         Test   P_value Test_Statistic
#>           Ebrahim-Farrington 0.9344997      -1.250567
#> X-squared    Hosmer-Lemeshow 0.9431075       2.855296
```

## New in 2.0.0: Directed test, ensemble, and the full battery

Version 2.0.0 turns the package into a full goodness-of-fit toolkit. The
**Directed EF (DEF)** test concentrates power on calibration-curve shape
directions,
[`def.ensemble.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md)
combines the DEF bases via the Cauchy combination test, and
[`run.all.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
runs a whole battery of GOF tests at once.

``` r

# Directed Ebrahim-Farrington test (takes the fitted model)
def.gof(model)                       # default poly3 basis
#>                          Test Basis Test_Statistic       df        Method
#> 1 Directed Ebrahim-Farrington poly3      0.1333555 2.058835 satterthwaite
#>     p_value
#> 1 0.9394923
def.gof(model, basis = "ensemble")   # combine all three bases (Cauchy)
#>           Test Combiner         Components k   p_value
#> 1 DEF ensemble      cct poly2+poly3+stukel 3 0.8921128

# Ensemble of the three DEF bases
def.ensemble.gof(model)
#>           Test Combiner         Components k   p_value
#> 1 DEF ensemble      cct poly2+poly3+stukel 3 0.8921128
def.ensemble.gof(model, add_ef = TRUE)   # add the omnibus EF
#>           Test Combiner            Components k   p_value
#> 1 DEF ensemble      cct poly2+poly3+stukel+EF 4 0.9070102
```

[`run.all.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
returns one tidy data frame, one row per test (here the fast battery;
the slow tests are shown below):

``` r

run.all.gof(model, include_slow = FALSE)
#> 
#> Goodness-of-fit battery: 21 tests  (0 reject at 0.05)
#> ========================================================
#>  Test                        Statistic   df  p-value    
#>  --- Global --------------------------------------------
#>  Pearson                         500.2  498   0.4640    
#>  Deviance                        548.4  498   0.0587 .  
#>  Information-Matrix            0.01797    2   0.9911    
#>  --- Standardized --------------------------------------
#>  Osius-Rojek                    0.1242        0.9012    
#>  McCullagh                     0.04384        0.4825    
#>  Copas-RSS                    0.001485    1   0.9988    
#>  EF                             -1.251    8   0.9345    
#>  EF-normal [a]                  -1.251    8   0.8945    
#>  --- Partition -----------------------------------------
#>  HL                              2.855    8   0.9431    
#>  HL-equalwidth                    4.86    8   0.7724    
#>  Pigeon-Heyse                    2.877    9   0.9690    
#>  F-test [b]                       1.11    9   0.3536    
#>  --- Covariate-space -----------------------------------
#>  Tsiatis                         5.778    9   0.7619    
#>  Xie                             8.312  8.5   0.4535    
#>  Pulkstenis-Robinson [c]                           -    
#>  --- Directed ------------------------------------------
#>  DEF.poly2                     0.06176 1.08   0.8227    
#>  DEF.poly3                      0.1334 2.06   0.9395    
#>  DEF.stukel                     0.3191    2   0.8313    
#>  Stukel                        0.03057    2   0.9848    
#>  --- Ensemble ------------------------------------------
#>  Ensemble.Vote(3DEF) [d]                      0.8921    
#>  Ensemble.Univ(3DEF+EF) [d]                   0.9070    
#> --------------------------------------------------------
#>  Signif.:  *** <.001   ** <.01   * <.05   . <.1
#>  Notes:
#>    [a] normal reference (thesis)
#>    [b] deviance residuals ~ groups (ANOVA F)
#>    [c] Not applicable: needs a categorical covariate
#>    [d] Cauchy combination of the directed tests
```

Add `include_slow = TRUE` to also run the opt-in slow tests (le Cessie,
the GAM-based tests, Stute-Zhu, eHL, BAGofT, and the Lai & Liu
standardized-power test), or pass `tests = c("EF", "DEF.poly3", "HL")`
to run a chosen subset.

### Powerful, but not liberal

Most GOF tests for logistic regression are **partition-based** (they
group the data and compare observed with expected counts), and that is
the family
[`ef.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md)
and
[`def.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md)
belong to. A key property of the directed tests is that they gain power
**without** inflating the type I error rate. In a Monte Carlo study (n =
500, 1000 replications, α = 0.05), the partition tests compare as
follows:

| Test | Size (null) | Power: quadratic | Power: wrong link |
|----|:--:|:--:|:--:|
| Hosmer–Lemeshow (decile) | 0.060 | 0.588 | 0.179 |
| Hosmer–Lemeshow (equal-width) | 0.053 | 0.332 | 0.244 |
| Pigeon–Heyse | 0.035 | 0.535 | 0.133 |
| EF (omnibus) | 0.058 | 0.480 | 0.218 |
| Tsiatis | 0.056 | 0.574 | 0.162 |
| Xie | 0.042 | 0.557 | 0.147 |
| DEF (poly3) | 0.060 | 0.709 | 0.404 |
| DEF (ensemble, vote) | 0.066 | 0.767 | 0.468 |

DEF and its vote ensemble are the most powerful in the family while
keeping the size near the nominal 0.05 — they are not liberal — and they
roughly double the power of Hosmer–Lemeshow, Tsiatis, and Xie on the
wrong-link misfit.

Partition-based tests are intuitive and work for sparse data and
continuous covariates (where the Pearson and deviance chi-square tests
fail), but their result depends on the grouping choice and the simpler
members (HL) can have low power. DEF keeps the intuitive
fitted-probability grouping but *directs* the test at calibration-curve
shapes, which is why it tops the table without losing size control.

> Note: as of 2.0.0,
> [`ef.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md)
> defaults to the chi-square reference (`method = "chisq"`); use
> `method = "normal"` for the version 1.0.0 behaviour.

## Power Analysis Example

How much does directing the test buy? The comparison below runs each
test on the *same* simulated datasets, so the difference between the
columns is the tests’ behaviour and not the noise between two separate
runs. The data carry a quadratic term; the fitted model omits it.

``` r

set.seed(2026)

power_paired <- function(n, beta_quad = 0.15, n_sims = 200, G = 10) {
  have_rs <- requireNamespace("ResourceSelection", quietly = TRUE)
  rej <- c(EF = 0, HL = 0)
  for (i in seq_len(n_sims)) {
    x    <- runif(n, -2, 2)
    y    <- rbinom(n, 1, plogis(x + beta_quad * x^2))
    fit  <- glm(y ~ x, family = binomial())      # misspecified: no quadratic term
    p    <- fitted(fit)

    # both tests see this same dataset
    if (ef.gof(y, p, G = G)$p_value < 0.05) rej["EF"] <- rej["EF"] + 1
    if (have_rs &&
        ResourceSelection::hoslem.test(y, p, g = G)$p.value < 0.05) rej["HL"] <- rej["HL"] + 1
  }
  c(n = n, EF = unname(rej["EF"]) / n_sims,
    HL = if (have_rs) unname(rej["HL"]) / n_sims else NA_real_)
}

power_results <- as.data.frame(do.call(rbind, lapply(c(200, 500, 1000), power_paired)))
print(power_results, row.names = FALSE)
#>     n    EF    HL
#>   200 0.070 0.065
#>   500 0.165 0.155
#>  1000 0.275 0.275
```

Read the table rather than a headline: on this misfit the two are close
at the smaller sample sizes, and any separation appears only as `n`
grows. The directed tests of the next section are where the gain against
a *shape* departure shows up; the omnibus EF test is built for breadth,
not for beating Hosmer-Lemeshow on a departure that grouping already
exposes.

## Handling Grouped Data (Original Farrington Test)

For datasets with grouped observations (multiple trials per covariate
pattern), you can use the original Farrington test:

``` r

# Simulate grouped data
set.seed(456)
n_groups <- 30
m_trials <- sample(5:20, n_groups, replace = TRUE)
x_grouped <- rnorm(n_groups)
prob_grouped <- plogis(0.2 + 0.8 * x_grouped)
y_grouped <- rbinom(n_groups, m_trials, prob_grouped)

# Create data frame and fit model
data_grouped <- data.frame(
  successes = y_grouped,
  trials = m_trials,
  x = x_grouped
)

model_grouped <- glm(
  cbind(successes, trials - successes) ~ x,
  data = data_grouped,
  family = binomial()
)

predicted_probs_grouped <- fitted(model_grouped)

# Original Farrington test for grouped data
result_grouped <- ef.gof(
  y_grouped,
  predicted_probs_grouped,
  model = model_grouped,
  m = m_trials,
  G = NULL  # No automatic grouping for original test
)

print(result_grouped)
#>                  Test Test_Statistic   p_value
#> 1 Farrington-Original      -1.476122 0.9300444
```

## Practical Guidelines

### When to Use Each Test Mode

1.  **Ebrahim-Farrington mode** (`G` specified):
    - Binary response data (0/1)
    - Want automatic grouping
    - Computationally efficient
    - Recommended for most applications
2.  **Original Farrington mode** (`m` provided, `G = NULL`):
    - Grouped binomial data
    - Multiple trials per covariate pattern
    - Requires fitted model object

### Choosing the Number of Groups

- **G = 10**: Standard choice, comparable to Hosmer-Lemeshow
- **G = 4-8**: For smaller datasets (n \< 200)
- **G = 20+**: For larger datasets (n \> 20000)
- **Rule of thumb**: Ensure each group sample size is large enough.

### Interpreting Results

- **p-value \> 0.05**: No evidence of lack of fit (fail to reject H₀)
- **p-value ≤ 0.05**: Evidence of model misspecification (reject H₀)
- **Very small p-values**: Strong evidence against model adequacy

## Advantages over Hosmer-Lemeshow Test

1.  **Better Power**: More sensitive to model misspecification
2.  **Theoretical Foundation**: Based on rigorous asymptotic theory
3.  **Sparse Data Handling**: Specifically designed for fully sparse
    data
4.  **Computational Efficiency**: Simplified calculations for binary
    data

## Limitations and Considerations

1.  **Group Selection**: Results can vary with different numbers of
    groups
2.  **Sample Size**: More reliable with larger sample sizes (n ≥ 100)
3.  **Model Complexity**: Performance with highly complex models needs
    further study

## References

1.  Farrington, C. P. (1996). On Assessing Goodness of Fit of
    Generalized Linear Models to Sparse Data. *Journal of the Royal
    Statistical Society. Series B (Methodological)*, 58(2), 349-360.

2.  Ebrahim, Khaled Ebrahim (2025). Goodness-of-Fits Tests and
    Calibration Machine Learning Algorithms for Logistic Regression
    Model with Sparse Data. *Master’s Thesis*, Alexandria University.

3.  Hosmer, D. W., & Lemeshow, S. (2000). Applied Logistic Regression,
    Second Edition. New York: Wiley.

4.  Hosmer, D. W., & Lemeshow, S. (1980). A goodness-of-fit test for the
    multiple logistic regression model. *Communications in Statistics -
    Theory and Methods*, 9(10), 1043–1069.
    <https://doi.org/10.1080/03610928008827941>

The Ebrahim-Farrington test provides a powerful and practical tool for
assessing goodness-of-fit in logistic regression, particularly for
binary data and sparse datasets. Its simplified implementation makes it
accessible for routine use while maintaining strong theoretical
foundations.
