# Cubic Calibration Likelihood-Ratio Test for Logistic Regression

`cubic.calib.gof()` asks whether the calibration curve of a binary risk
model bends on the logit scale. It fits two logistic regressions of the
outcome on the logit of the model's own predicted risk, \\\hat\eta\\:
one with \\\hat\eta\\ alone, and one that adds \\\hat\eta^2\\ and
\\\hat\eta^3\\. Twice the drop in deviance is referred to a chi-squared
distribution on 2 degrees of freedom.

## Usage

``` r
cubic.calib.gof(object, predicted_probs = NULL)
```

## Arguments

- object:

  Either a fitted `glm` with `family = binomial` and a 0/1 response, or
  a numeric 0/1 outcome vector, in which case `predicted_probs` is
  required.

- predicted_probs:

  Predicted probabilities for `object` when it is an outcome vector, for
  example from an external model. Ignored when `object` is a `glm`.

## Value

A one-row `data.frame` with columns `Test` (`"Cubic-LR"`),
`Test_Statistic` (the likelihood-ratio statistic), `df` (2, or fewer if
a column cannot be estimated), `Method` (`"likelihood ratio"`) and
`p_value`. A fit that fails or does not converge gives an `NA` p-value
with a warning.

## Details

**Where the test comes from.** This is not a new test, and the package
does not claim it. Testing a fitted model by adding powers of its own
fitted values goes back to Tukey's (1949) one degree of freedom for
non-additivity, which adds the square of the fitted values, and Ramsey
(1969) made it a general specification test for linear regression,
RESET, by adding their square and cube. Pregibon (1980) carried the idea
to generalized linear models as a goodness-of-link test: refit the model
with constructed variables computed from its fitted linear predictor,
and read the change in deviance. Morgan (1985) studied the logistic
model whose linear predictor is a cubic polynomial, the cubic logistic
model, whose extra terms bend the two tails of the logistic curve.
`cubic.calib.gof()` is this family's test for a logistic model: the
constructed variables are \\\hat\eta^2\\ and \\\hat\eta^3\\, and the
test is the likelihood-ratio test of their two coefficients. The GiViTI
calibration belt (Finazzi et al., 2011; Nattino et al., 2014) fits the
same kind of polynomial in the logit of the predicted risk, choosing its
degree from the data and testing all of its coefficients; this function
fixes the degree at three and tests only the curvature.

**What it tests.** When the model was fitted by maximum likelihood to
the same data, the logistic regression of `y` on \\\hat\eta\\ reproduces
the fit exactly (intercept 0, slope 1), so the test asks only whether
the calibration curve is curved on the logit scale: the quadratic term
picks up an asymmetric bend, as from a complementary log-log truth
fitted by a logit, and the cubic term a symmetric one, as from tails
that are too heavy or too light. With the predictions of an external
model (`predicted_probs`) the intercept and slope are free as well, so
the test still asks only about curvature; calibration-in-the-large and
the calibration slope need their own tests (Cox, 1958).

**What to know before using it.** The test weighs every record by its
own linear predictor. A few records with corrupted covariates, whose
predictions are pushed to the extremes, can therefore drive its
false-alarm rate far above the nominal level; the grouped directed test
[`edge.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.gof.md)
is protected against that. The polynomial is fitted on the standardized
\\\hat\eta\\, which spans the same polynomials and leaves the statistic
unchanged but keeps the fit well conditioned.

## References

Cox, D. R. (1958). Two further applications of a model for binary
regression. *Biometrika*, 45(3-4), 562-565.
[doi:10.1093/biomet/45.3-4.562](https://doi.org/10.1093/biomet/45.3-4.562)

Finazzi, S., Poole, D., Luciani, D., Cogo, P. E. and Bertolini, G.
(2011). Calibration belt for quality-of-care assessment based on
dichotomous outcomes. *PLoS ONE*, 6(2), e16110.
[doi:10.1371/journal.pone.0016110](https://doi.org/10.1371/journal.pone.0016110)

Morgan, B. J. T. (1985). The cubic logistic model for quantal assay
data. *Journal of the Royal Statistical Society, Series C (Applied
Statistics)*, 34(2), 105-113.
[doi:10.2307/2347362](https://doi.org/10.2307/2347362)

Nattino, G., Finazzi, S. and Bertolini, G. (2014). A new calibration
test and a reappraisal of the calibration belt for the assessment of
prediction models based on dichotomous outcomes. *Statistics in
Medicine*, 33(14), 2390-2407.
[doi:10.1002/sim.6100](https://doi.org/10.1002/sim.6100)

Pregibon, D. (1980). Goodness of link tests for generalized linear
models. *Journal of the Royal Statistical Society, Series C (Applied
Statistics)*, 29(1), 15-23.
[doi:10.2307/2346405](https://doi.org/10.2307/2346405)

Ramsey, J. B. (1969). Tests for specification errors in classical linear
least-squares regression analysis. *Journal of the Royal Statistical
Society, Series B*, 31(2), 350-371.
[doi:10.1111/j.2517-6161.1969.tb00796.x](https://doi.org/10.1111/j.2517-6161.1969.tb00796.x)

Tukey, J. W. (1949). One degree of freedom for non-additivity.
*Biometrics*, 5(3), 232.
[doi:10.2307/3001938](https://doi.org/10.2307/3001938)

## See also

[`edge.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.gof.md)
for the grouped directed test,
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md),
where this test is the row `"Cubic-LR"`.

## Author

Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>

## Examples

``` r
set.seed(1)
x <- runif(1000, -3, 3)
y <- rbinom(1000, 1, plogis(0.6 * x))
fit <- glm(y ~ x, family = binomial())
cubic.calib.gof(fit)                             # correct model
#>       Test Test_Statistic df           Method   p_value
#> 1 Cubic-LR      0.4651523  2 likelihood ratio 0.7924894

y2 <- rbinom(1000, 1, 1 - exp(-exp(0.6 * x)))       # complementary log-log truth
cubic.calib.gof(glm(y2 ~ x, family = binomial()))   # fitted by a logit: the curve bends
#>       Test Test_Statistic df           Method    p_value
#> 1 Cubic-LR       14.32408  2 likelihood ratio 0.00077547

## frozen predictions, e.g. an external validation
cubic.calib.gof(y, predicted_probs = fitted(fit))
#>       Test Test_Statistic df           Method   p_value
#> 1 Cubic-LR      0.4651523  2 likelihood ratio 0.7924894
```
