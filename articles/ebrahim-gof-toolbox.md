# A goodness-of-fit and calibration toolbox for logistic regression

## Why a toolbox

Before the inferences or predictions of a logistic regression model can
be trusted, its *goodness of fit* has to be checked. Many tests exist
for this, but in practice they are **scattered across different packages
with different interfaces**, base R ships only the Hosmer–Lemeshow idea,
and most of the classical tests **lose power on sparse data** — the
common situation where continuous covariates give almost every
observation its own covariate pattern.

`ebrahim.gof` addresses this in one place. It:

- **introduces the author’s own sparse-data tests** — the omnibus
  Ebrahim–Farrington (EF) test, the *directed* EF / **EDGE** test that
  targets smooth calibration-shape departures, and a Cauchy-combination
  ensemble; and
- **aggregates a wide range of classical and modern tests**
  (Hosmer–Lemeshow, McCullagh, Osius–Rojek, le Cessie–van Houwelingen,
  Stute–Zhu, the binary-adaptive `BAGofT` test, and the `givitiR`
  calibration test) behind **one call**,
  [`run.all.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md),
  so they can be compared head to head. Each aggregated test is obtained
  from its own package and is credited to its authors; they are provided
  for comparison, not claimed as original here.

## The example data

The bundled `gof_demo` dataset was generated from a model whose true
linear predictor is **quadratic in age**. If we fit a model that treats
age *linearly*, it is therefore subtly misspecified through a smooth,
low-dimensional calibration distortion — exactly the misfit a *directed*
test is built to detect and an *omnibus* test tends to miss.

``` r

data("gof_demo", package = "ebrahim.gof")
str(gof_demo)
#> 'data.frame':    800 obs. of  5 variables:
#>  $ outcome  : int  1 0 0 0 0 0 0 0 0 1 ...
#>  $ age      : num  54.9 47.8 27 34.3 47.8 21.3 43.3 63.1 32.6 49 ...
#>  $ bmi      : num  32.3 29.2 22.8 26.2 23.4 19 34.6 33.3 27.6 21.3 ...
#>  $ sex      : int  1 0 1 0 0 1 0 0 1 0 ...
#>  $ treatment: int  0 0 1 0 0 1 1 0 0 0 ...

fit <- glm(outcome ~ age + bmi + sex + treatment,
           data = gof_demo, family = binomial)
```

## One call: the whole battery

[`run.all.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
runs the battery and returns a tidy `gof_battery` table. Here we use the
fast subset (`include_slow = FALSE`) and `install = "no"` so nothing is
installed or resampled while the vignette builds; tests whose optional
package is not present are simply skipped.

``` r

set.seed(1)
battery <- run.all.gof(fit, include_slow = FALSE, install = "no")
battery
#> 
#> Goodness-of-fit battery: 21 tests  (16 reject at 0.05)
#> ========================================================
#>  Test                        Statistic   df  p-value    
#>  --- Global --------------------------------------------
#>  Pearson                           760  795   0.8087    
#>  Deviance                        839.1  795   0.1352    
#>  Information-Matrix              26.39    5  7.5e-05 ***
#>  --- Standardized --------------------------------------
#>  Osius-Rojek                    -2.388        0.0169 *  
#>  McCullagh                      -2.725        0.9968    
#>  Copas-RSS                       3.352    1   0.0008 ***
#>  EF                              1.988    8   0.0431 *  
#>  EF-normal [a]                   1.988    8   0.0234 *  
#>  --- Partition -----------------------------------------
#>  HL                              15.52    8   0.0499 *  
#>  HL-equalwidth                   22.11    7   0.0024 ** 
#>  Pigeon-Heyse                    15.58    9   0.0761 .  
#>  F-test [b]                      2.088    9   0.0283 *  
#>  --- Covariate-space -----------------------------------
#>  Tsiatis                         63.53    9  2.8e-10 ***
#>  Xie                              51.9    7  6.1e-09 ***
#>  Pulkstenis-Robinson [c]         7.282    4   0.1217    
#>  --- Directed ------------------------------------------
#>  DEF.poly2                       8.299 1.08   0.0036 ** 
#>  DEF.poly3                        8.33 2.07   0.0150 *  
#>  DEF.stukel                        6.5 1.56   0.0118 *  
#>  Stukel                          14.24    2   0.0008 ***
#>  --- Ensemble ------------------------------------------
#>  Ensemble.Vote(3DEF) [d]                      0.0070 ** 
#>  Ensemble.Univ(3DEF+EF) [d]                   0.0088 ** 
#> --------------------------------------------------------
#>  Signif.:  *** <.001   ** <.01   * <.05   . <.1
#>  Notes:
#>    [a] normal reference (thesis)
#>    [b] deviance residuals ~ groups (ANOVA F)
#>    [c] split on categorical: sex, treatment
#>    [d] Cauchy combination of the directed tests
```

Read the table by *family*. The **global omnibus** statistics (Pearson,
deviance) spread their power over hundreds of covariate patterns and
tend to **miss** this smooth misfit. The **directed** tests —
`DEF.poly2`, `DEF.poly3` (EDGE), `DEF.stukel`, and Stukel’s score test —
concentrate their few degrees of freedom on the calibration-curve shape
and **detect** it, as do the EF test and the ensemble. That contrast is
the whole point of a directed test.

## Sparse or grouped? Let the package show you both

Three of the battery’s tests are *pattern-sensitive*: `Pearson`,
`Deviance`, and `McCullagh` have a reference distribution that depends
on whether they are computed per observation (the “sparse” form) or per
distinct *covariate pattern* (the “grouped” form). With all-continuous
covariates, as in `gof_demo`, every observation is its own pattern (800
patterns for 800 observations), so the two forms **coincide**; the
package recognises this degenerate case and omits the redundant
`(grouped)` rows, showing only the sparse form in the table above. (This
is also the regime where the per-observation Pearson and deviance
statistics lose their chi-squared reference entirely.) The paired
`(grouped)` rows appear only when covariate patterns actually repeat —
as they do next.

The bundled `gof_demo_grouped` dataset makes the distinction real. It
uses the *same* data-generating process, but `age` is recorded in
10-year bins and `bmi` as integers, so the 800 observations share only
328 covariate patterns. The battery detects this automatically and the
two forms now disagree — on this data the sparse deviance says “fine” (p
about 0.17) while the grouped deviance rejects (p about 0.04), with the
degrees of freedom dropping from 795 to 323:

``` r

data("gof_demo_grouped", package = "ebrahim.gof")
fit_g <- glm(outcome ~ age + bmi + sex + treatment,
             data = gof_demo_grouped, family = binomial)
set.seed(1)
run.all.gof(fit_g, include_slow = FALSE, install = "no")
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

Same model formula, same test family, opposite verdicts — so handling
the grouping wrong silently changes conclusions. The package’s answer is
to never make that choice for you invisibly: both forms are printed side
by side whenever they differ, and in truly sparse data (where no grouped
“correction” can rescue Pearson or the deviance) the partition and
directed families below are the tools that remain valid.

## The directed tests on their own

Each test is also callable directly. The EDGE test projects the grouped
standardized residuals onto a small basis of smooth calibration shapes:

``` r

edge.gof(fit)                 # directed EDGE test (poly3 basis)
#>   Test Basis Test_Statistic       df        Method    p_value
#> 1 EDGE poly3       8.329512 2.072814 satterthwaite 0.01496412
def.gof(fit, basis = "poly2") # a lower-order directed basis
#>                          Test Basis Test_Statistic       df        Method
#> 1 Directed Ebrahim-Farrington poly2        8.29914 1.076512 satterthwaite
#>       p_value
#> 1 0.003576024
ef.gof(fit)                   # the omnibus EF test
#>                 Test Test_Statistic    p_value
#> 1 Ebrahim-Farrington       1.987687 0.04309049
```

The Cauchy-combination ensemble pools the directed bases into a single
p-value:

``` r

def.ensemble.gof(fit)
#>           Test Combiner         Components k     p_value
#> 1 DEF ensemble      cct poly2+poly3+stukel 3 0.006959379
```

## Beyond the battery

[`run.all.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
selects among the classical tests and the directed ones. Several of the
package’s procedures sit outside it, because they need something the
battery does not have: a penalty, a frozen set of weights, or a scorer
of your own. They are called directly on the same fitted model.

### Penalised (ridge) fits

Under a penalty the classical references are not merely weak but
**invalid**: shrinkage biases the fitted probabilities, so the grouped
residuals acquire a non-centrality and the usual chi-squared reference
is wrong. Two functions apply the same correction and differ only in
what they refer it to.

``` r

set.seed(4)
n <- 300; p <- 8
X <- matrix(rnorm(n * p), n, p)
y <- rbinom(n, 1, plogis(0.2 + X %*% c(0.9, -0.6, 0.4, rep(0, p - 3))))

# closed form: one fit, no resampling, and no Monte Carlo error in the p-value
calm.gof(X, y, lambda = 40)
#> 
#> CALM: closed-form goodness of fit for penalized logistic regression
#> lambda = 40 (theory scale) = 0.1333 on glmnet's scale | G = 10 | p/n = 0.027
#>   SC.HL             statistic    6.189   p = 0.4904
#>   SC.EDGE           statistic    3.453   p = 0.1017
#>   SC.EDGE.adaptive  statistic    3.453   p = 0.1017   [degree 3, rho_hat 0.93]
#> 
#> SC.EDGE.adaptive is the statistic to prefer; see ?calm.gof for the scope of the reference.
```

[`shrink.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/shrink.gof.md)
computes the same statistics and calibrates them by prepivoting instead,
which costs `B` penalised refits. Prefer
[`calm.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md)
unless you need `penalize` to choose which columns are shrunk, or
`p >= n`, which it refuses. Note the penalty scale: `lambda` here is `n`
times a **glmnet** lambda, and
[`calm.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md)
takes a `lambda_scale` argument while
[`shrink.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/shrink.gof.md)
does not.

[`?calm.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md)
states where the reference is validated. It is not validated for
strongly unbalanced outcomes once `p/n` is an appreciable fraction, and
the function warns when you call it there.

### Frozen-weight and pretrained tests

[`legoft()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/legoft.md)
combines eleven statistics with weights fixed offline and shipped
frozen, so two analysts obtain the same p-value;
[`legoft.localize()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/legoft.localize.md)
then says which domain of evidence carries the misfit, with familywise
error control.
[`deepgof1()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md)
is a pretrained convolutional statistic whose level comes from your own
parametric bootstrap rather than from the network. None of these is
retrained when you call it.

### Building your own

[`gof.features()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof.features.md)
turns a fitted model into a feature vector of test statistics, and
[`deploy.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deploy.gof.md)
refers a scorer of yours to a parametric bootstrap under the fitted
model.

## Going further

- **The full battery.** Set `include_slow = TRUE` to add the
  resampling/refit-based tests (le Cessie–van Houwelingen, Stute–Zhu,
  `BAGofT`, GAM-smooth tests) and `calibration_plot = TRUE` for the
  `givitiR` calibration test and belt. These need their respective
  optional packages;
  [`gof_install_suggests()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof_install_suggests.md)
  installs them.
- **Speed.** The EDGE test is calibrated in closed form (no bootstrap,
  no refit), so it costs a few milliseconds even at large `n`. An
  optional GPU backend for the heavier O(n^(2)/O(n)3) aggregated tests
  (le Cessie, projection) is available in the package’s replication
  materials.
- **Choosing a test.** For the structured misfit that dominates practice
  — a wrong link or an omitted smooth term — prefer a directed test
  (EDGE) over an omnibus one; keep an omnibus test in reserve for rough,
  high-frequency departures.

## References

Farrington CP (1996). “On Assessing Goodness of Fit of Generalized
Linear Models to Sparse Data.” *Journal of the Royal Statistical Society
B*, **58**(2), 349–360.

Hosmer DW, Lemeshow S (1980). “Goodness of Fit Tests for the Multiple
Logistic Regression Model.” *Communications in Statistics – Theory and
Methods*, **9**(10), 1043–1069.
