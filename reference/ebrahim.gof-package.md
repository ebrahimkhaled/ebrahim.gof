# ebrahim.gof: Goodness-of-Fit and Calibration Tests for Logistic Regression

A unified toolbox of goodness-of-fit and calibration tests for binary
logistic regression, callable in a single line via
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md).
The package is aimed particularly at *sparse* data, where the classical
Hosmer–Lemeshow test loses power, and at *penalized* fits, where it is
not merely weak but invalid.

## Choosing a test

The battery is the place to start; the table below is for when you
already know something about what you are looking for.

- No prior idea of what is wrong:

  [`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
  — runs the whole battery and groups the results by the departure each
  test detects.

- Sparse data, no direction in mind:

  [`ef.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md)
  — the omnibus Ebrahim–Farrington test, which groups automatically and
  needs no model object.

- Misfit expected in the shape of the calibration curve:

  [`edge.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.gof.md)
  — spends its few degrees of freedom on the smooth directions where
  structured misfit concentrates.

- Misfit expected in the covariates themselves:

  [`cdef.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/cdef.gof.md)
  — the covariate-space directed test.

- Unwilling to choose one direction:

  [`edges.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edges.gof.md)
  — a Cauchy combination over the directed bases, which pays little for
  the directions that turn out to be empty.

- Wanting one number, reproducible between analysts:

  [`legoft`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/legoft.md)
  — weights fixed offline and shipped frozen, so nothing is retrained
  when you call it;
  [`legoft.localize`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/legoft.localize.md)
  then says which domain of evidence carries the misfit, with familywise
  error control.

- Willing to spend a bootstrap for more power:

  [`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md)
  — a pretrained convolutional statistic whose level comes from your own
  parametric bootstrap rather than from the network.

- A penalized (ridge) fit:

  [`calm.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md)
  for a closed-form reference from a single fit, or
  [`shrink.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/shrink.gof.md)
  for the same correction referred to a prepivoting bootstrap. Under a
  penalty the usual chi-squared references are wrong, not just
  conservative.

- A scorer of your own:

  [`gof.features`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof.features.md)
  turns a fit into a feature vector and
  [`deploy.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deploy.gof.md)
  applies a scorer to it.

## Aggregated tests (for comparison)

[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
also runs, in one call, a wide range of classical and modern tests —
Hosmer–Lemeshow, McCullagh, Osius–Rojek, le Cessie–van Houwelingen,
Stute–Zhu, the binary-adaptive BAGofT test, and the givitiR calibration
test. Each aggregated test is obtained from its own package (where
installed) and is attributed to its authors; these are provided for
head-to-head comparison, not claimed as original to this package.
[`gof_install_suggests`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof_install_suggests.md)
installs the optional packages they need.

## Data

[`gof_demo`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof_demo.md)
is a bundled example dataset with a documented, reproducible misfit for
illustrating the battery.

## Citing the methods

Each of the author's tests has a paper behind it. Run
`citation("ebrahim.gof")` for the current references; they are kept
there rather than duplicated here, because several are moving from
preprint to journal.

## See also

The vignette
[`vignette("ebrahim-gof-toolbox", package = "ebrahim.gof")`](https://ebrahimkhaled.github.io/ebrahim.gof/articles/ebrahim-gof-toolbox.md).

## Author

**Maintainer**: Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>
([ORCID](https://orcid.org/0009-0006-7839-8778))

## Examples

``` r
set.seed(1)
n <- 200
x <- rnorm(n)
y <- rbinom(n, 1, plogis(0.3 + 0.9 * x))
fit <- glm(y ~ x, family = binomial())

# the omnibus test on the fitted probabilities
ef.gof(y, fitted(fit))
#>                 Test Test_Statistic   p_value
#> 1 Ebrahim-Farrington     -0.7714014 0.7666854

# the directed test, when misfit is expected in the calibration shape
edge.gof(fit)
#>   Test Basis Test_Statistic       df        Method   p_value
#> 1 EDGE poly3       2.642516 2.057006 satterthwaite 0.2717296
```
