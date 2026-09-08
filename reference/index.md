# Package index

## Start here

The battery runs everything and groups the results by the departure each
test detects. The package page has a table for choosing a single test.

- [`run.all.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
  : Run a Battery of Goodness-of-Fit Tests at Once
- [`ebrahim.gof-package`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ebrahim.gof-package.md)
  [`ebrahim.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ebrahim.gof-package.md)
  : ebrahim.gof: Goodness-of-Fit and Calibration Tests for Logistic
  Regression

## Sparse data

Where the classical Hosmer-Lemeshow test loses power. The omnibus test
asks whether anything is wrong; the directed tests spend their few
degrees of freedom on the shapes where structured misfit concentrates.

- [`ef.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/ef.gof.md)
  : Ebrahim-Farrington Goodness-of-Fit Test for Logistic Regression
- [`edge.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.gof.md)
  : EDGE: Directed Goodness-of-Fit Test for Binary Logistic Regression
- [`edges.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edges.gof.md)
  : EDGES: Cauchy-Combination Ensemble of Directed GOF Tests (Alias)
- [`cdef.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/cdef.gof.md)
  : Covariate-Space Directed Ebrahim-Farrington (CDEF) Goodness-of-Fit
  Test
- [`def.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md)
  : Directed Ebrahim-Farrington (DEF) Goodness-of-Fit Test
- [`def.ensemble.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.ensemble.gof.md)
  : Combine Directed GOF Tests into One Decision (Ensemble)

## Penalized (ridge) fits

Under a penalty the usual chi-squared references are not merely weak but
invalid, because shrinkage displaces the grouped residuals. These two
apply the same correction and differ only in the reference they refer it
to.

- [`calm.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/calm.gof.md)
  : Closed-form goodness-of-fit test for penalized logistic regression
  (CALM)
- [`shrink.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/shrink.gof.md)
  : Shrinkage-corrected goodness-of-fit test for penalized (ridge)
  logistic regression

## Pretrained and frozen-weight tests

Weights fixed offline and shipped frozen, so two analysts obtain the
same p-value. Nothing is retrained at call time.

- [`legoft()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/legoft.md)
  : LEGofT: frozen-weight combination goodness-of-fit test for binary
  logistic regression
- [`legoft.localize()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/legoft.localize.md)
  : Localize misspecification with familywise error control
- [`deepgof1()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md)
  : DeepGOF-1: a pretrained goodness-of-fit test for logistic regression

## Building your own scorer

- [`gof.features()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof.features.md)
  : Goodness-of-fit evidence features for a fitted model
- [`deploy.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deploy.gof.md)
  : Deployable learned-ensemble GOF test via parametric bootstrap

## Data and setup

- [`gof_demo`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof_demo.md)
  : Synthetic binary outcome data with a smooth calibration misfit

- [`gof_demo_grouped`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof_demo_grouped.md)
  :

  Grouped-covariate companion to `gof_demo` (replicated covariate
  patterns)

- [`gof_install_suggests()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof_install_suggests.md)
  :

  Install the optional packages used by
  [`run.all.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)

## Printing and plotting

- [`print(`*`<gof_battery>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/print.gof_battery.md)
  : Print a goodness-of-fit battery
- [`plot(`*`<gof_battery>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/plot.gof_battery.md)
  : Plot the GiViTI calibration belt from a goodness-of-fit battery
