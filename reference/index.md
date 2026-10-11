# Package index

## Start here

The battery runs everything and groups the results by the departure each
test detects. The package page has a table for choosing a single test.

- [`run.all.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
  : Run a Battery of Goodness-of-Fit Tests at Once
- [`run.all.external()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.external.md)
  : Run the External-Validation Tests at Once
- [`edge.stream()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.stream.md)
  [`update(`*`<edge_stream>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.stream.md)
  [`summary(`*`<edge_stream>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.stream.md)
  [`print(`*`<edge_stream>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.stream.md)
  : EDGE for Streaming Data: Calibration Monitoring Without
  Recomputation
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
- [`edge.belt()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.belt.md)
  [`print(`*`<edge_belt>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.belt.md)
  [`plot(`*`<edge_belt>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/edge.belt.md)
  : EDGE-Belt: the grouped calibration picture behind the EDGE test
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
- [`deepgof1.external()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.external.md)
  : DeepGOF-1 for frozen predictions: external validation of a risk
  model
- [`gof_certificate()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/gof_certificate.md)
  : Can DeepGOF-1 detect this departure? A certificate from the frozen
  weights

## Localizing the misfit

Which part of the misfit is present – the intercept, the slope, the link
or the covariates – with familywise error control, and the update it
calls for. For frozen predictions on validation data, or a fitted
logistic model.

- [`localize.external()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md)
  : Localize the misfit of frozen predictions, with error control
- [`localize.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.gof.md)
  : Localize the misfit of a fitted logistic model, with error control

## Comparison tests

Published tests implemented here so the battery can be compared with
them on the same data: the projection test of Liu et al. (2024), a fast
implementation of the adaptive-partition test BAGofT, and the cubic
calibration likelihood-ratio test.

- [`projection.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/projection.gof.md)
  : Projection Goodness-of-Fit Test for Binary Regression
- [`bagoft.fast()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/bagoft.fast.md)
  : Fast BAGofT: the Binary Adaptive Goodness-of-Fit Test
- [`cubic.calib.gof()`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/cubic.calib.gof.md)
  : Cubic Calibration Likelihood-Ratio Test for Logistic Regression

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
- [`summary(`*`<gof_battery>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/summary.gof_battery.md)
  [`print(`*`<summary.gof_battery>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/summary.gof_battery.md)
  : Summarise a goodness-of-fit battery by family
- [`plot(`*`<gof_battery>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/plot.gof_battery.md)
  : Plot the GiViTI calibration belt from a goodness-of-fit battery
- [`plot(`*`<edge_gof>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/plot.edge_gof.md)
  : Plot an EDGE result: the grouped residuals, the shape tested and the
  bands
- [`plot(`*`<deepgof1>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/plot.deepgof1.md)
  : Plot the residual map of a DeepGOF-1 test
- [`plot(`*`<gof_localize>`*`)`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/plot.gof_localize.md)
  : Plot a localization verdict
