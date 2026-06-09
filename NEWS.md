# ebrahim.gof 2.0.0

## New features

* `def.gof()` - the Directed Ebrahim-Farrington (DEF) goodness-of-fit test.
  Projects grouped standardized residuals onto a smooth calibration-shape basis
  (`"poly2"`, `"poly3"`, `"stukel"`) and calibrates the statistic as a weighted
  sum of chi-square_1 variables (Satterthwaite by default; Imhof via the
  suggested `CompQuadForm`). `basis = "ensemble"` is a shortcut to
  `def.ensemble.gof()`.
* `def.ensemble.gof()` - combines the three DEF bases (optionally the omnibus EF,
  or extra p-values) into one decision via the Cauchy combination test (default),
  with `minp` and `fisher` offered for comparison.
* `ef.gof()`, `def.gof()`, and `def.ensemble.gof()` now accept **either** a fitted
  `glm` **or** `(y, predicted_probs)` as input. For `def.gof`, supplying the design
  matrix `X` (with the `y`/`predicted_probs` form) gives the exact calibration;
  without it a conservative chi-square reference is used and a warning is issued.

## Breaking changes

* `ef.gof()` now defaults to the **chi-square** reference (`method = "chisq"`):
  the grouped statistic is referred to a chi-square_{G-2} distribution. Use
  `method = "normal"` to reproduce the previous (standardized-normal) p-value.

* `run.all.gof()` - a one-shot runner that returns a tidy `data.frame`, one row
  per test. Pass a fitted `glm` for the whole battery, or `(y, predicted_probs)`
  for the prediction-only tests. One failing test never aborts the run. This
  build bundles Pearson, Deviance, Osius-Rojek, Copas-RSS, Hosmer-Lemeshow
  (deciles and equal-width), Pigeon-Heyse, EF, the three DEF bases, Stukel, the
  covariate-space tests Tsiatis, Xie, and Pulkstenis-Robinson, and the two
  Cauchy-combination ensemble rows. Osius-Rojek, Copas-RSS, Pigeon-Heyse,
  Tsiatis, and Pulkstenis-Robinson were verified to match their original
  implementations to ~1e-15 (Xie's statistic also matches).
* All `run.all.gof()` tests were verified to reproduce the implementations used
  in the original thesis simulation. In particular `Osius-Rojek` and `Stukel` now
  follow `LogisticDx::gof.glm` (Stukel via `statmod::glm.scoretest`; `statmod`
  added to Suggests), matching it numerically; `Copas-RSS` matches `rms`'s gof
  residual; `HL` matches `ResourceSelection::hoslem.test`; and `HL-equalwidth`,
  Pigeon-Heyse, Tsiatis, Xie, and Pulkstenis-Robinson match their source scripts.
* A second EF row, `EF-normal`, reports the omnibus EF test with the normal
  reference used in the thesis simulation (the `EF` row uses the chi-square
  default).
* More opt-in slow (`include_slow = TRUE`) tests: the GAM-based `HL-GAM`,
  `PR-GAM`, and `Xie-GAM` (Xie et al. 2021; need `mgcv`; HL-GAM and PR-GAM match
  the source `gam_gof_tests` exactly, Xie-GAM uses a fixed clustering seed), and
  `Stute-Zhu` (cumulative-residual parametric-bootstrap test; sequential, set
  reps via `control = list("Stute-Zhu" = list(B = ...))`; statistic matches the
  source exactly).
* `Lai-Liu-HL` (Lai & Liu 2018, standardized-power procedure for the
  Hosmer-Lemeshow test). It has no p-value: it resamples to a target size, fits
  the model, estimates the HL rejection rate ("standardized power"), and returns
  a randomized accept/reject decision. The standardized power is reported as the
  statistic and the decision in the `Note` (set `n0`/`k` via `control`). Verified
  to match the source `lai_liu_test` exactly.
* Two further opt-in slow tests: `eHL` (the e-value Hosmer-Lemeshow test of Henzi
  et al. 2024; base-R reimplementation, with attribution, of the marius-cp/eHL
  code, matching it to ~1e-11; reported as `p = min(1, 1/e)`), and `BAGofT` (the
  binary-adaptive GOF test, wrapping the `BAGofT` package; set `nsim` via
  `control = list(BAGofT = list(nsim = ...))`).
* An opt-in slow test, `le-Cessie` (le Cessie-van Houwelingen 1995, general
  multivariate smoothed-residual test), runs when `include_slow = TRUE`. It is
  O(n^2)-O(n^3). Adapted with attribution from the USGS `smwrStats` package
  (public domain); verified to match it exactly.
* The `Xie` test uses the corrected degrees of freedom `G - k/2 - 1` with `k` the
  number of predictors. (Earlier thesis runs used `df = G - 0.5`, an artifact of
  `coef()` returning `NULL` on a predicted-probability list; the statistic is the
  same, only the p-value differs.)

* Added the `Information-Matrix` test (White 1982 / Orme 1988), the closed-form
  IM test; verified to match the thesis `IMtest_fast` exactly.

## Pending for the 2.0.0 release

* The remaining thesis tests are all slow / third-party and will be added as
  opt-in `include_slow = TRUE` tests in a later build: the GAM-based tests
  (HL-GAM, PR-GAM, Xie-GAM; need `mgcv`), the bootstrap tests (Hosmer bootstrap,
  Stute-Zhu), the e-value HL (`eHL`; needs `isotone`), and `BAGofT`. (McCullagh
  is not added: it appears only in the unused `goflogit` macro, not in the thesis
  simulation.)

# ebrahim.gof 1.0.0

## Initial Release

This is the first release of the ebrahim.gof package, implementing the Ebrahim-Farrington goodness-of-fit test for logistic regression models.

### Features

* **Main Function**: `ef.gof()` - Performs the Ebrahim-Farrington goodness-of-fit test
* **Dual Mode Support**:
  - Ebrahim-Farrington test with automatic grouping for binary data
  - Original Farrington test for grouped binomial data
* **Comprehensive Documentation**: Detailed help files and vignette
* **Robust Testing**: Extensive test suite with edge case handling
* **Input Validation**: Thorough parameter checking and error messages

### Key Capabilities

* **Binary Data**: Automatic grouping of binary (0/1) responses
* **Grouped Data**: Support for binomial data with multiple trials
* **Flexible Grouping**: User-specified number of groups (G)
* **Statistical Rigor**: Based on Farrington's (1996) theoretical framework
* **Sparse Data**: Optimized for sparse and challenging datasets

### Advantages over Existing Tests

* **Better Power**: More sensitive than Hosmer-Lemeshow test
* **Simplified Implementation**: Easy-to-use interface
* **Theoretical Foundation**: Rigorous asymptotic properties
* **Computational Efficiency**: Fast execution for binary data

### Technical Details

* **Test Statistic**: Uses modified Pearson chi-square with correction term
* **Distribution**: Standard normal under null hypothesis
* **Expected Value**: G - 2 for grouped binary data
* **Variance**: 2(G - 2) for grouped binary data

### References

* Farrington, C. P. (1996). On Assessing Goodness of Fit of Generalized Linear Models to Sparse Data. *Journal of the Royal Statistical Society. Series B (Methodological)*, 58(2), 349-360.
* Ebrahim, Khaled Ebrahim (2025). Goodness-of-Fits Tests and Calibration Machine Learning Algorithms for Logistic Regression Model with Sparse Data. *Master's Thesis*, Alexandria University.

### Author

Ebrahim Khaled Ebrahim (Alexandria University)
Email: ebrahimkhaled@alexu.edu.eg 