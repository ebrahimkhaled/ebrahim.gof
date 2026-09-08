# ebrahim.gof 2.7.0

## Bug fix

* The grouped-data example on `?ef.gof` now passes `G = NULL`, and runs. Without it the
  call fell through to the automatic-grouping branch, which ignores `m` and `model` and
  refers binomial counts to the binary statistic, so the documented example reported a
  p-value of zero on data drawn from the fitted model. The example is no longer wrapped
  in `\dontrun{}`, so `R CMD check` exercises the original Farrington branch as well.
  The `@note` bullet that said supplying `m` selects the original test has been corrected
  to say that `G` must be set to `NULL`.

* `gof.features()` honoured its documented `tests` argument only after the fact: it ran
  the whole battery and subset the result, so a caller asking for nine p-values paid for
  all twenty-five. The returned feature vector is unchanged; only the cost was wrong.

* `print.shrink.gof()` was defined but never registered, so `shrink.gof()` results printed
  as a raw nested list. It is now an S3 method, and reports `p/n` as `print.calm.gof()`
  already did.

## Documentation

* The `Description` field now names `calm.gof()`, `edges.gof()` and `cdef.gof()`, and gives
  the penalized case its own sentence rather than filing `shrink.gof()` under sparse data.
  The field is frozen for the life of a release and is the text CRAN's own search indexes,
  so an omission there costs months; this is the same slip recorded at 2.5.0.

* `Authors@R` used positional arguments, which made "Khaled Ebrahim" the family name, so
  `citation()` rendered "Khaled Ebrahim E" and BibTeX filed the package under K. The
  arguments are now named and the surname agrees with `inst/CITATION`.

* `?shrink.gof` had two roxygen drafts merged, which left every parameter documented twice
  and the real title and description buried inside `\value`, so the rendered page's
  description was its own title repeated. Rewritten.

* `?calm.gof` now documents its class-balance scope, and `calm.gof()` warns when it is
  called on a strongly unbalanced outcome, where the reference is not validated.

* `?legoft` no longer carries a title indistinguishable from `?deepgof1`.

* The package landing page, `?ebrahim.gof`, is no longer marked internal and now opens
  with a table for choosing among the tests.

* Every function page gained `\concept` entries, which is what `??` and the documentation
  mirrors search.

* `inst/CITATION` took its version from a hardcoded string, which had been stale since
  2.4.0; it now reads `meta$Version`. The benchmark paper is recorded as in press at the
  Journal of Intelligent Computing and Data Science rather than as a preprint.

* The full-battery example forwards `nsim = 20` to BAGofT, which cuts
  `R CMD check --run-donttest` from 213 to 54 seconds without changing what it shows.


## New features

* `calm.gof()` -- a closed-form reference distribution for the shrinkage-corrected
  goodness-of-fit statistics, so that they can be used without a bootstrap. It is the
  companion of `shrink.gof()`: the statistics are the same, only the reference differs,
  and where `shrink.gof()` needs several hundred penalized refits `calm.gof()` needs one
  fit and returns in a fraction of a second.

  The reference is built by replacing the fitted Bernoulli variance, which shrinkage
  inflates towards one quarter, by a de-noised estimate of the null variance obtained
  from the fit alone through the observable adjustments of Bellec (2025), and then
  reading the exact tail of the resulting weighted chi-squared law (Davies 1980).

  Three statistics are returned. `SC.EDGE.adaptive` is the one to prefer: it projects
  the corrected residual onto orthogonal polynomials in the group-mean fitted
  probability and chooses the degree from the data, keeping the cubic direction only
  when the fitted index is accurate enough to carry a cubic signal. `SC.HL` is the
  shrinkage-corrected Hosmer-Lemeshow statistic, reported for continuity with that
  tradition; its reference relies on a constant calibrated by simulation, which does
  not transfer to every design, and `?calm.gof` says where it fails.

  The reference is validated for aspect ratios p/n up to 0.4 and for penalties that
  shrink towards zero. The lasso is not covered.

## Dependency change

* `CompQuadForm` moves from Suggests to Imports. `calm.gof()` cannot produce a p-value
  without it, so it is no longer optional.

# ebrahim.gof 2.6.0

## New features

* `deepgof1()` -- DeepGOF-1, a pretrained goodness-of-fit test for binary logistic
  regression whose statistic is a convolutional network. The network reads the fitted
  model's *residual map* (a 6x6 grid of standardized residual sums over the ranks of the
  two strongest covariates), so misfit is detected by its spatial pattern: an omitted
  quadratic paints a stripe, an omitted interaction a saddle, a local departure a bump.

  The network was trained once, offline, on simulated departures and ships frozen as
  18,273 numbers. Nothing is trained when you call the function, and the p-value is the
  rank of the observed score inside your own parametric bootstrap -- so the level is a
  property of the calibration, not of what the network learned. A bootstrap refit that
  fails is scored `+Inf`, counting against rejection.

  Implemented in pure base R (the forward pass is a few small matrix multiplies), so the
  package still needs no Python and no additional dependency. The shipped weights are
  verified against the training framework in `tests/testthat/test-deepgof1.R`; agreement
  is to 5e-11.

  It is a small-sample instrument: on a published benchmark it outpowers every classical
  partition test at every sample size, most clearly at n = 50 to 200, while tests that use
  the whole covariate space rather than a two-covariate grid do better overall. The
  training corpus, training code, benchmark harness and every per-replicate p-value behind
  those numbers are archived separately from this package.

# ebrahim.gof 2.5.0

## New features

* `legoft()` -- a pretrained goodness-of-fit test for binary logistic regression. It
  combines eleven classical and directed statistics with weights that were fixed offline
  and ship frozen, and calibrates the combination by a parametric bootstrap at the fitted
  parameters. Nothing is retrained when you call it, so two analysts running it on the
  same data obtain the same p-value.

  On twelve held-out misspecification families the rule attained higher mean power than
  each of the fourteen individual tests in the study, at n = 500 and n = 1000. It does
  **not** dominate them family by family: against the two covariate-space directed tests
  the advantage is in the mean rather than in per-family consistency, which is the
  dilution a combination pays for having no catastrophic blind spot.

  It does **not** beat an equal-weight Cauchy combination of the same eleven members
  (Liu and Xie, 2020). The difference is +0.0045 at n = 500 (sign test p = 0.23) and
  -0.0005 at n = 1000. Users who want the simpler rule lose nothing measurable; `legoft()`
  reduces to it exactly at the boundary of its weight family, so it cannot do worse.

* `shrink.gof()` -- goodness of fit for **penalized (ridge) logistic regression**. The
  Hosmer--Lemeshow test is not valid when the coefficients are shrunk: penalization biases
  the fitted probabilities, the grouped residuals acquire a non-centrality, and the usual
  chi-squared reference is wrong. This removes the estimated non-centrality and refers the
  corrected statistic to a bootstrap built from the debiased generator (Beran prepivoting),
  on either the decile or the EDGE basis.

  Note the scale: `lambda` is on the theory scale, `lambda = n * lambda_glmnet`. Passing a
  \pkg{glmnet} lambda unchanged is the commonest way to misuse this function.

  The implementation is vendored byte-identical from the source of Ebrahim (2026),
  "Shrinkage invalidates the Hosmer--Lemeshow test", so the paper and the package compute
  the same numbers. It depends only on base R and stats.

* `legoft.localize()` -- reports which of two domains of evidence carries the misfit,
  using closed testing, so the probability of implicating any collection of domains whose
  pooled evidence is jointly exchangeable with the reference draws' is at most `alpha`.

  The domains are defined by what the members can see rather than by taxonomy. `INDEX`
  holds the seven statistics that read the linear predictor -- the grouping tests, which
  stratify on fitted risk, and the directed tests, which examine bends in that same index.
  `COV` holds the four that read directions orthogonal to it. Calibration and link
  readings are pooled deliberately: they are **not separately identifiable**, because
  fitted risk is a monotone transform of the index, so both read one axis.

  A verdict locates the evidence. It does not name the repair: a departure of one kind
  can move members assigned to the other domain, and an unimplicated domain is not
  thereby certified correct.

## Documentation

* The help for `run.all.gof()` now describes each test individually. It previously named
  around twenty-five tests in a single paragraph with a parenthetical each, which was
  enough to tell you a test existed but not enough to decide whether to run it or what a
  rejection from it meant. Entries are now grouped by the mechanism the test uses --
  global and standardized, partition, directed, covariate-space, smoothing, resampling,
  calibration, combinations -- and each says what the test computes, which departure it
  is built to notice, and where it fails.

  The failure modes are stated because they are the part users get wrong. Pigeon--Heyse is
  conservative in sparse designs; Stukel's two-parameter form does not always hold its
  nominal level; `HL-equalwidth` is frequently not computable when fitted risks are
  concentrated in a narrow band; GiViTI's internal and external forms are not
  interchangeable; and the chi-squared reference with `G - 2` degrees of freedom was
  established by simulation rather than derived, so `G = 10` is a convention and worth
  varying.

* Ten classical tests were implemented but never cited. Pigeon--Heyse, Copas, White,
  Orme, Kuss, Lai--Liu, Nattino (GiViTI), Zhang (BAGofT), Liu--Xie (Cauchy combination)
  and Hosmer et al. (1997) now appear in the references with resolved DOIs.

* The package's own methods now point somewhere a reader can follow: the arXiv preprints
  for the Ebrahim--Farrington test, the benchmark study and the thesis, and the archived
  reproduction materials for the EDGE test, the EDGES ensemble, the detection-subspaces
  framework and `shrink.gof()`.

* The examples show how to *read* the returned panel -- subsetting by p-value, counting
  by family, contrasting a correctly specified model against one with an omitted quadratic
  term, varying `G` -- rather than only how to call the function.

* The `Description` field predated `legoft()` and `shrink.gof()` and mentioned neither.

## Performance

* `gof_lecessie()` is much faster at the sample sizes where it was unusable, with
  **identical numerics**.
  Nothing about the test changed: not the statistic, not the moment reference, not the
  degrees of freedom, not the p-value. Two exact algebraic identities replaced two
  matrix products that were costing O(n^3):

  - The moment reference (I-H)'R(I-H) is now assembled from the rank-p factors of
    H = VX(X'VX)^{-1}X' instead of forming H as an n-by-n matrix and multiplying it
    out, which is O(n^2 p).
  - The variance term 2 tr(MVMV) is evaluated as 2 sum_ij M_ij^2 mu2_i mu2_j, which is
    O(n^2). This is not a new identity: le Cessie and van Houwelingen (1995) state the
    variance in exactly that elementwise form in eq. (A.6) before collapsing it to the
    trace. The package had been computing the trace version.

  Measured against the previous implementation, both byte-compiled and called through
  the installed package: 2.5x at n = 200, 4.5x at n = 500, 11.0x at n = 1000, 42.6x at
  n = 2000 and 87.7x at n = 3000 (48.1 s down to 1.1 s, and 139.9 s down to 1.6 s). The
  ratio grows with n because the change is O(n^3) to O(n^2 p) rather than a constant
  factor; at small n the shared cost of building the kernel matrix, which is unchanged,
  still dominates. The largest relative difference in the statistic was 6e-14 across
  null models, misspecified models, factor covariates and designs whose fitted
  probabilities reach machine zero. Rejection rates agree to four decimal places under
  the null and under alternatives, so no previously reported result moves.
  `tests/testthat/test-lecessie-algebra.R` pins the agreement at 1e-12 relative against
  the previous code, kept verbatim.

  Memory is unchanged: the test still builds the n-by-n kernel matrix and is still
  O(n^2) in space. This change buys time, not memory.

* The comment attached to the 2.4.1 transpose fix was wrong about the source and has
  been corrected. It said `smwrStats::leCessie.test()` was written for ordinary least
  squares, where the hat matrix is symmetric and the transpose is a no-op. It is not:
  that function computes the weighted, non-symmetric H correctly and then omits the
  transpose, so it departs from le Cessie and van Houwelingen (1995), Section 4, which
  prescribes (I-H)'R(I-H) in words. The omission is a size bug rather than a power bug.
  At n = 200 over 4000 replicates, the smwrStats form rejects a true null 7.05% of the
  time at the 5% level against 5.55% for the corrected form, and 7.62% against 5.15%
  when the fitted probabilities are extreme; power against a quadratic or an
  interaction departure is 1.000 either way.

## Notes on what is *not* here

* A ridge rule shrunk toward the equal-weight combination ("shrink-to-CCT") was one of
  six candidates evaluated during development and was **not selected**: it ties the
  equal-weight rule on most training folds and loses on oscillating departures (0.820
  against 0.990). Its appeal was that it reduces exactly to the equal-weight rule in the
  limit; `legoft()` has that same property at the boundary of its own weight family, so
  nothing is lost by omitting it. The candidate results are archived with the manuscript.

# ebrahim.gof 2.4.1

## Bug fix

* `le-Cessie` (the le Cessie--van Houwelingen smoothed-residual test) computed its
  moment reference as `(I - H) R (I - H)`, where the residual expansion requires
  `(I - H)' R (I - H)`. In a weighted fit `H = VX(X'VX)^{-1}X'` is idempotent but
  **not symmetric**, so the transpose matters; the inherited implementation came from
  the unweighted linear-model setting, where the two forms coincide. Only the
  reference moments (`E[Q]`, `Var[Q]`, and hence the matched degrees of freedom and
  the p-value) were affected -- the raw quadratic form `Q` was always correct.

  **Wording corrected in 2.5.0.** The 2.4.1 note said "the statistic was always correct".
  That is true of the raw `Q` and misleading about the `Statistic` column users actually
  see: that column reports `Q` rescaled by the matched moments, `Q * 2E[Q]/Var[Q]`, so it
  moves with the correction exactly as the degrees of freedom and the p-value do. On the
  package's own regression example the reported statistic goes from 16.794627 to
  15.895427 and the p-value from 0.149928 to 0.177213 (the matched degrees of freedom
  are 11.648746 after the fix). Expect the printed statistic and degrees of freedom to
  change too, not only the p-value.
  The consequence was a mildly liberal test: at `n = 1000` with the default bandwidth
  the empirical null size was about 0.063 instead of the exact-moment value 0.054.
  Users comparing results across versions should expect slightly larger p-values
  (hence slightly lower power) from `le-Cessie` after this fix.

# ebrahim.gof 2.4.0

## Print method

* `print.gof_battery()` redesigned for readability and width-robustness: tests are
  grouped under `--- Family ---` separators, per-test notes moved to compact `[a]`
  footnotes (so the table no longer wraps on narrow terminals), numbers are
  right-aligned, and the header reports how many tests reject at 0.05.

## Battery behaviour

* `run.all.gof()` now omits a `"(grouped)"` row (Pearson / deviance / McCullagh)
  when it is identical to the sparse form -- i.e., when no covariate pattern
  repeats -- so fully sparse data no longer shows duplicate rows. The grouped row
  is still reported when patterns repeat and the two forms differ.
* The `BAGofT` control entry now forwards the test's key tuning parameters:
  `nsim`, `nsplits`, `ne`, and the random-forest partitioner's `Kmax` (maximum
  number of adaptive partition cells), `ntree`, `nmin`, `mtry`, `maxnodes`.


## New features

* `run.all.gof()` gains `parallel = FALSE` and `ncores = NULL` arguments. With
  `parallel = TRUE`, the resampling loops of the two slow bootstrap tests --
  **Stute-Zhu** (parametric bootstrap) and **Lai-Liu-HL** (stratified
  resampling) -- run on a local PSOCK cluster via base
  `parallel::parLapply()`, which works on every platform including Windows.
  `ncores` defaults to `max(1, parallel::detectCores() - 1)`. Everything else
  in the battery is untouched, and the sequential default reproduces previous
  versions exactly.
* Parallel runs are reproducible: the workers' RNG streams are initialized
  with `parallel::clusterSetRNGStream()`, seeded deterministically from the
  session RNG, so two runs from the same `set.seed()` state (and the same
  `ncores`) give identical bootstrap p-values. The parallel L'Ecuyer-CMRG
  streams necessarily differ from the serial stream, so `parallel = TRUE`
  and `parallel = FALSE` results differ within Monte-Carlo error at the same
  seed (standard behaviour, documented in `?run.all.gof`).
* Base `parallel` added to `Imports` (no new external dependencies).
* `edges.gof()` — a brand-name alias for `def.ensemble.gof()`, matching the
  ensemble's published name **EDGES** (the Cauchy-combination ensemble of the
  EDGE directed bases). It takes the same arguments and returns the same value;
  `def.ensemble.gof()` is retained unchanged for backward compatibility.

## New data

* New bundled dataset **`gof_demo_grouped`** — a companion to `gof_demo` built
  from the same data-generating process, but with `age` recorded in 10-year
  bins and `bmi` as integers *before* the outcome is generated, so the 800
  observations share only 328 covariate patterns. It demonstrates the
  sparse-versus-grouped distinction that `run.all.gof()` reports side by side:
  on this data the per-observation ("sparse") deviance and the
  covariate-pattern ("grouped") deviance reach *opposite* verdicts at the 5%
  level, while on the all-continuous `gof_demo` the two forms coincide (the
  degenerate case). Generated reproducibly in `data-raw/make_gof_demo.R`.
* New vignette subsection "Sparse or grouped? Let the package show you both"
  walking through the demo.

# ebrahim.gof 2.3.0

## Documentation and data (toolbox reframe)

* The package is reframed as a **goodness-of-fit and calibration toolbox** for
  logistic regression (Title/Description/README updated), foregrounding the
  one-call `run.all.gof()` battery and the author's own sparse-data tests
  (EF / DEF / EDGE / ensemble), with the aggregated classical and modern tests
  credited to their authors.
* New bundled dataset **`gof_demo`** — a small synthetic dataset with a
  documented, reproducible smooth (quadratic-in-age) misfit, for illustrating
  the battery (see `data-raw/make_gof_demo.R`).
* New vignette **"A goodness-of-fit and calibration toolbox for logistic
  regression"** walking through the battery on `gof_demo`.
* Added a package-level help page (`?ebrahim.gof`).

# ebrahim.gof 2.2.0

## New features

* `edge.gof()` — the new primary interface to the directed test, matching the
  method's published name **EDGE** (Ebrahim Directed Goodness-of-fit
  Evaluation). It computes exactly the same statistic as `def.gof()` and
  returns `Test = "EDGE"` in the result row. `def.gof()` is retained,
  unchanged, as a fully supported legacy name — existing code keeps working.
* Includes the grouped-form additions developed since 2.1.0 (see 2.1.1 notes
  below): grouped (covariate-pattern) forms of Pearson / Deviance / McCullagh
  reported alongside the sparse forms, and the Osius-Rojek test computed on
  covariate patterns per its classical definition.

# ebrahim.gof 2.1.1

## New features

* `gof_install_suggests()` - installs the optional ("Suggests") packages that
  the slow tests in `run.all.gof()` rely on (`givitiR` + `callr` for GiViTI,
  `mgcv` for the GAM tests, `BAGofT`, `ResourceSelection`). It only installs the
  ones that are missing (anything already present is left untouched), and asks
  for confirmation first (in interactive sessions); pass `ask = FALSE` to skip
  the prompt in a setup script. With `update = TRUE` it also checks
  `old.packages()` and offers to update any that are out of date.
* `run.all.gof()` gains an `install` argument (`"ask"` by default). When a test
  in the run needs an optional package that is not installed, an interactive
  session offers to install it (with `"ask"`) or installs it silently (with
  `"yes"`); `"no"` keeps the previous behaviour of just skipping the test with a
  note. Per CRAN policy, a non-interactive session (scripts, `R CMD check`)
  never installs anything, regardless of the setting -- the library is never
  modified without the user's explicit consent.

# ebrahim.gof 2.1.0

## New features

* `cdef.gof()` - the Covariate-Space Directed Ebrahim-Farrington test. Like
  `def.gof()` but the directed basis lives in covariate space (polynomials and
  pairwise products, natural splines, or a `"combined"` basis that also includes
  fitted-probability bends), with the same Farrington Omega-projection calibration.
  It detects omitted interactions and local/oscillatory misfit that
  fitted-probability grouping can miss; rank-deficient bases are reduced
  automatically.
* `gof.features()` - the goodness-of-fit evidence vector (one-sided z-scores from
  a panel of tests plus the covariate-space directed tests), the input to a
  learned-ensemble GOF test.
* `deploy.gof()` - a deployable learned-ensemble test: given a pre-trained scorer,
  it calibrates the p-value by a per-dataset parametric bootstrap from the fitted
  model, so it is valid on any data set without knowing the truth.

## `run.all.gof()` additions and improvements

* `McCullagh` - the McCullagh (1985) exact-conditional-moments standardization of
  the Pearson statistic (SAS GOFLOGIT / Kuss 2002 algorithm). Verified to
  reproduce the thesis low-birth-weight result (p = 0.937) to machine precision.
* `GiViTI` - the GiViTI polynomial calibration test (Nattino, Finazzi &
  Bertolini), wrapping `givitiR` run inside an isolated `callr` subprocess so a
  crash in givitiR's compiled dependencies returns `NA` instead of aborting the
  session. Verified against the thesis result (internal p = 0.586). Opt-in slow;
  `control = list(GiViTI = list(devel = "internal"))`. Adds `givitiR` and `callr`
  to Suggests.
* `BAGofT` now runs on single-predictor models: its random-forest partitioner
  needs at least two predictors, so a constant helper column is added to the data
  (not the formula) - the workaround documented in Kuss (2002) / the thesis -
  instead of failing.
* Reworked output: `run.all.gof()` returns an object of class `gof_battery` (still
  a `data.frame`) with a dedicated `print` method - rows grouped by test family,
  p-values formatted (four decimals or scientific, `-` when not available), and a
  significance flag. All `Note` messages were rewritten to clear, human-readable
  phrases.
* `include_slow` now defaults to `TRUE`, so the full battery runs by default; a
  one-time message notes which slow tests are included and that
  `include_slow = FALSE` gives a quick fast-tests-only run.
* New `calibration_plot` argument: with `calibration_plot = TRUE` (and `GiViTI`
  among the tests) the GiViTI calibration belt is computed, stored on the result,
  and drawn; a `plot()` method (`plot.gof_battery`) redraws the stored belt.
* `F-test` - the modified Hosmer-Lemeshow F-test (deviance residuals
  ANOVA-F-tested across deciles), following `LogisticDx::gof.glm`.
* `GiViTI-external` - the GiViTI calibration test under the external development
  assumption, so the internal and external forms now run side by side (matching
  the thesis, which reported both).
* BAGofT's per-simulation console output ("Calculating results from N th ...") is
  now suppressed, so the printed battery stays clean.

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
  Stute-Zhu), the e-value HL (`eHL`; needs `isotone`), and `BAGofT`.

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