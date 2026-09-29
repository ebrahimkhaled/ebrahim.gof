# ebrahim.gof 2.9.0

## Bug fix

* `deepgof1()` refitted each bootstrap sample by evaluating the model formula on the model frame.
  For a term that transforms a covariate, such as `log(x)`, `ns(x, 3)` or `poly(x, 2)`, the model
  frame holds the transformed column and not `x`, so every refit failed, every replicate was scored
  `+Inf`, and the p-value was 1 whatever the data. The bootstrap now refits on the fitted model's
  design matrix with `glm.fit()`, which also keeps a spline basis fixed, as a parametric bootstrap
  under the fitted model requires. For models without such terms the p-value for a given seed is
  the same as in 2.8.0.

* `deepgof1()` now stops with a clear message for a grouped (`cbind(successes, failures)`) or
  weighted binomial fit. The bootstrap draws one Bernoulli outcome per row, so such fits were never
  served correctly.

## New features

* `deepgof1(reading = "allpairs")` scores the residual map of every pair of covariates and takes
  the largest score as its statistic. The same maximum is taken in every bootstrap replicate, so the
  p-value needs no correction for the choice of pair. The default rule (`reading = "v1"`, unchanged)
  reads one map, over the two columns with the largest |b| * sd; it reads covariates by their linear
  effects, so it can pass over a covariate whose effect is U-shaped. The two readings suit different
  designs. With few covariates the maximum over pairs has more power: on the benchmark of Liu et al.
  (2024, Statistics and Computing 34:175), whose settings have at most three covariates, with B = 199
  and matched level .05, it has power .640 against .580 for the default, and a null rejection rate
  of .050. With many covariates the maximum pays for the number of pairs: with two active covariates
  among ten, the default rule finds the active pair in 93 to 100 per cent of datasets and has more
  power than the all-pairs reading in all twelve settings studied (.453 against .286 on average).

* In the all-pairs reading the maps are drawn over the covariates, not over the columns of the model
  matrix: a covariate that enters as `ns(x, 3)`, `poly(x, 2)` or `I(x^2)` is ranked by `x` itself, and
  an interaction adds no axis of its own. Transformed covariates are read from the data the model was
  fitted to, on the rows the fit kept. A factor is one covariate, entered by its level codes. The
  argument `covariates` restricts the pairs to a chosen set.

* A model with one covariate is now accepted: the map is 36 quantile cells along its ranks, the
  construction of the benchmark above for its one-covariate setting, where it has power .735 at
  matched level .05. Earlier versions stopped with an error.

* The result now also holds `map`, the 6 x 6 map that gave the statistic (rows follow the first
  axis), and, for the all-pairs reading, `pairs`, the observed score of every pair, so the test says
  where the misfit lies as well as whether it is there.

# ebrahim.gof 2.8.0

## Bug fix

* `deepgof1()` put covariate values that are tied into cells of its residual map by row order.
  When the rows of the data are sorted by the outcome, as clinical files often are, tied rows were
  then placed by their outcome, while the bootstrap outcomes were not sorted, so a correct model
  could be rejected. On `aplore3::glow500`, sorted with all fractures last, the model with the
  interactions AGE x PRIORFRAC and MOMFRAC x ARMASSIST gave p = .005 in file order and a median
  p of .54 over 30 random row orders. Ties among covariate values are now broken at random. The
  random order is drawn once per call and used for the observed map and for every bootstrap map, so
  the p-value no longer depends on the order of the rows. For rows in random order the test has
  the same distribution as before. When no covariate column has a tie nothing extra is drawn, and
  the p-value for a given seed is the same as in 2.7.0.

* The `Stukel` row of `run.all.gof()` squared and summed the two marginal score statistics of
  Stukel's two tail directions and referred the sum to chi-squared on 2 degrees of freedom. That is
  the statistic of `LogisticDx::gof.glm()` ("SstBoth"), which the battery was built to reproduce,
  but it is not chi-squared on 2 degrees of freedom: once the model is fitted the two directions are
  correlated, and the sum is liberal: for example, a post-fit correlation of about -0.71 in a typical
  design, and a rejection rate of up to about 7 per cent at the 5 per cent level. Against symmetric
  departures from the logit it also loses much of its power. The row now reports the joint score
  statistic, which is the score test for adding both directions to the model and agrees with
  `anova(..., test = "Rao")` up to glm's convergence tolerance. The joint test gives up a little power
  against one-sided (cloglog-type) departures. The old statistic is still available, for reproducing
  earlier results only, with `control = list(Stukel = list(form = "marginal"))`. The 2.0.0 entry
  below, which says the row matches `LogisticDx`, now holds only for that form. Aliased columns of the
  model are left out of the joint and likelihood-ratio forms.

* When every fitted risk lay on one side of one half, one of Stukel's directions was identically
  zero and the row returned `NaN` without a note. It now reports the one-degree-of-freedom score test
  on the remaining direction, and says so in `Note`.

* `def.gof()` and `edge.gof()` with `basis = "stukel"` stopped with "system is computationally
  singular" when one group's mean fitted risk lay just above one half. The basis column for risks at
  or above one half was then about 1e-7: large enough to pass the old filter, small enough to break
  the solve. The kept columns are now scaled to unit length before any solve, which does not change
  the statistic.

## New features

* `control = list(Stukel = list(form = "lr"))` gives the likelihood-ratio test for the same two
  directions: the model is refitted with them added, and the drop in deviance is referred to
  chi-squared on the number of added columns the refit can estimate (one when every fitted risk lies
  on one side of one half). If the refit fails or does not converge, as it can under separation, the
  row returns `NA` and says why. The joint and likelihood-ratio forms need an unweighted logit fit to
  binary data, and return `NA` with a note otherwise.

* `edge.gof()` and `def.gof()` gain the basis `"sym"`: one column, eta|eta|, at the logit of each
  group's mean fitted risk. It is Stukel's symmetric direction in grouped form and is aimed at tails
  that are too heavy or too light on both sides, such as a probit or cauchit truth fitted by a logit.
  The battery reports it as `DEF.sym`, so the fast battery has one more row.
  `def.ensemble.gof()` accepts `"sym"` as a component; its default components are unchanged.

* `edge.gof()`, `def.gof()` and `def.ensemble.gof()` gain `weights = c("unit", "score")`. `"unit"`
  is the published statistic and stays the default. `"score"` multiplies each basis column by the
  square root of its group's variance. For a logit fit the statistic then becomes the Rao score test
  for adding the grouped columns to the model (for other links it is a score-type test). It is
  referred to chi-squared on the rank of its information matrix, which is the number of columns
  unless one is redundant. When group variances differ strongly, as they do at high discrimination,
  this keeps shapes on the logit scale from losing their signal. The result keeps its six columns, with `Method = "score"` and an integer `df`. In the
  battery, use `control = list(DEF.sym = list(weights = "score"))`. The two ensemble rows of the
  battery still combine the unit form of `DEF.poly2`, `DEF.poly3` and `DEF.stukel` at the battery's
  `G`; when `control` gives those rows other `weights` or another `G`, the ensemble rows say
  "unit form" in `Note`.

* `edge.gof()`, `def.gof()` and `def.ensemble.gof()` accept `G = "auto"`, which uses
  `max(10, round(n / 25))` groups, the partition rule of the EDGE paper. A number is used as before.
  In the battery, `control = list(DEF.poly3 = list(G = "auto"))` applies it to one directed row, and
  `Note` records the number of groups used. `run.all.gof(G = "auto")` resolves the rule once, so
  every row, the EF and ensemble rows included, uses the same number of groups, and the directed rows
  record it in `Note`. Any other non-numeric `G` is now an error rather than a failure of each row.

* `projection.gof()` is the projection test of Escanciano (2006), as defined for logistic
  regression by Liu et al. (2024, Statistics and Computing 34:175), also in the battery as the slow
  row `"Projection"` (family "Bootstrap"). The residual process is taken along every direction of the
  covariate space rather than only along the fitted linear predictor, as in `Stute-Zhu`, so the test
  also sees departures such as an omitted interaction. The weight is the closed form of the integral
  over the sphere, (pi - angle) / (2 pi), on the raw covariates, and the p-value comes from Liu et
  al.'s model-based bootstrap with B = 1000 refits by default, the number they use in their data
  examples. The weight costs O(n^3) time and O(n^2) memory, in pure R; the battery row is skipped
  above n = 3000 unless `control = list(Projection = list(max_n = ...))` raises it. The tests check
  the weight against a Monte Carlo over random directions, and exactly in the one-covariate case.
  Timing with two covariates and B = 1000: 1.1 s at n = 200 and 6.4 s at n = 500, of which the
  weight takes 0.2 s and 4.4 s.

* `bagoft.fast()` computes the BAGofT test of Zhang, Ding and Yang (2023, JASA) and returns the same
  p-values as the BAGofT package 1.0.0 for the same seed, identical to the last bit, including the
  distance-correlation pre-selection that BAGofT runs above five covariates. It draws every random
  number by the same call, in the same order, and removes the package's per-split overhead (formula
  parsing, model frames, `xtabs()`, `cut()` labels, the `predict()` wrappers). Timing with the
  package defaults (nsplits = 100, nsim = 100), two covariates: 99 s against 239 s for the package at
  n = 200, and 188 s at n = 500, where the package is estimated at about 360 s (39 s against 20 s
  measured with nsim = 10). The rest of the time is the forests themselves. It also runs where
  BAGofT 1.0.0 stops: a single covariate, and formulas with transformed terms. The code is adapted
  from BAGofT (GPL-3); its authors are credited in `DESCRIPTION` as contributors and copyright
  holders. It needs `randomForest`, and `dcov` above five covariates, both now in Suggests.

* The `BAGofT` row of `run.all.gof()` now uses `bagoft.fast()` by default when `randomForest` is
  installed. The p-value is the same as before for the same seed; only the time changes.
  `control = list(BAGofT = list(engine = "package"))` calls the BAGofT package as before, and `Note`
  says which engine ran.

## Behaviour changes

* `def.gof()`, and so `edge.gof()` and `def.ensemble.gof()`, now warn when there are fewer events,
  or fewer non-events, than groups. The p-value is still returned, but some groups then hold almost
  no events and the grouped reference distribution is unreliable. `def.ensemble.gof()` warns once
  rather than once per basis, and in `run.all.gof()` the directed rows report it in `Note` instead.

* A sample with no event, or no non-event, has no maximum-likelihood fit. `def.gof()` in either form,
  `edge.gof()` and `def.ensemble.gof()` now return `NA` for it, with a warning of class
  `def_degenerate`, instead of stopping with an error or, for the score form, returning a p-value of
  zero. `def.ensemble.gof()` warns once. In `run.all.gof()` the directed rows say so in `Note`, and
  every form of the `Stukel` row, `"marginal"` included, returns `NA` with a note.

* The score form of `def.gof()` (and `edge.gof()`) and the joint form of the `Stukel` row now leave
  out a column that the fitted model already spans: one whose information after the fit is below
  1e-10 of its information before the fit. This happens when the fitted logit is constant or nearly
  so, for example in a sample with two or three events. The information of such a column is rounding
  noise, and inverting it could give a p-value of zero. When no column is left, `def.gof()` returns
  `NA` with a warning of class `def_no_information`, and the `Stukel` row returns `NA` with a note.

* `deepgof1()`: when some covariate column has tied values, the random tie-breaking draws from the
  random-number stream, so for a given seed such calls give a different p-value than in 2.7.0.
  Calls without ties give the same p-value as before. The map still uses the two covariates selected
  by the fitted coefficients, as in 2.7.0.

* The full battery of `run.all.gof()` (`include_slow = TRUE`) has one more slow row, `Projection`.
  The `BAGofT` row no longer needs the BAGofT package when `randomForest` is installed, and
  `gof_install_suggests()` now offers `randomForest` and `dcov` for it.

## Documentation

* `?run.all.gof` said that Stukel's two-parameter form does not always hold its nominal level. That
  was a property of the summed statistic, not of the test, and the sentence has been replaced.

* The reference to Zhang, Ding and Yang in `?run.all.gof` gave the wrong issue and pages; it is now
  JASA 118(542), 1115-1125 (2023), and a second, incomplete entry for the same paper has been
  removed.

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

* The catalogue entry for `Deviance` called it "the more conservative of the two" against
  `Pearson`. Measured on correctly specified models it rejects 100 per cent where Pearson
  rejects none: on ungrouped binary data its expectation depends on the fitted risks alone,
  so the direction of the error is set by the risk profile and not by the fit. Corrected.

* `?deepgof1` said to see `run.all.gof()` to run it alongside the battery. It is not one of
  the tests the battery selects, so naming it in `tests` errors. Both that sentence and the
  battery's own text now say which procedures sit outside the panel.

* `?shrink.gof` gained a section on choosing between it and `calm.gof()`, including the
  trap that `calm.gof()` takes `lambda_scale` and `shrink.gof()` does not, so the same
  number means penalties a factor of n apart.

* The Details headings of `?run.all.gof` now carry the literal `Family` labels that appear
  on every row of the output, so a reader can go from a label to its section.

* The examples on `?def.gof` and `?def.ensemble.gof` fitted correctly specified models, so
  no test ever fired. They now use `gof_demo`, whose documented misfit the directed tests
  detect (p = 0.015), and then the corrected model, which they do not (p = 0.797).

* Both vignettes predated 2.5.0 and mentioned none of the methods added since. The toolbox
  vignette gained a section on the penalised, frozen-weight and pretrained procedures.

* A pkgdown site is configured (`_pkgdown.yml` and a workflow), giving the reference pages
  and vignettes indexable HTML addresses.

* `ef.gof()` now warns, rather than messages, when `model` or `m` is supplied while `G` is
  left at its default, since both are then ignored.

* Smoke tests were added for ten exports that had none, `calm.gof()` among them.

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