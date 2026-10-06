#' Run a Battery of Goodness-of-Fit Tests at Once
#'
#' @description
#' Runs several goodness-of-fit tests for a binary logistic regression in one
#' call and returns one tidy \code{data.frame}, one row per test. Pass a fitted
#' \code{glm} to run the whole battery; pass \code{(y, predicted_probs)} to run
#' the tests that need only predictions. Each test is wrapped so that a failure of
#' one test never aborts the whole run.
#'
#' @details
#' \strong{How to read the battery.} Every test here answers the same question --
#' does the fitted model describe the data it was fitted to -- but they differ in
#' the departure each is built to notice. That is why the panel is more
#' informative than any single p-value: agreement across families is evidence of
#' fit, and disagreement tells you \emph{what kind} of misfit is present. A test
#' that rejects points at the departure its own construction is sensitive to.
#'
#' The one thing the panel cannot do is rescue an invalid reference
#' distribution. Under sparse data -- almost every covariate pattern unique --
#' the classical chi-square references fail, which is the situation this package
#' was written for.
#'
#' Every row the battery returns carries a \code{Family} label, and the sections below are
#' those labels: find the label in the output, then the section of the same name here.
#'
#' \strong{Global and Standardized statistics} (\code{Family} "Global", "Standardized"). These compare observed and
#' fitted responses over the whole sample without grouping.
#' \itemize{
#'   \item \code{Pearson} -- the sum of squared Pearson residuals. Its
#'     chi-square reference assumes many observations per covariate pattern. On
#'     sparse data that assumption fails and the test is unreliable in both
#'     directions; it is reported for completeness and comparison, not for use.
#'     It is printed without a significance flag and is not counted among the
#'     rejections of the header or of \code{summary()}; the \code{Pearson (grouped)}
#'     row, shown when covariate patterns repeat, is the one to read.
#'   \item \code{Deviance} -- the likelihood-ratio statistic against the
#'     saturated model. On ungrouped binary data its expectation depends on the
#'     fitted risks alone and not on whether they agree with the responses, so it
#'     sits above \eqn{n - p} when the fitted risks are near one half and far below
#'     it when they are extreme; the direction of the error is set by the risk
#'     profile rather than by the fit. On a correctly specified model with
#'     mid-range risks it therefore rejects almost always. It is reported for
#'     completeness and comparison, not for use, and like \code{Pearson} it is not
#'     counted among the rejections.
#'   \item \code{Osius-Rojek} -- rescues the Pearson statistic by standardizing
#'     it with its asymptotic mean and variance computed under increasing sample
#'     size rather than increasing cell counts, giving a normal reference that
#'     stays valid when patterns are unique; the deviate is referred two-sided.
#'     Fast, parameter-free, and a
#'     sensible default. It can be anti-conservative at small samples and
#'     conservative under a strongly skewed covariate.
#'   \item \code{McCullagh} -- standardizes the Pearson statistic by its
#'     \emph{exact} conditional moments rather than asymptotic ones (computed by
#'     the Kuss 2002 algorithm). Typically the most powerful of the global
#'     statistics. The deviate is referred to the upper tail of the normal, as
#'     in the GOFLOGIT macro, so a Pearson statistic below its conditional mean
#'     is not counted as misfit. The moments are computed from the p x p
#'     information matrix, so the cost is linear in the sample size.
#'   \item \code{Copas-RSS} -- the unweighted sum of squares of the raw
#'     residuals, standardized by its own moments and referred two-sided to the
#'     normal (\code{df} is \code{NA}, as for the other normal deviates).
#'     Weighting each residual equally rather than by its variance makes it
#'     comparatively sensitive to misfit in the middle of the risk range.
#'   \item \code{Information-Matrix} -- the White/Orme test. It compares two
#'     estimators of the information matrix that agree only if the model is
#'     correctly specified, so it is an omnibus check on the whole
#'     specification rather than on calibration alone.
#' }
#'
#' \strong{Partition tests} (\code{Family} "Partition"). These sort observations by fitted risk, group them,
#' and compare observed with expected counts group by group.
#' \itemize{
#'   \item \code{HL} -- the Hosmer-Lemeshow test, the field's default: G groups
#'     cut at percentiles of fitted risk (the "deciles of risk" when
#'     \code{G = 10}), referred to chi-square on \eqn{G - 2} degrees of freedom.
#'     Note that this reference was established by simulation, not derivation.
#'     Familiar and cheap, but modest in power, and its result depends on the
#'     grouping.
#'   \item \code{HL-equalwidth} -- the same statistic with groups cut at fixed
#'     probability intervals instead of percentiles. When fitted risks are
#'     concentrated in a narrow range, intervals can come out empty and the test
#'     may not be computable at all, which is why percentile grouping is usually
#'     preferred.
#'   \item \code{Pigeon-Heyse} -- a variance-corrected Hosmer-Lemeshow statistic
#'     that accounts for the variability of fitted probabilities within each
#'     group. The correction is conservative in sparse designs, so a
#'     non-rejection carries less weight than the nominal level suggests.
#'   \item \code{F-test} -- deviance residuals compared across the risk groups by
#'     a one-way analysis of variance. It is markedly liberal under sparsity and
#'     grows more so with sample size; it is included for comparison, printed
#'     without a significance flag and not counted among the rejections.
#'   \item \code{EF}, \code{EF-normal} -- the Ebrahim-Farrington test,
#'     Farrington's correction on the Hosmer-Lemeshow risk groups. \code{EF} uses
#'     the chi-square reference on \eqn{G - 2} degrees of freedom at the battery's
#'     \code{G}; \code{EF-normal} uses the plug-in normal reference, which is the
#'     reference for many groups, at \code{G = "auto"} (groups of about 25 records)
#'     unless \code{control = list("EF-normal" = list(G = ...))} sets another
#'     number; it needs the model or \code{X}. See \code{\link{ef.gof}}. Tied fitted
#'     risks are put in random order.
#'   \item \code{Lai-Liu-HL} -- Lai and Liu's procedure for using the
#'     Hosmer-Lemeshow test in large samples, where any test eventually rejects.
#'     It standardizes power to a reference sample size \code{n0} and so returns
#'     no p-value: the statistic is the standardized power and the accept/reject
#'     decision appears in the \code{Note} column. Tune with
#'     \code{control = list("Lai-Liu-HL" = list(n0 = ..., k = ...))}.
#'   \item \code{HL-largeN} -- Nattino, Pennell and Lemeshow's (2020) large-sample
#'     Hosmer-Lemeshow test. Instead of perfect fit it tests whether the misfit,
#'     measured by \eqn{\epsilon = \sqrt{\lambda/n}}, exceeds a tolerance \eqn{\epsilon_0}:
#'     the ordinary statistic is referred to a noncentral \eqn{\chi^2_{G-2}} with
#'     noncentrality \eqn{\epsilon_0^2 n}. By their convention \eqn{\epsilon_0} is the
#'     misfit that would be just significant at \code{n0 = 1e6}; change it with
#'     \code{control = list("HL-largeN" = list(n0 = ...))}. The \code{Note} gives
#'     \eqn{\epsilon_0} and the estimate \eqn{\hat\epsilon}. Meant for samples in the tens of
#'     thousands and above; in small samples it is conservative.
#' }
#'
#' \strong{Directed tests} (\code{Family} "Directed"). Rather than asking whether anything is wrong, these
#' ask whether a \emph{particular} shape of departure is present, which buys
#' power when the guess is right.
#' \itemize{
#'   \item \code{EDGE}, \code{EDGE.G10} -- the EDGE test (the cubic basis of
#'     \code{\link{edge.gof}}) at its default partition, \code{G = "auto"}, and at
#'     ten groups, whatever the battery's \code{G}. When the default partition has
#'     ten groups (\eqn{n \le 250}) the two are one test and only \code{EDGE} is
#'     shown. The statistic \eqn{S} is a weighted sum of chi-squares: the note
#'     says when its p-value comes from \eqn{S/c \sim \chi^2_{df}}
#'     (Satterthwaite), so \code{Statistic} and \code{df} alone do not give it.
#'   \item \code{DEF.poly2}, \code{DEF.poly3}, \code{DEF.stukel}, \code{DEF.sym} --
#'     the directed forms, each aiming the test at a smooth departure in the shape
#'     of the calibration curve: a quadratic or cubic drift in the linear
#'     predictor, Stukel's asymmetry-and-tail family, or Stukel's symmetric
#'     direction, tails too heavy or too light on both sides. Powerful when the
#'     misfit resembles the chosen basis, weaker when it does not. Each row takes
#'     \code{weights} and \code{G} through \code{control}, for example
#'     \code{control = list(DEF.sym = list(weights = "score", G = "auto"))}; see
#'     \code{\link{def.gof}}. In the full battery \code{DEF.poly3} is left out
#'     when it would repeat \code{EDGE.G10} (the battery's \code{G} is 10 and
#'     \code{control} does not move either row); name it in \code{tests} to get it. On a sample with no event, or no non-event, there is
#'     no fitted model; these rows and the \code{Stukel} row are then \code{NA}, and
#'     \code{Note} says why.
#'   \item \code{Stukel} -- a score test against Stukel's generalized logistic
#'     link, which nests the logit and lets the two tails bend independently. It
#'     is aimed squarely at link misspecification. The two tail directions are
#'     tested jointly on 2 degrees of freedom (1 when every fitted risk lies on
#'     one side of one half, which \code{Note} then says). Testing them jointly
#'     gives up a little power against one-sided (cloglog-type) departures. Up to
#'     2.7.0 this row summed two marginal statistics and was liberal; see NEWS.
#'     \code{control = list(Stukel = list(form = "lr"))} gives the likelihood-ratio
#'     test for the same two directions, and \code{form = "marginal"} the old sum,
#'     for reproducing earlier results only. The likelihood-ratio refit can fail
#'     to converge under separation; the row is then \code{NA}, with a note. The
#'     joint form leaves out a direction whose information after the fit is below
#'     \eqn{10^{-10}} times its information before the fit, as when the fitted logit
#'     is constant; when no direction is left the row is \code{NA}, with a note.
#' }
#'
#' \strong{Covariate-space tests} (\code{Family} "Covariate-space"). These partition the covariates themselves
#' rather than the fitted risk, so they can see structure that risk-ordering
#' averages away -- an omitted interaction, for instance, need not disturb the
#' marginal calibration at all.
#' \itemize{
#'   \item \code{Tsiatis} -- a score test that adds indicator variables for
#'     regions of the covariate space and asks whether they improve the fit.
#'   \item \code{Xie} -- clusters the covariate space and compares observed with
#'     expected counts per cluster, using the corrected degrees of freedom
#'     \eqn{G - k/2 - 1} with \eqn{k} the number of predictors.
#'   \item \code{Pulkstenis-Robinson} -- crosses a categorical covariate with
#'     risk groups, so it needs at least one categorical predictor. The package
#'     auto-detects one (any factor, character or logical, or a numeric with few
#'     distinct values, controlled by
#'     \code{getOption("ebrahim.gof.pr.maxlev", 6)}) and returns \code{NA} with a
#'     note when there is none.
#' }
#'
#' \strong{Smoothing and GAM tests} (\code{Family} "Smoothing", "GAM"). These replace grouping with a smoother, so nothing
#' is lost to an arbitrary choice of bin edges.
#' \itemize{
#'   \item \code{le-Cessie} -- the le Cessie-van Houwelingen score test, which
#'     smooths the residuals over the covariate space with a kernel and asks
#'     whether the smoothed surface is further from zero than chance allows.
#'     Sensitive to local structure that omnibus statistics average away, and one
#'     of the most expensive tests here, with cost growing sharply in the sample
#'     size.
#'   \item \code{HL-GAM}, \code{PR-GAM}, \code{Xie-GAM} -- variants that fit a
#'     deliberately overfitted generalized additive model and use it to define
#'     the grouping, letting the data rather than the analyst choose where the
#'     boundaries fall. These need \pkg{mgcv}.
#' }
#'
#' \strong{Resampling tests} (\code{Family} "Bootstrap"). When a statistic has no usable closed-form
#' reference, these build one by simulation.
#' \itemize{
#'   \item \code{Stute-Zhu} -- a cumulative-residual test: residuals are
#'     accumulated along the fitted linear predictor and the largest excursion of
#'     that path is compared with a parametric bootstrap. It needs no binning and
#'     no bandwidth, and in practice it is the best-behaved test in the battery
#'     on size. The p-value is \eqn{(1 + k)/(1 + B)}, with \eqn{k} the number of
#'     bootstrap statistics at least as large as the observed one, so it is never
#'     0. Set the number of resamples with
#'     \code{control = list("Stute-Zhu" = list(B = ...))}.
#'   \item \code{BAGofT} -- the binary adaptive test, which splits the data, uses
#'     one part to learn a partition that separates fitted from observed, and
#'     tests on the other. Its behaviour depends strongly on how many splits and
#'     resamples it is given; set them with
#'     \code{control = list(BAGofT = list(nsim = ...))} and be aware that the
#'     published default is far more expensive than a single split. By default
#'     the row is computed by \code{\link{bagoft.fast}}, which returns the same
#'     p-value as the \pkg{BAGofT} package for the same seed and needs only
#'     \pkg{randomForest} (and \pkg{dcov} above five covariates);
#'     \code{control = list(BAGofT = list(engine = "package"))} calls
#'     \pkg{BAGofT} itself. When no simulated statistic is as extreme as the
#'     observed one, \pkg{BAGofT} returns 0; the row then reports
#'     \eqn{1/(nsim + 1)}, the smallest Monte Carlo p-value, and says so.
#'   \item \code{Projection} -- the projection test of Escanciano (2006) as
#'     defined for logistic regression by Liu et al. (2024): the cumulative
#'     residual process is taken along every direction of the covariate space,
#'     not only along the fitted linear predictor as in \code{Stute-Zhu}, so it
#'     also sees departures such as an omitted interaction. Model-based
#'     bootstrap with \code{B = 1000} refits by default; see
#'     \code{\link{projection.gof}}. Its weight matrix costs \eqn{O(n^3)} time and
#'     \eqn{O(n^2)} memory, so the row is skipped above \eqn{n = 3000}. Set
#'     \code{control = list(Projection = list(B = ..., max_n = ...))}.
#' }
#'
#' \strong{Calibration tests} (\code{Family} "Calibration"). These come from clinical prediction, and ask
#' directly whether predicted risks match observed frequencies.
#' \itemize{
#'   \item \code{GiViTI}, \code{GiViTI-external} -- the GiViTI polynomial
#'     calibration test, which fits a polynomial of the fitted risk and tests
#'     whether it departs from the identity, under the internal and external
#'     development assumptions respectively. The two are \emph{not}
#'     interchangeable: on data used to fit the model the internal form is the
#'     appropriate one, so the full battery runs \code{GiViTI} only, and
#'     \code{GiViTI-external} runs when it is named in \code{tests}. For frozen
#'     predictions on new data use \code{\link{run.all.external}}. It also produces the calibration belt, which shows
#'     \emph{where} on the risk scale a model drifts; see
#'     \code{calibration_plot}. Wraps \pkg{givitiR} in an isolated \pkg{callr}
#'     subprocess, so a failure inside its compiled dependencies returns
#'     \code{NA} instead of ending your session. Select with
#'     \code{control = list(GiViTI = list(devel = "internal"))}.
#'   \item \code{eHL} -- an e-value form of the Hosmer-Lemeshow test, reported as
#'     \eqn{p = \min(1, 1/e)}. E-values are safe under optional stopping, which
#'     conventional p-values are not, but the conversion shown here is
#'     conservative.
#'   \item \code{Cubic-LR} -- the cubic calibration likelihood-ratio test: the
#'     square and cube of the logit of the fitted risk are added to a logistic
#'     regression of the outcome on it, and the drop in deviance is referred to
#'     chi-squared on 2 degrees of freedom. It is the goodness-of-link idea of
#'     Pregibon (1980) with the cubic logistic model of Morgan (1985), and a
#'     fixed-degree relative of the GiViTI polynomial; fast, and needs only the
#'     predictions. See \code{\link{cubic.calib.gof}}, which documents its origin.
#' }
#'
#' \strong{Combinations} (\code{Family} "Ensemble"). Rather than choosing one test, these pool several.
#' \itemize{
#'   \item \code{Ensemble.Vote(3DEF)} and \code{Ensemble.Univ(3DEF+EF)} -- Cauchy
#'     combinations of the directed tests, and of those plus the omnibus EF. The
#'     Cauchy combination is valid without knowing how the members correlate,
#'     which is what makes pooling dependent tests possible at all. The point is
#'     to avoid having to guess the departure in advance, at the cost of being
#'     slightly less powerful than the single best member would have been. Both
#'     rows always combine the unit form of \code{DEF.poly2}, \code{DEF.poly3} and
#'     \code{DEF.stukel} at the battery's \code{G}. When \code{control} gives those
#'     rows other \code{weights} or another \code{G}, the ensemble rows do not
#'     follow, and their \code{Note} says "unit form".
#'   \item See also \code{\link{legoft}}, a pretrained combination whose weights
#'     are fixed offline and ship frozen, so two analysts running it on the same
#'     data obtain the same p-value.
#' }
#'
#' \strong{Implementation notes.} \code{Tsiatis} and \code{Xie} cluster the
#' covariate space with k-means using a fixed internal seed, so results are
#' reproducible and your own random stream is left untouched. The equal-frequency
#' groups of \code{HL}, \code{F-test}, \code{EF} and the \code{DEF} rows split tied
#' fitted risks by row order, so with many ties (grouped data, or a model on
#' discrete covariates) their results can depend on the order of the rows;
#' randomising the row order is advised. Every bundled test
#' reproduces the implementation used in the original simulation study:
#' \code{Osius-Rojek} follows \pkg{LogisticDx}'s \code{gof.glm},
#' \code{Copas-RSS} follows the \pkg{rms} gof residual, and
#' \code{HL} follows \code{ResourceSelection::hoslem.test}. The exception is
#' \code{Stukel}: its default joint score statistic agrees with
#' \code{anova(..., test = "Rao")} for the augmented model, up to glm's
#' convergence tolerance, and only
#' \code{form = "marginal"} reproduces \pkg{LogisticDx} (through
#' \code{statmod::glm.scoretest} when \pkg{statmod} is installed).
#'
#' \strong{Procedures the battery does not select.} Some of the package's own methods are
#' not part of the panel and are called directly on the fitted model, their p-values read
#' beside it: \code{\link{deepgof1}} and \code{\link{legoft}} are of that kind. Naming
#' them in \code{tests} does not reach them.
#'
#' For goodness of fit after \emph{penalized} fitting, where none of the above references
#' are valid because the coefficients are shrunk, use \code{\link{calm.gof}} or
#' \code{\link{shrink.gof}}. These take the design, the response and the penalty rather
#' than a fitted \code{glm}, so they are called separately from the battery.
#'
#' @param object A fitted binary logistic \code{\link[stats]{glm}}, with one row per
#'   binary record and no prior weights, or a binary (0/1) response vector \code{y}
#'   (then supply \code{predicted_probs}). A grouped (\code{cbind}) or weighted fit is
#'   refused: expand it to one row per record. Another link gives a warning, since
#'   most tests assume the logit. A warning is also given when some fitted risks are 0
#'   or 1 (separation) or when there are fewer than 10 events or non-events.
#' @param predicted_probs Numeric predicted probabilities; required when
#'   \code{object} is a \code{y} vector.
#' @param X Optional design matrix; lets the directed (DEF) and covariate-space tests
#'   run from the \code{(y, predicted_probs)} form. It should be the model matrix of
#'   the fit that gave \code{predicted_probs}: a data frame is coded with
#'   \code{model.matrix(~ ., X)}, and a matrix without a constant column gets an
#'   intercept column in front.
#' @param tests Either \code{"all"} (default) or a character vector of test names
#'   to run (e.g. \code{c("EF","DEF.poly3","HL")}).
#' @param G Integer number of groups passed to the grouping tests (default 10), or
#'   \code{"auto"} for \code{max(10, ceiling(n / 25))} (see \code{\link{def.gof}}).
#'   \code{"auto"} is resolved once, so every row, the ensemble rows included, uses
#'   the same number of groups; the directed rows record it in \code{Note}.
#' @param include_slow Logical; when \code{TRUE} (the default) the full battery
#'   runs, including the slow tests: le Cessie-van Houwelingen smoothing
#'   (O(n^2)-O(n^3)), the GAM tests, Stute-Zhu, eHL, BAGofT, Projection, and GiViTI. Set
#'   \code{FALSE} for a quick run with the fast tests only. A one-time message
#'   notes this whenever slow tests are included.
#' @param parallel Logical; when \code{TRUE}, the resampling loops of the slow
#'   bootstrap tests (\code{Stute-Zhu} and \code{Lai-Liu-HL}) are run on a local
#'   PSOCK cluster via \code{\link[parallel]{parLapply}} (works on all
#'   platforms, including Windows). All other tests are unaffected. The default
#'   \code{FALSE} keeps every loop sequential, exactly as in previous versions.
#' @param ncores Integer; the number of worker processes used when
#'   \code{parallel = TRUE}. The default \code{NULL} uses
#'   \code{max(1, parallel::detectCores() - 1)}. Values below 2 fall back to
#'   the sequential path.
#' @param calibration_plot Logical; when \code{TRUE} and \code{GiViTI} is among
#'   the tests, also compute and draw the GiViTI calibration belt and store it on
#'   the result (retrievable with \code{plot()}). Default \code{FALSE}.
#' @param install One of \code{"ask"} (default), \code{"no"}, or \code{"yes"},
#'   controlling what happens when a test in the run needs an optional package
#'   that is not installed. In an \emph{interactive} session, \code{"ask"} lists
#'   the missing packages and asks before installing, and \code{"yes"} installs
#'   them without asking; \code{"no"} never installs (the test is just skipped
#'   with a note). In a non-interactive session (scripts, \code{R CMD check})
#'   nothing is ever installed, regardless of this setting. See
#'   \code{\link{gof_install_suggests}}.
#' @param control Optional named list of per-test options. Recognized entries:
#'   \code{"Stute-Zhu" = list(B = ...)} (bootstrap replicates);
#'   \code{GiViTI = list(devel = "internal"/"external")};
#'   \code{"Lai-Liu-HL" = list(n0 = ..., k = ..., alpha = ...)};
#'   \code{Stukel = list(form = "joint"/"lr"/"marginal")} (the joint score test by
#'   default, the likelihood-ratio refit, or the pre-2.8.0 marginal sum);
#'   \code{DEF.poly2}, \code{DEF.poly3}, \code{DEF.stukel} and \code{DEF.sym}
#'   \code{= list(weights = "unit"/"score", G = ...)}, where \code{G} may be
#'   \code{"auto"} (see \code{\link{def.gof}}); and
#'   \code{BAGofT = list(...)} which forwards to the binary adaptive test --
#'   \code{nsim} (resampling iterations; default 100), \code{nsplits}, \code{ne}
#'   (the estimation-split size), and the random-forest partitioner's tuning
#'   \code{Kmax} (maximum number of adaptive partition cells), \code{ntree},
#'   \code{nmin}, \code{mtry}, \code{maxnodes}, and \code{engine}
#'   (\code{"auto"}, \code{"fast"} or \code{"package"}; see \code{\link{bagoft.fast}}).
#'   Example: \code{list(BAGofT = list(nsim = 200, Kmax = 8, ntree = 500))}; and
#'   \code{Projection = list(B = ..., scale = FALSE, max_n = 3000)} (bootstrap refits,
#'   standardized covariates, and the sample size above which the row is skipped; see
#'   \code{\link{projection.gof}}).
#'
#' @return A \code{data.frame} (of class \code{gof_battery}) with columns
#'   \code{Test}, \code{Family}, \code{Statistic}, \code{df}, \code{p_value},
#'   and \code{Note}, one row per test. A dedicated \code{print} method shows the
#'   rows grouped by family with formatted p-values and significance flags; the
#'   underlying columns remain available for programmatic use.
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' x <- runif(n, -3, 3)
#' y <- rbinom(n, 1, 1 / (1 + exp(-(0.6 * x))))
#' fit <- glm(y ~ x, family = binomial())
#'
#' ## The fast tests. Every covariate pattern here is unique, so this is the
#' ## sparse case the package is written for.
#' res <- run.all.gof(fit, include_slow = FALSE)
#' res
#'
#' ## The return value is a plain data.frame, so the panel can be read
#' ## programmatically as well as printed.
#' res[res$p_value < 0.05, c("Test", "Family", "p_value")]   # what rejected
#' table(res$Family)                                        # coverage by family
#'
#' ## A correctly specified model: the panel should mostly agree, and any
#' ## isolated rejection is the false positive you expect at the 5 percent level.
#' mean(res$p_value < 0.05, na.rm = TRUE)
#'
#' ## Now a model that is genuinely wrong -- the quadratic term is omitted.
#' y2  <- rbinom(n, 1, 1 / (1 + exp(-(0.6 * x + 0.5 * x^2))))
#' bad <- glm(y2 ~ x, family = binomial())
#' run.all.gof(bad, include_slow = FALSE)
#'
#' ## Pick specific tests, for instance one per family, which is the pairing the
#' ## package recommends over relying on any single statistic.
#' run.all.gof(fit, tests = c("McCullagh", "HL", "Stukel", "Tsiatis"))
#'
#' ## Grouping tests take the number of groups; the choice is a convention
#' ## rather than a derived optimum, so it is worth varying.
#' for (g in c(5, 10, 20))
#'   print(run.all.gof(fit, tests = "HL", G = g))
#'
#' \donttest{
#' ## The full battery (include_slow = TRUE by default). The slow tests need the
#' ## suggested packages mgcv, randomForest, givitiR and callr; in an interactive
#' ## session run.all.gof() offers to install any that are missing
#' ## (install = "ask"). See also gof_install_suggests().
#' ## The control= list forwards options to the individual tests; the reductions
#' ## here keep the example quick without changing what it demonstrates.
#' run.all.gof(fit, install = "no",
#'             control = list("Stute-Zhu" = list(B = 50),
#'                            BAGofT = list(nsim = 20),
#'                            Projection = list(B = 99)))
#'
#' ## The GiViTI calibration belt shows WHERE on the risk scale a model drifts,
#' ## which a single p-value cannot.
#' res2 <- run.all.gof(fit, tests = c("McCullagh", "GiViTI"),
#'                     calibration_plot = TRUE)
#' plot(res2)   # redraw the stored belt
#'
#' ## Run the bootstrap loops on a PSOCK cluster. Seeds are handled internally,
#' ## so a parallel run reproduces a serial one.
#' set.seed(1)
#' run.all.gof(fit, tests = "Stute-Zhu", parallel = TRUE, ncores = 2,
#'             control = list("Stute-Zhu" = list(B = 50)))
#' }
#'
#' @note \strong{Grouped vs sparse forms.} \code{Pearson}, \code{Deviance} and
#'   \code{McCullagh} are reported in two forms: the default (sparse / one-trial)
#'   form and a \code{"(grouped)"} form computed on the distinct covariate
#'   patterns (each a Binomial\eqn{(m_g, P_g)}). The two are identical when every
#'   covariate pattern is unique (fully sparse data, as in the simulation) and
#'   differ only when patterns repeat (\eqn{m_g > 1}). To avoid clutter, the
#'   \code{"(grouped)"} row is shown \emph{only} when it actually differs from the
#'   sparse form (i.e., when some pattern repeats); on fully sparse data it is a
#'   duplicate and is omitted. \code{Osius-Rojek} is always computed on covariate
#'   patterns, matching its classical (LogisticDx) definition.
#'
#'   \strong{Farrington vs EF.} The \emph{original} Farrington (1996) test is a
#'   grouped (covariate-pattern) test. The Ebrahim-Farrington (\code{EF}) test is
#'   its \emph{sparse-data} counterpart: it does not group by covariate pattern
#'   but forms \code{G} data-dependent bins of the predicted risk, so it applies
#'   directly to fully sparse data. Use \code{EF} for sparse binary data; the
#'   grouped Farrington form is appropriate only when covariate patterns repeat.
#'
#'   \strong{Reproducibility of the parallel path.} With \code{parallel = TRUE}
#'   the cluster's random-number streams are initialized with
#'   \code{\link[parallel]{clusterSetRNGStream}}, seeded deterministically from
#'   the session's current RNG state. Two runs from the same
#'   \code{\link{set.seed}} state (and the same \code{ncores}) therefore give
#'   identical bootstrap p-values. Note that the parallel L'Ecuyer-CMRG streams
#'   necessarily differ from the serial RNG stream, so \code{parallel = TRUE}
#'   results differ (within Monte-Carlo error) from \code{parallel = FALSE}
#'   results at the same seed; this is standard and both are valid. Results also
#'   depend on \code{ncores}, because the replicates are split across workers.
#' @seealso \code{\link{ef.gof}}, \code{\link{def.gof}}, \code{\link{def.ensemble.gof}},
#'   \code{\link{gof_install_suggests}}.
#'
#' @references
#' The aggregated tests are due to their original authors; they are provided
#' here for comparison and credited as follows.
#'
#' Farrington CP (1996). "On Assessing Goodness of Fit of Generalized Linear
#' Models to Sparse Data." \emph{Journal of the Royal Statistical Society B},
#' \strong{58}(2), 349--360. \doi{10.1111/j.2517-6161.1996.tb02086.x}
#'
#' Hosmer DW, Lemeshow S (1980). "Goodness of Fit Tests for the Multiple
#' Logistic Regression Model." \emph{Communications in Statistics -- Theory and
#' Methods}, \strong{9}(10), 1043--1069. \doi{10.1080/03610928008827941}
#'
#' McCullagh P (1985). "On the Asymptotic Distribution of Pearson's Statistic in
#' Linear Exponential Family Models." \emph{International Statistical Review},
#' \strong{53}(1), 61--67. \doi{10.2307/1402880}
#'
#' Osius G, Rojek D (1992). "Normal Goodness-of-Fit Tests for Multinomial Models
#' with Large Degrees of Freedom." \emph{Journal of the American Statistical
#' Association}, \strong{87}(420), 1145--1152.
#' \doi{10.1080/01621459.1992.10476271}
#'
#' le Cessie S, van Houwelingen JC (1991). "A Goodness-of-Fit Test for Binary
#' Regression Models, Based on Smoothing Methods." \emph{Biometrics},
#' \strong{47}(4), 1267--1282. \doi{10.2307/2532385}
#'
#' Stukel TA (1988). "Generalized Logistic Models." \emph{Journal of the
#' American Statistical Association}, \strong{83}(402), 426--431.
#' \doi{10.1080/01621459.1988.10478613}
#'
#' Stute W, Zhu LX (2002). "Model Checks for Generalized Linear Models."
#' \emph{Scandinavian Journal of Statistics}, \strong{29}(3), 535--545.
#' \doi{10.1111/1467-9469.00304}
#'
#' Escanciano JC (2006). "A Consistent Diagnostic Test for Regression Models
#' Using Projections." \emph{Econometric Theory}, \strong{22}(6), 1030--1051.
#' \doi{10.1017/S0266466606060506}
#'
#' Liu H, Li X, Chen F, Haerdle W, Liang H (2024). "A Comprehensive Comparison of
#' Goodness-of-Fit Tests for Logistic Regression Models." \emph{Statistics and
#' Computing}, \strong{34}, 175. \doi{10.1007/s11222-024-10487-5}
#'
#' Tsiatis AA (1980). "A Note on a Goodness-of-Fit Test for the Logistic
#' Regression Model." \emph{Biometrika}, \strong{67}(1), 250--251.
#' \doi{10.1093/biomet/67.1.250}
#'
#' Xie XJ, Pendergast J, Clarke W (2008). "Increasing the Power: A Practical
#' Approach to Goodness-of-Fit Test for Logistic Regression Models with
#' Continuous Predictors." \emph{Computational Statistics & Data Analysis},
#' \strong{52}(5), 2703--2713. \doi{10.1016/j.csda.2007.09.027}
#'
#' Pulkstenis E, Robinson TJ (2002). "Two Goodness-of-Fit Tests for Logistic
#' Regression Models with Continuous Covariates." \emph{Statistics in Medicine},
#' \strong{21}(1), 79--93. \doi{10.1002/sim.943}
#'
#' Nattino G, Finazzi S, Bertolini G (2014). "A New Calibration Test and a
#' Reappraisal of the Calibration Belt for the Assessment of Prediction Models
#' Based on Dichotomous Outcomes." \emph{Statistics in Medicine}, \strong{33}(14),
#' 2390--2407. \doi{10.1002/sim.6100}
#'
#'
#' Pigeon JG, Heyse JF (1999). "An Improved Goodness of Fit Statistic for
#' Probability Prediction Models." \emph{Biometrical Journal}, \strong{41}(1),
#' 71--82. \doi{10.1002/(SICI)1521-4036(199903)41:1<71::AID-BIMJ71>3.0.CO;2-O}
#'
#' Copas JB (1989). "Unweighted Sum of Squares Test for Proportions."
#' \emph{Journal of the Royal Statistical Society C}, \strong{38}(1), 71--80.
#' \doi{10.2307/2347682}
#'
#' White H (1982). "Maximum Likelihood Estimation of Misspecified Models."
#' \emph{Econometrica}, \strong{50}(1), 1--25. \doi{10.2307/1912526}
#'
#' Orme C (1988). "The Calculation of the Information Matrix Test for Binary
#' Data Models." \emph{The Manchester School}, \strong{56}(4), 370--376.
#' \doi{10.1111/j.1467-9957.1988.tb01339.x}
#'
#' Kuss O (2002). "Global Goodness-of-Fit Tests in Logistic Regression with
#' Sparse Data." \emph{Statistics in Medicine}, \strong{21}(24), 3789--3801.
#' \doi{10.1002/sim.1421}
#'
#' Lai X, Liu L (2018). "A Simple Test Procedure in Standardizing the Power of
#' Hosmer-Lemeshow Test in Large Data Sets." \emph{Journal of Statistical
#' Computation and Simulation}, \strong{88}(13), 2463--2472.
#' \doi{10.1080/00949655.2018.1467912}
#'
#' Nattino G, Finazzi S, Bertolini G (2014). "A New Calibration Test and a
#' Reappraisal of the Calibration Belt for the Assessment of Prediction Models
#' Based on Dichotomous Outcomes." \emph{Statistics in Medicine}, \strong{33}(14),
#' 2390--2407. \doi{10.1002/sim.6100}
#'
#' Zhang J, Ding J, Yang Y (2023). "Is a Classification Procedure Good Enough?
#' A Goodness-of-Fit Assessment Tool for Classification Learning."
#' \emph{Journal of the American Statistical Association}, \strong{118}(542),
#' 1115--1125. \doi{10.1080/01621459.2021.1979010}
#'
#' Liu Y, Xie J (2020). "Cauchy Combination Test: A Powerful Test with Analytic
#' p-Value Calculation under Arbitrary Dependency Structures." \emph{Journal of
#' the American Statistical Association}, \strong{115}(529), 393--402.
#' \doi{10.1080/01621459.2018.1554485}
#'
#' Hosmer DW, Hosmer T, le Cessie S, Lemeshow S (1997). "A Comparison of
#' Goodness-of-Fit Tests for the Logistic Regression Model." \emph{Statistics in
#' Medicine}, \strong{16}(9), 965--980.
#' \doi{10.1002/(sici)1097-0258(19970515)16:9<965::aid-sim509>3.0.co;2-o}
#'
#' The methods introduced by this package, and the studies that evaluate them,
#' are reported in the following. Reproduction materials for each are archived
#' and citable.
#'
#' Ebrahim EK, Khattab IG, El-Kotory A (2026). "A Modified Hosmer-Lemeshow
#' Goodness-of-Fit Test for Asymmetric Links: Second-Order Power and Robustness."
#' arXiv:2607.15454 [stat.ME]. \doi{10.48550/arXiv.2607.15454}
#'
#' Ebrahim EK, El-Kotory A (2026). "Benchmarking Goodness-of-Fit and Calibration
#' Algorithms for Logistic Regression Classifiers: A Large-Scale Simulation
#' Study under Sparse Data." \emph{Journal of Intelligent Computing and Data
#' Science}, \bold{2}(2), 97-128. \doi{10.21608/jicds.2026.505341.1048}
#' Reproduction materials: \doi{10.5281/zenodo.21286171}
#'
#' Ebrahim EK (2026). "Goodness-of-Fit Tests and Calibration Machine-Learning
#' Algorithms for Logistic Regression with Sparse Data." M.Sc. thesis,
#' Alexandria University. arXiv:2608.11140 [stat.ME].
#' \doi{10.48550/arXiv.2608.11140}
#'
#' Ebrahim EK, Hussein OAE-A, El-Kotory A (2026). "A Grouped Calibration Test for Logistic
#' Regression That Tolerates a Few Corrupted Records." arXiv:2608.20511 [stat.ME].
#' \doi{10.48550/arXiv.2608.20511} Reproduction materials (first version):
#' \doi{10.5281/zenodo.21247541}
#'
#' Ebrahim EK (2026). "EDGES: A Selection-Free Ensemble Goodness-of-Fit Test."
#' Reproduction materials. \doi{10.5281/zenodo.21320865}
#'
#' Ebrahim EK (2026). "Detection Subspaces: A Theory of Goodness-of-Fit Tests."
#' Reproduction materials. \doi{10.5281/zenodo.21687498}
#'
#' Ebrahim EK (2026). "Shrinkage Invalidates the Hosmer-Lemeshow Test: Goodness of Fit
#' for Penalized Logistic Regression, with an Application to Glaucoma Diagnosis."
#' arXiv:2609.06413 [stat.ME]. \doi{10.48550/arXiv.2609.06413} Reproduction materials
#' for \code{\link{shrink.gof}}: \doi{10.5281/zenodo.21900114}
#' @importFrom stats fitted predict model.matrix model.frame coef deviance pchisq binomial glm.fit kmeans median dist anova lm pnorm
#' @importFrom utils capture.output
#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept test battery
#' @concept Hosmer-Lemeshow
#' @concept binary classification
#' @export
run.all.gof <- function(object, predicted_probs = NULL, X = NULL,
                        tests = "all", G = 10, include_slow = TRUE,
                        parallel = FALSE, ncores = NULL,
                        calibration_plot = FALSE, install = c("ask", "no", "yes"),
                        control = list()) {

  install <- match.arg(install)
  if (!identical(G, "auto") && (!is.numeric(G) || length(G) != 1))
    stop("run.all.gof: 'G' must be a single number or 'auto'.")
  ctx <- .gof_context(object, predicted_probs, X, G = G)
  all_tests <- identical(tests, "all")
  sel <- if (all_tests) {
    # rows run only on request: GiViTI-external assumes the model was developed elsewhere
    names(.GOF_REGISTRY)[!vapply(.GOF_REGISTRY, function(e) isTRUE(e$on_request), logical(1))]
  } else intersect(tests, names(.GOF_REGISTRY))
  if (length(sel) == 0) stop("run.all.gof: no known tests selected.")
  if (all_tests) sel <- setdiff(sel, .gof_repeated_rows(ctx, control))

  # Optional PSOCK cluster for the slow bootstrap loops (Stute-Zhu, Lai-Liu-HL).
  # Reproducible: clusterSetRNGStream is seeded deterministically from the
  # session RNG, so same set.seed() + same ncores => identical results. The
  # sequential default path is completely untouched.
  if (isTRUE(parallel) && isTRUE(include_slow) && ctx$has_model &&
      length(intersect(sel, c("Stute-Zhu", "Lai-Liu-HL"))) > 0) {
    nc <- if (!is.null(ncores)) max(1L, as.integer(ncores)[1]) else {
      dc <- parallel::detectCores()
      if (is.na(dc)) 2L else max(1L, dc - 1L)
    }
    if (nc >= 2L) {
      cl <- tryCatch(parallel::makePSOCKcluster(nc), error = function(e) NULL)
      if (is.null(cl)) {
        message("run.all.gof: could not start a PSOCK cluster; ",
                "running sequentially.")
      } else {
        on.exit(parallel::stopCluster(cl), add = TRUE)
        parallel::clusterCall(cl, function(lp) { .libPaths(lp); NULL }, .libPaths())
        parallel::clusterSetRNGStream(cl, iseed = sample.int(.Machine$integer.max, 1L))
        ctx$cl <- cl
      }
    }
  }
  # If a test in this run needs an optional package that isn't installed, offer
  # to install it (interactive sessions only, with confirmation -- see
  # gof_install_suggests). Non-interactive runs are never touched.
  .gof_maybe_install(sel, include_slow, install)
  if (all_tests && isTRUE(include_slow))
    message("run.all.gof: running the full battery, including the slow tests ",
            "(le-Cessie, the GAM tests, Stute-Zhu, eHL, BAGofT, Projection, GiViTI). ",
            "For a quick run with the fast tests only, set include_slow = FALSE.")

  rows <- list(); skipped_model <- FALSE
  for (nm in sel) {
    e <- .GOF_REGISTRY[[nm]]
    if (isTRUE(e$slow) && !include_slow) next
    if (isTRUE(e$needs_model) && !ctx$has_model && is.null(ctx$X)) {
      skipped_model <- TRUE
      rows[[nm]] <- data.frame(Test = nm, Family = e$family, Statistic = NA_real_,
                               df = NA_real_, p_value = NA_real_,
                               Note = "Not applicable: needs a fitted glm", stringsAsFactors = FALSE)
      next
    }
    res <- tryCatch(e$fn(ctx, control[[nm]]),
                    error = function(err) list(Statistic = NA, df = NA, p_value = NA,
                                               Note = paste("error:", conditionMessage(err))))
    rows[[nm]] <- data.frame(Test = nm, Family = e$family,
                             Statistic = .gof_num(res$Statistic), df = .gof_num(res$df),
                             p_value = .gof_num(res$p_value),
                             Note = if (is.null(res$Note)) "" else res$Note,
                             stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)

  # Drop redundant "(grouped)" rows when covariate patterns do not repeat. On
  # fully sparse data (every observation its own pattern) the grouped and
  # per-observation forms of Pearson / deviance / McCullagh are the same test, so
  # the "(grouped)" row is a duplicate. It is kept only when some pattern holds
  # more than one record; the numerical comparison is the fallback when the
  # patterns cannot be counted.
  .eq_gof <- function(x, y) (is.na(x) && is.na(y)) ||
    (!is.na(x) && !is.na(y) && isTRUE(all.equal(x, y, tolerance = 1e-8)))
  gi <- grep(" \\(grouped\\)$", out$Test)
  J  <- if (length(gi) && !is.null(ctx$X))
          tryCatch(.gof_patterns(ctx$X, ctx$ph, ctx$y)$J, error = function(e) NA_integer_)
        else NA_integer_
  drop_rows <- integer(0)
  for (k in gi) {
    base <- sub(" \\(grouped\\)$", "", out$Test[k])
    bi <- which(out$Test == base)
    if (isTRUE(J == ctx$n) ||
        (length(bi) == 1L && .eq_gof(out$Statistic[k], out$Statistic[bi]) &&
         .eq_gof(out$p_value[k], out$p_value[bi])))
      drop_rows <- c(drop_rows, k)
  }
  if (length(drop_rows)) out <- out[-drop_rows, , drop = FALSE]

  # ensemble rows (only when a model is available and running the full set)
  if (ctx$has_model && identical(tests, "all")) {
    # a few-events or no-fit warning is already in the Note of the directed rows, and EF's
    # size-rule warning in the Note of its own row
    quiet <- function(expr)
      withCallingHandlers(expr, def_few_events = function(w) invokeRestart("muffleWarning"),
                          def_degenerate = function(w) invokeRestart("muffleWarning"),
                          ef_sparse_groups = function(w) invokeRestart("muffleWarning"))
    v3 <- tryCatch(quiet(def.ensemble.gof(ctx$model, G = ctx$G))$p_value, error = function(e) NA_real_)
    vu <- tryCatch(quiet(def.ensemble.gof(ctx$model, add_ef = TRUE, G = ctx$G))$p_value,
                   error = function(e) NA_real_)
    # these rows combine the unit form at the battery G; say so when control has moved the
    # three DEF rows they combine to other weights or another G
    moved <- any(vapply(control[c("DEF.poly2", "DEF.poly3", "DEF.stukel")], function(o) {
      g <- if (identical(o$G, "auto")) .def_auto_G(ctx$n) else o$G
      (!is.null(o$weights) && !identical(o$weights, "unit")) || (!is.null(g) && !isTRUE(all(g == ctx$G)))
    }, logical(1)))
    out <- rbind(out, data.frame(
      Test = c("Ensemble.Vote(3DEF)", "Ensemble.Univ(3DEF+EF)"), Family = "Ensemble",
      Statistic = NA_real_, df = NA_real_, p_value = c(v3, vu),
      Note = paste0("Cauchy combination of the directed tests",
                    if (moved) sprintf("; unit form, G = %s", format(ctx$G))),
      stringsAsFactors = FALSE))
  }

  if (skipped_model)
    message("Only prediction-based tests were run. Pass the fitted glm (or X) ",
            "to also run the directed and refit-based tests.")
  rownames(out) <- NULL
  class(out) <- c("gof_battery", "data.frame")

  # Optional GiViTI calibration belt: compute it in an isolated subprocess,
  # store it on the result, and draw it.
  if (isTRUE(calibration_plot) && "GiViTI" %in% sel) {
    belt <- .giviti_belt(ctx, control[["GiViTI"]])
    if (is.null(belt)) {
      message("run.all.gof: GiViTI calibration belt unavailable ",
              "(needs the 'givitiR' and 'callr' packages).")
    } else {
      attr(out, "giviti_belt") <- belt
      try(plot(belt), silent = TRUE)
    }
  }
  out
}

# Compute the GiViTI calibration belt object in an isolated callr subprocess
# (same crash-safety as the GiViTI test). Returns NULL if unavailable.
.giviti_belt <- function(ctx, opts = list()) {
  if (!requireNamespace("givitiR", quietly = TRUE) ||
      !requireNamespace("callr", quietly = TRUE))
    return(NULL)
  devel <- if (is.null(opts$devel)) "internal" else opts$devel
  o <- as.numeric(ctx$y); e <- as.numeric(ctx$ph)
  rng <- .rng_snapshot(); on.exit(.rng_restore(rng), add = TRUE)   # callr draws from the session stream
  tryCatch(
    callr::r(function(o, e, devel) {
      givitiR::givitiCalibrationBelt(o = o, e = e, devel = devel)
    }, args = list(o = o, e = e, devel = devel), timeout = 120),
    error = function(err) NULL)
}

#' Plot the GiViTI calibration belt from a goodness-of-fit battery
#'
#' Draws the GiViTI calibration belt stored on a \code{\link{run.all.gof}} result
#' that was produced with \code{calibration_plot = TRUE}. The belt shows the
#' fitted calibration curve with a confidence region against the 45-degree line.
#'
#' @param x A \code{gof_battery} object from \code{\link{run.all.gof}}.
#' @param ... Passed to the \pkg{givitiR} plot method.
#' @return \code{x}, invisibly.
#' @examples
#' \donttest{
#' fit <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial)
#' res <- run.all.gof(fit, tests = "GiViTI", calibration_plot = TRUE, install = "no")
#' plot(res)
#' }
#' @importFrom graphics plot
#' @exportS3Method plot gof_battery
plot.gof_battery <- function(x, ...) {
  belt <- attr(x, "giviti_belt")
  if (is.null(belt)) {
    message("No GiViTI calibration belt is stored on this result. ",
            "Re-run run.all.gof(..., calibration_plot = TRUE).")
    return(invisible(NULL))
  }
  if (!requireNamespace("givitiR", quietly = TRUE)) {
    message("Install the 'givitiR' package to draw the calibration belt.")
    return(invisible(NULL))
  }
  plot(belt, ...)
  invisible(x)
}

#' Print a goodness-of-fit battery
#'
#' Formats the \code{\link{run.all.gof}} result as a compact, readable table:
#' rows grouped by test family, p-values shown to four decimals (or scientific
#' for very small values, \code{"-"} when not available), and a significance
#' flag. The object is still a plain \code{data.frame} underneath, so all the
#' raw columns remain available for programmatic use.
#'
#' The header counts the rows that reject at the 0.05 level. The rows printed for
#' comparison only (the sparse \code{Pearson} and \code{Deviance} statistics and the
#' \code{F-test}; see \code{\link{run.all.gof}}) carry no significance flag and are not
#' counted. Notes are wrapped to \code{getOption("width")}.
#'
#' @param x A \code{gof_battery} object returned by \code{\link{run.all.gof}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @examples
#' fit <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial)
#' print(run.all.gof(fit, include_slow = FALSE, install = "no"))
#' @exportS3Method print gof_battery
print.gof_battery <- function(x, ...) {
  d <- as.data.frame(x)
  ord <- match(d$Family, .battery_family_order); ord[is.na(ord)] <- 99L
  d <- d[order(ord, seq_len(nrow(d))), , drop = FALSE]

  ## distinct notes -> footnote markers [a], [b], ... (keeps the table narrow)
  notes <- ifelse(is.na(d$Note), "", d$Note)
  uniq  <- unique(notes[notes != ""])
  marks <- letters[seq_along(uniq)]; names(marks) <- uniq
  mk    <- ifelse(notes == "", "", paste0(" [", marks[notes], "]"))

  fp <- vapply(d$p_value, function(p)
    if (is.na(p)) "-" else if (p < 1e-4) sprintf("%.1e", p)
    else formatC(p, format = "f", digits = 4), character(1))
  counted <- .gof_counted(d)
  sig <- vapply(d$p_value, function(p)
    if (is.na(p)) "" else if (p < .001) "***" else if (p < .01) "**"
    else if (p < .05) "*" else if (p < .1) "." else "", character(1))
  sig[!counted] <- ""                     # comparison rows carry no flag
  fnum <- function(v, dig) vapply(v, function(z)
    if (is.na(z)) "" else formatC(z, format = "g", digits = dig), character(1))
  st <- fnum(d$Statistic, 4)
  df <- vapply(d$df, function(z)           # whole-number df in full (2997, not 3e+03)
    if (is.na(z)) "" else if (abs(z - round(z)) < 1e-8) sprintf("%.0f", z)
    else formatC(z, format = "g", digits = 3), character(1))

  tname <- paste0(d$Test, mk)
  wT <- max(nchar(tname), nchar("Test"))
  wS <- max(nchar(st), nchar("Statistic"))
  wD <- max(nchar(df), 2L)
  wP <- max(nchar(fp), nchar("p-value")); wG <- 3L
  tw <- 1L + wT + 2L + wS + 1L + wD + 2L + wP + 1L + wG   # table width

  padL <- function(s, w) formatC(s, width = -w, flag = " ")  # left-justified
  padR <- function(s, w) formatC(s, width = w)               # right-justified
  row  <- function(t, s, dd, p, g, tag = "")
    cat(" ", padL(t, wT), "  ", padR(s, wS), " ", padR(dd, wD), "  ",
        padR(p, wP), " ", padL(g, wG), tag, "\n", sep = "")

  nrej <- sum(counted & !is.na(d$p_value) & d$p_value < .05)
  cat(sprintf("\nGoodness-of-fit battery: %d %s  (%d reject at 0.05)\n",
              nrow(d), ngettext(nrow(d), "test", "tests"), nrej))
  cat(strrep("=", tw), "\n", sep = "")
  row("Test", "Statistic", "df", "p-value", "")
  fam_now <- ""
  for (i in seq_len(nrow(d))) {
    if (!identical(d$Family[i], fam_now)) {
      fam_now <- d$Family[i]
      lab <- paste0("--- ", fam_now, " ")
      cat(" ", lab, strrep("-", max(0L, tw - nchar(lab) - 1L)), "\n", sep = "")
    }
    row(tname[i], st[i], df[i], fp[i], sig[i])
  }
  cat(strrep("-", tw), "\n", sep = "")
  cat(" Signif.:  *** <.001   ** <.01   * <.05   . <.1\n")
  if (length(uniq)) {
    cat(" Notes:\n")
    w <- max(30L, getOption("width", 80L) - 1L)
    for (u in uniq)
      cat(strwrap(u, width = w - 7L, initial = paste0("   [", marks[u], "] "),
                  prefix = strrep(" ", 7L)), sep = "\n")
  }
  invisible(x)
}

# Coerce a possibly-NULL/empty scalar to a numeric (NA if missing).
.gof_num <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x)[1]

# Rows of the full battery that would repeat another row exactly: DEF.poly3 is the cubic
# basis at the battery's G, which is EDGE.G10 when that G is 10, and EDGE.G10 is EDGE
# itself when the default partition has ten groups (n <= 250). Rows that control moves to
# other weights or another G are kept.
.gof_repeated_rows <- function(ctx, control) {
  untouched <- function(nm) is.null(control[[nm]]$weights) && is.null(control[[nm]]$G)
  drop <- character(0)
  if (isTRUE(ctx$G == 10) && untouched("DEF.poly3") && untouched("EDGE.G10"))
    drop <- c(drop, "DEF.poly3")
  if (.def_auto_G(ctx$n) == 10 && untouched("EDGE.G10") && untouched("EDGE"))
    drop <- c(drop, "EDGE.G10")
  drop
}

# Rows printed for comparison only: on ungrouped binary data the sparse Pearson and
# deviance statistics have no chi-square reference, and the F-test on deviance residuals
# is liberal. They are not counted among the rejections of print() and summary().
.battery_comparison_only <- c("Pearson", "Deviance", "F-test")
.gof_counted <- function(d) !(d$Test %in% .battery_comparison_only)

# The X of the (y, predicted_probs, X) form must be the model matrix of the fit: numeric,
# intercept first. A data frame is coded with model.matrix(), and a matrix without a
# constant column gets one, so the directed and covariate-space rows see the same design
# as from the fitted glm.
.gof_design <- function(X, n) {
  if (is.data.frame(X)) X <- stats::model.matrix(~ ., data = X)
  X <- as.matrix(X)
  if (!is.numeric(X) || NROW(X) != n)
    stop("run.all.gof: 'X' must be numeric (or a data frame) with one row per record.", call. = FALSE)
  if (anyNA(X)) stop("run.all.gof: 'X' has missing values.", call. = FALSE)
  const <- which(apply(X, 2L, function(v) all(abs(v - 1) < 1e-12)))
  if (!length(const)) X <- cbind(`(Intercept)` = 1, X)
  else if (const[1] != 1L) X <- X[, c(const[1], setdiff(seq_len(ncol(X)), const[1])), drop = FALSE]
  X
}

# Warn when the data cannot support the battery: fitted risks of 0 or 1 (separation)
# make the standardized and partition statistics explode, and very few events or
# non-events leave the grouped references without support.
.gof_check_support <- function(y, ph) {
  if (any(ph < 1e-8 | ph > 1 - 1e-8))
    warning("run.all.gof: some fitted risks are 0 or 1 to machine precision; the data may be ",
            "separated, and the standardized, partition and directed rows are unreliable.",
            call. = FALSE)
  k <- min(sum(y), length(y) - sum(y))
  if (k < 10)
    warning(sprintf(paste0("run.all.gof: only %d %s; the large-sample references of the ",
                           "battery are unreliable at this size."),
                    k, if (sum(y) <= length(y) - sum(y)) "events" else "non-events"),
            call. = FALSE)
  invisible(NULL)
}

# Build the one context object every test reads from.
.gof_context <- function(object, predicted_probs = NULL, X = NULL, G = 10) {
  ctx <- if (inherits(object, "glm")) {
    if (object$family$family != "binomial")
      stop("run.all.gof: the model must be a binomial glm.")
    y  <- as.numeric(object$y)
    if (!all(y %in% c(0, 1)) || any(object$prior.weights != 1))
      stop("run.all.gof: the tests need one row per binary (0/1) record and no prior weights; ",
           "expand grouped (cbind) or weighted data to one row per record first.", call. = FALSE)
    if (!identical(object$family$link, "logit"))
      warning("run.all.gof: the model uses the '", object$family$link, "' link; most tests ",
              "here assume a logit fit (the Stukel row says so in its note).", call. = FALSE)
    .gof_check_support(y, as.numeric(object$fitted.values))
    ph <- pmin(pmax(as.numeric(object$fitted.values), 1e-6), 1 - 1e-6)
    list(y = y, ph = ph, X = stats::model.matrix(object),
         data = tryCatch(stats::model.frame(object), error = function(e) NULL),
         model = object, G = G, n = length(y), has_model = TRUE,
         p = length(stats::coef(object)))
  } else {
    if (!is.numeric(object))
      stop("run.all.gof: 'object' must be a fitted binomial glm or a numeric (0/1) y vector.")
    y <- as.numeric(object)
    if (is.null(predicted_probs))
      stop("run.all.gof: supply 'predicted_probs' when 'object' is not a glm.")
    ph <- pmin(pmax(as.numeric(predicted_probs), 1e-6), 1 - 1e-6)
    if (length(y) != length(ph))
      stop("run.all.gof: 'object' (y) and 'predicted_probs' lengths differ.")
    if (!all(y %in% c(0, 1)))
      stop("run.all.gof: 'object' must be a binary (0/1) vector or a glm.")
    if (!is.null(X)) X <- .gof_design(X, length(y))
    .gof_check_support(y, as.numeric(predicted_probs))
    list(y = y, ph = ph, X = X, data = NULL, model = NULL, G = G, n = length(y),
         has_model = FALSE, p = if (is.null(X)) NA_integer_ else ncol(X))
  }
  if (identical(G, "auto")) {                      # resolved once, so every row uses the same G
    ctx$G <- .def_auto_G(ctx$n)
    ctx$G_auto <- TRUE
  }
  ctx
}

# Equal-frequency grouping of the predicted probabilities into G groups.
.gof_groups_ef <- function(ph, G) {
  n <- length(ph)
  pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
}

# Hosmer-Lemeshow statistic for a given grouping; drops empty/degenerate groups.
.gof_hl_stat <- function(y, ph, grp) {
  idx <- split(seq_along(y), grp)
  O  <- vapply(idx, function(I) sum(y[I]),  numeric(1))
  E  <- vapply(idx, function(I) sum(ph[I]), numeric(1))
  ng <- vapply(idx, length, numeric(1))
  keep <- ng > 0 & E > 1e-8 & E < ng - 1e-8
  hl <- sum((O[keep] - E[keep])^2 / (E[keep] * (1 - E[keep] / ng[keep])))
  list(stat = hl, df = sum(keep) - 2)
}

# ---- test wrappers (internal): each returns list(Statistic, df, p_value, Note) ----

gof_pearson <- function(ctx, opts = list()) {
  V  <- ctx$ph * (1 - ctx$ph)
  X2 <- sum((ctx$y - ctx$ph)^2 / V)
  df <- if (ctx$has_model) ctx$model$df.residual else ctx$n - 1
  list(Statistic = X2, df = df, p_value = stats::pchisq(X2, df, lower.tail = FALSE),
       Note = if (ctx$has_model) "sparse form; comparison only, not counted"
              else "df = n-1 (no model); comparison only, not counted")
}

gof_deviance <- function(ctx, opts = list()) {
  if (ctx$has_model) {
    D <- ctx$model$deviance; df <- ctx$model$df.residual
  } else {
    D  <- -2 * sum(ctx$y * log(ctx$ph) + (1 - ctx$y) * log(1 - ctx$ph)); df <- ctx$n - 1
  }
  list(Statistic = D, df = df, p_value = stats::pchisq(D, df, lower.tail = FALSE),
       Note = if (ctx$has_model) "sparse form; comparison only, not counted"
              else "df = n-1 (no model); comparison only, not counted")
}

# Osius-Rojek normal-approximation test, matching LogisticDx::gof.glm exactly.
# The statistic is defined on the J DISTINCT COVARIATE PATTERNS, not the raw rows:
# z = (PrG - (J - p)) / sqrt(A1 + RSS1), where PrG is the grouped Pearson statistic,
# A1 = 2*(J - sum(1/n_g)) is the excess-variance correction over the pattern sizes
# n_g, and RSS1 is the residual SS of the WLS regression of (1-2P_g)/(n_g P_g(1-P_g))
# on the pattern design matrix with weights n_g P_g(1-P_g). For fully sparse data
# (every pattern unique, n_g = 1) this reduces to A1 = 0 and J = N, so simulation
# results are unchanged; only data with repeated covariate patterns differ.
gof_osius <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y; p <- ncol(X)
  # aggregate observations sharing an identical covariate pattern (row of X)
  key  <- apply(round(X, 10), 1L, paste, collapse = "\r")
  ug   <- !duplicated(key)
  nm   <- key[ug]
  ng   <- as.numeric(tapply(rep(1, length(y)), key, sum)[nm])   # trials per pattern
  ysum <- as.numeric(tapply(y,  key, sum)[nm])                  # events per pattern
  Pg   <- as.numeric(tapply(ph, key, `[`, 1L)[nm])              # fitted prob per pattern
  Xu   <- X[ug, , drop = FALSE]                                 # one row per pattern
  J    <- length(ng)
  Vg   <- ng * Pg * (1 - Pg)
  PrG  <- sum((ysum - ng * Pg)^2 / Vg)                          # grouped Pearson
  A1   <- 2 * (J - sum(1 / ng))                                 # = 0 only if all n_g = 1
  cvar <- (1 - 2 * Pg) / Vg
  fit  <- tryCatch(stats::lm.wfit(Xu, cvar, Vg), error = function(e) NULL)
  if (is.null(fit))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: weighted least-squares step failed"))
  RSS1 <- sum(Vg * fit$residuals^2)
  varz <- A1 + RSS1
  if (!is.finite(varz) || varz <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: non-positive variance"))
  z <- (PrG - (J - p)) / sqrt(varz)
  list(Statistic = z, df = NA_real_, p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
       Note = if (J < length(y)) sprintf("grouped to %d covariate patterns", J) else "normal deviate, two-sided")
}

# Copas (1989) unweighted residual-sum-of-squares test. Binary expansion is
# trivial (one trial per observation). Ported from goflogit (7 tests.R).
# The variance c'(I - W X (X'WX)^{-1} X') W c of goflogit is evaluated as
# c'Wc - (c'WX)(X'WX)^{-1}(X'Wc): O(n p^2) time and O(n p) memory, where the
# n x n form cost O(n^2 p) and O(n^2). The statistic is a normal deviate,
# referred two-sided (the chi-square(1) of z^2 that goflogit reports).
gof_copas <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  V <- ph * (1 - ph)
  copas    <- sum((y - ph)^2)
  meacopas <- sum(V)
  c12 <- 1 - 2 * ph
  b   <- crossprod(X, V * c12)                          # X'Wc
  sb  <- tryCatch(solve(crossprod(X, V * X), b), error = function(e) NULL)
  if (is.null(sb))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: singular information matrix"))
  varcopas <- sum(V * c12^2) - sum(b * sb)
  if (!is.finite(varcopas) || varcopas <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: non-positive variance"))
  z <- (copas - meacopas) / sqrt(varcopas)
  list(Statistic = z, df = NA_real_, p_value = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
       Note = "normal deviate, two-sided")
}

# White/Orme information-matrix (IM) test. Explained sum of squares from
# regressing the Pearson residuals on the auxiliary regressors [sqrt(V) X |
# sqrt(V)(1-2p) X^2]; statistic ~ chi-square_{ncol(X)}. From IM (infromation
# Matrix).R (IMtest_fast); matches the thesis simulation.
gof_im <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  if (any(ph < 1e-8) || any(ph > 1 - 1e-8))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: fitted probabilities too extreme (near 0 or 1)"))
  r       <- (y - ph) / sqrt(ph * (1 - ph))
  w_sqrt  <- sqrt(ph * (1 - ph))
  W_aux   <- cbind(w_sqrt * X, (w_sqrt * (1 - 2 * ph)) * (X^2))
  xtx_inv <- tryCatch(solve(crossprod(W_aux)), error = function(e) NULL)
  if (is.null(xtx_inv))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: singular auxiliary matrix"))
  Wr <- crossprod(W_aux, r)
  im <- as.numeric(crossprod(Wr, xtx_inv %*% Wr))
  df <- ncol(X)
  list(Statistic = im, df = df, p_value = stats::pchisq(im, df, lower.tail = FALSE), Note = "")
}

gof_hl <- function(ctx, opts = list()) {
  h <- .gof_hl_stat(ctx$y, ctx$ph, .gof_groups_ef(ctx$ph, ctx$G))
  if (h$df < 1) return(list(Statistic = h$stat, df = h$df, p_value = NA_real_,
                            Note = "too few non-empty groups"))
  list(Statistic = h$stat, df = h$df,
       p_value = stats::pchisq(h$stat, h$df, lower.tail = FALSE), Note = "")
}

# Nattino, Pennell & Lemeshow (2020, Biometrics 76:549) large-sample Hosmer-Lemeshow test. In large samples
# the ordinary test rejects any misfit, however small, because its noncentrality grows with n. They test
# H0: eps <= eps0 instead, where eps = sqrt(lambda / n) does not grow with n, by referring the ordinary
# statistic to a noncentral chi-square on G - 2 df with noncentrality eps0^2 n. eps0 is the misfit that
# would be just significant at n0 = 10^6 (their convention): eps0^2 = (qchisq(.95, G - 2) - (G - 2)) / n0.
# Reproduces their application: C = 25.35, n = 315,828, G = 10 gives p = 0.010 (0.001 for the ordinary test).
.hl_largeN_p <- function(C, n, df, n0 = 1e6) {
  eps0_sq <- (stats::qchisq(0.95, df) - df) / n0
  list(p = stats::pchisq(C, df, ncp = eps0_sq * n, lower.tail = FALSE), eps0 = sqrt(eps0_sq),
       eps_hat = sqrt(max(C - df, 0) / n))
}

gof_hl_largeN <- function(ctx, opts = list()) {
  n0 <- if (is.null(opts$n0)) 1e6 else opts$n0
  h <- .gof_hl_stat(ctx$y, ctx$ph, .gof_groups_ef(ctx$ph, ctx$G))
  if (h$df < 1) return(list(Statistic = h$stat, df = h$df, p_value = NA_real_,
                            Note = "too few non-empty groups"))
  r <- .hl_largeN_p(h$stat, length(ctx$y), h$df, n0)
  list(Statistic = h$stat, df = h$df, p_value = r$p,
       Note = sprintf("H0: eps <= %.2e (n0 = %g); eps_hat = %.2e", r$eps0, n0, r$eps_hat))
}

gof_hlw <- function(ctx, opts = list()) {
  br <- seq(0, 1, length.out = ctx$G + 1); br[1] <- -Inf; br[length(br)] <- Inf
  grp <- cut(ctx$ph, breaks = br, labels = FALSE)
  h <- .gof_hl_stat(ctx$y, ctx$ph, grp)
  if (h$df < 1) return(list(Statistic = h$stat, df = h$df, p_value = NA_real_,
                            Note = "too few non-empty equal-width bins"))
  list(Statistic = h$stat, df = h$df,
       p_value = stats::pchisq(h$stat, h$df, lower.tail = FALSE), Note = "")
}

# Pigeon-Heyse J2: Hosmer-Lemeshow-type statistic with a per-group variance
# correction phi_k and df = g - 1. Quantile (equal-frequency) groups. From
# pigeonheyse.R; needs only the response and fitted probabilities.
gof_ph_test <- function(ctx, opts = list()) {
  ph <- ctx$ph; y <- ctx$y; g <- ctx$G
  br  <- stats::quantile(ph, probs = seq(0, 1, length.out = g + 1))
  grp <- tryCatch(cut(ph, breaks = br, labels = FALSE, include.lowest = TRUE),
                  error = function(e) NULL)
  if (is.null(grp) || length(unique(grp[!is.na(grp)])) < 2)
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "could not form groups (ties in fitted probabilities)"))
  idx  <- split(seq_along(y), grp)
  Ok   <- vapply(idx, function(I) sum(y[I]),   numeric(1))
  nk   <- vapply(idx, length,                  numeric(1))
  pbar <- vapply(idx, function(I) mean(ph[I]), numeric(1))
  Vk   <- nk * pbar * (1 - pbar)
  phik <- vapply(idx, function(I) sum(ph[I] * (1 - ph[I])), numeric(1)) / Vk
  ok   <- is.finite(Vk) & Vk > 0 & is.finite(phik) & phik > 0
  J2   <- sum(((Ok[ok] - nk[ok] * pbar[ok])^2 / Vk[ok]) / phik[ok])
  df   <- length(idx) - 1
  list(Statistic = J2, df = df, p_value = stats::pchisq(J2, df, lower.tail = FALSE), Note = "")
}

gof_ef <- function(ctx, opts = list()) {
  note <- ""                                       # the size-rule warning goes to Note
  r <- withCallingHandlers(ef.gof(ctx$y, ctx$ph, G = ctx$G),   # chisq reference (package default)
                           ef_sparse_groups = function(w) {
                             note <<- "fewer than one half expected event or non-event in a group"
                             invokeRestart("muffleWarning")
                           })
  list(Statistic = r$Test_Statistic, df = r$df, p_value = r$p_value, Note = note)
}

# EF with the plug-in normal reference, the reference for many groups: it runs at
# G = "auto" (groups of about 25 records) whatever the battery's G, unless
# control = list("EF-normal" = list(G = ...)) says otherwise.
gof_ef_normal <- function(ctx, opts = list()) {
  if (!ctx$has_model && is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model or X"))
  G <- if (is.null(opts$G)) "auto" else opts$G
  r <- if (ctx$has_model) ef.gof(ctx$model, G = G, reference = "normal")
       else ef.gof(ctx$y, ctx$ph, G = G, reference = "normal", X = ctx$X)
  list(Statistic = r$Test_Statistic, df = NA, p_value = r$p_value,
       Note = sprintf("normal reference, many groups (G = %d)", as.integer(r$G)))
}

gof_def <- function(ctx, opts = list()) {
  if (!ctx$has_model && is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model or X"))
  b  <- if (is.null(opts$basis))   "poly3" else opts$basis
  wt <- if (is.null(opts$weights)) "unit"  else opts$weights
  G  <- if (is.null(opts$G)) ctx$G else if (identical(opts$G, "auto")) .def_auto_G(ctx$n) else opts$G
  few <- NULL                                      # def.gof's few-events and no-fit warnings go to Note
  keep_few <- function(w) {
    few <<- c(few, sub("^def.gof: ", "", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }
  r  <- if (ctx$has_model)
          withCallingHandlers(def.gof(ctx$model, G = G, basis = b, weights = wt), def_few_events = keep_few,
                              def_degenerate = keep_few, def_no_information = keep_few)
        else suppressWarnings(withCallingHandlers(
          def.gof(ctx$y, ctx$ph, X = ctx$X, G = G, basis = b, weights = wt), def_few_events = keep_few,
          def_degenerate = keep_few, def_no_information = keep_few))
  auto <- identical(opts$G, "auto") || (is.null(opts$G) && isTRUE(ctx$G_auto))
  note <- c(if (r$Method == "score") "score form",
            if (r$Method == "satterthwaite") "S/c ~ chi-square(df), Satterthwaite",
            if (auto) sprintf("G = %d (auto)", as.integer(G)), few)
  list(Statistic = r$Test_Statistic, df = r$df, p_value = r$p_value,
       Note = paste(note, collapse = "; "))
}

# Stukel (1988) test against the generalized logistic link. The two directions are the
# sign-split squared logits of LogisticDx::gof.glm ("SstBoth"): za = eta^2/2 where the
# fitted risk is at least one half, zb = -eta^2/2 where it is below.
#
# form = "joint" (the default from 2.8.0) is u'I^-1 u, with u = Z'(y - p) and I the
# information of Z after adjusting for the fitted coefficients: the Rao score test for
# adding both columns, which anova(..., test = "Rao") gives up to glm's convergence
# tolerance. Aliased columns of X are left out. When no fitted risk lies
# on one side of one half that column is identically zero, and the test is the 1-df score
# test on the other; Note says so. A column whose information after the fit is below 1e-10 of
# its information before the fit is one the model already spans and is also left out.
#
# form = "lr" refits the model with the non-zero columns added and refers the drop in
# deviance to chi-square on the number of columns the refit could estimate. Under separation
# the refit may not converge, and the row is then NA with a note.
#
# form = "marginal" is the statistic of 2.7.0 and earlier, and of LogisticDx: the two
# marginal score z's squared and summed. Once the model is fitted the two directions are
# correlated, so the sum is not chi-square(2) and is liberal: for example, a post-fit
# correlation of about -0.71 in a typical design, and a rejection rate of up to about 7% at
# the 5% level. It is kept only so that earlier results can be reproduced.
gof_stukel <- function(ctx, opts = list()) {
  form <- if (is.null(opts$form)) "joint" else match.arg(opts$form, c("joint", "marginal", "lr"))
  if (!ctx$has_model)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  ph  <- ctx$ph; y <- ctx$y; X <- ctx$X
  if (min(sum(y), length(y) - sum(y)) == 0)        # no maximum-likelihood fit, in any form
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = paste("Not run: no", if (sum(y) == 0) "events" else "non-events",
                             "in the sample, so the model has no maximum-likelihood fit")))
  eta <- as.numeric(ctx$model$linear.predictors)
  za  <- 0.5 * eta^2 * (ph >= 0.5)                 # Stukel direction, p >= 0.5
  zb  <- -0.5 * eta^2 * (ph < 0.5)                 # Stukel direction, p < 0.5

  if (form == "marginal") {                        # the 2.7.0 code, unchanged
    if (requireNamespace("statmod", quietly = TRUE)) {
      Z <- abs(statmod::glm.scoretest(ctx$model, cbind(za, zb)))
      chi <- sum(Z^2)
    } else {
      za_z <- .gof_score_z(za, y, ph, X)
      zb_z <- .gof_score_z(zb, y, ph, X)
      if (!is.finite(za_z) || !is.finite(zb_z))
        return(list(Statistic = NA, df = NA, p_value = NA, Note = "score test undefined"))
      chi <- za_z^2 + zb_z^2
    }
    return(list(Statistic = chi, df = 2, p_value = stats::pchisq(chi, 2, lower.tail = FALSE),
                Note = "marginal sum (LogisticDx SstBoth), not chi-square(2)"))
  }

  if (!identical(ctx$model$family$link, "logit") || any(ctx$model$prior.weights != 1) ||
      !all(y %in% c(0, 1)))
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "Not run: needs an unweighted logit fit to binary data"))
  Z    <- cbind(za, zb)
  keep <- colSums(Z != 0) > 0
  if (!any(keep))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: the linear predictor is zero"))
  Z    <- Z[, keep, drop = FALSE]
  side <- if (all(keep)) "" else if (keep[1]) "1 df: no fitted risk below 0.5"
          else "1 df: no fitted risk at or above 0.5"
  X    <- X[, !is.na(stats::coef(ctx$model)), drop = FALSE]   # aliased columns carry no information

  if (form == "lr") {
    off <- if (is.null(ctx$model$offset)) rep(0, length(y)) else ctx$model$offset
    f1  <- tryCatch(suppressWarnings(stats::glm.fit(cbind(X, Z), y, family = stats::binomial(),
                                                    offset = off)),
                    error = function(e) e)
    if (inherits(f1, "error"))
      return(list(Statistic = NA, df = NA, p_value = NA,
                  Note = paste("Not run: augmented fit failed:", conditionMessage(f1))))
    if (!isTRUE(f1$converged))
      return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: augmented fit did not converge"))
    k <- f1$rank - ctx$model$rank
    if (k < 1)
      return(list(Statistic = NA, df = NA, p_value = NA,
                  Note = "Not run: the Stukel columns are aliased with the model"))
    lr <- max(ctx$model$deviance - f1$deviance, 0)
    return(list(Statistic = lr, df = k, p_value = stats::pchisq(lr, k, lower.tail = FALSE), Note = side))
  }

  ph  <- as.numeric(ctx$model$fitted.values)       # unclamped, as the fit itself uses
  W   <- ph * (1 - ph)
  u   <- colSums(Z * (y - ph))
  ok  <- rep(TRUE, ncol(Z))
  chi <- tryCatch({
    ZWX <- crossprod(Z, W * X)
    I   <- crossprod(Z, W * Z) - ZWX %*% solve(crossprod(X, W * X), t(ZWX))
    # a direction the model already spans is left out: its information after the fit is below
    # 1e-10 of its information before the fit (every direction, when the fitted logit is constant)
    ok  <- diag(I) > 1e-10 * colSums(W * Z^2)
    if (any(ok)) as.numeric(crossprod(u[ok], solve(I[ok, ok, drop = FALSE], u[ok]))) else NA_real_
  }, error = function(e) NA_real_)
  if (!isTRUE(any(ok)))
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "Not run: no Stukel direction has information left after the fit"))
  if (!is.finite(chi) || chi < 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: singular score information"))
  k <- sum(ok)
  if (!all(ok)) side <- "1 df: the other Stukel direction has no information left after the fit"
  list(Statistic = chi, df = k, p_value = stats::pchisq(chi, k, lower.tail = FALSE), Note = side)
}

# Rao score-test z for adding one column to a fitted binomial glm. Equals
# statmod::glm.scoretest for the binomial family (dispersion = 1).
.gof_score_z <- function(z, y, ph, X) {
  W   <- ph * (1 - ph)
  U   <- sum(z * (y - ph))
  WzX <- crossprod(X, W * z)
  Vv  <- sum(W * z^2) - as.numeric(t(WzX) %*% solve(crossprod(X, W * X)) %*% WzX)
  if (!is.finite(Vv) || Vv <= 0) return(NA_real_)
  U / sqrt(Vv)
}

# Moore-Penrose pseudo-inverse via SVD (matches MASS::ginv; avoids a dependency).
.gof_ginv <- function(M, tol = sqrt(.Machine$double.eps)) {
  s <- svd(M)
  pos <- s$d > max(tol * s$d[1], 0)
  if (!any(pos)) return(matrix(0, ncol(M), nrow(M)))
  s$v[, pos, drop = FALSE] %*% (t(s$u[, pos, drop = FALSE]) / s$d[pos])
}

# k-means with a fixed seed (123) that does not disturb the caller's RNG state.
.gof_kmeans <- function(mat, centers, nstart = 1L) {
  has <- exists(".Random.seed", envir = .GlobalEnv)
  old <- if (has) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(123)
  cl <- stats::kmeans(mat, centers = centers, nstart = nstart)$cluster
  if (has) assign(".Random.seed", old, envir = .GlobalEnv)
  cl
}

# Tsiatis (1980) clustering score test: cluster the covariate space, then score-
# test the cluster indicators added to the model. Ported from Tsiatis.R.
gof_tsiatis <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  cov_mat <- X[, -1, drop = FALSE]                       # drop intercept column
  if (ncol(cov_mat) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no covariates to cluster"))
  cl <- tryCatch(.gof_kmeans(cov_mat, ctx$G), error = function(e) NULL)
  if (is.null(cl))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "clustering failed"))
  Xc <- stats::model.matrix(~ factor(cl))[, -1, drop = FALSE]
  if (ncol(Xc) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "only one non-empty cluster"))
  W   <- ph * (1 - ph)
  U   <- colSums(Xc * (y - ph))
  V11 <- crossprod(X, W * X)
  V12 <- crossprod(X, W * Xc)
  V22 <- crossprod(Xc, W * Xc)
  V   <- V22 - t(V12) %*% .gof_ginv(V11) %*% V12
  Tstat <- as.numeric(t(U) %*% .gof_ginv(V) %*% U)
  rankV <- sum(eigen(V, symmetric = TRUE, only.values = TRUE)$values > 1e-8)
  list(Statistic = Tstat, df = rankV,
       p_value = stats::pchisq(Tstat, rankV, lower.tail = FALSE), Note = "")
}

# Xie covariate-space grouped chi-square (own group rule, fractional df). From Xie.R.
gof_xie <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y
  k <- ncol(X) - 1
  cov_mat <- X[, -1, drop = FALSE]
  if (ncol(cov_mat) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no covariates to cluster"))
  G <- if (k < 5) 10 else k + 5
  cl <- tryCatch(.gof_kmeans(cov_mat, G, nstart = 25L), error = function(e) NULL)
  if (is.null(cl))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "clustering failed"))
  stat <- 0
  for (I in split(seq_along(y), cl)) {
    ng <- length(I); pbar <- mean(ph[I])
    if (pbar > 0 && pbar < 1)
      stat <- stat + (sum(y[I]) - ng * pbar)^2 / (ng * pbar * (1 - pbar))
  }
  df <- G - k / 2 - 1
  if (df <= 0) return(list(Statistic = stat, df = df, p_value = NA, Note = "df <= 0"))
  list(Statistic = stat, df = df, p_value = stats::pchisq(stat, df, lower.tail = FALSE), Note = "")
}

# Pulkstenis-Robinson: covariate patterns from categorical vars, split by median
# fitted prob, chi-square on the 2*M subgroups. Base-R port of PR_test_only.R.
gof_pr <- function(ctx, opts = list()) {
  if (!ctx$has_model || is.null(ctx$data))
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "needs a glm model (for categorical covariates)"))
  # model.frame stores the response in column 1; the rest are the model terms, and
  # "(offset)" or "(weights)" columns, which are not covariates.
  cand <- names(ctx$data)[-1]
  cand <- cand[!startsWith(cand, "(")]
  maxlev <- getOption("ebrahim.gof.pr.maxlev", 6)
  cat_vars <- if (!is.null(opts$cat_var)) opts$cat_var else
    cand[vapply(cand, function(v) {
      col <- ctx$data[[v]]
      is.factor(col) || is.character(col) || is.logical(col) ||
        (is.numeric(col) && length(unique(col)) <= maxlev)
    }, logical(1))]
  if (length(cat_vars) == 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs a categorical covariate"))

  y <- ctx$y; ph <- ctx$ph
  patt <- do.call(paste, c(lapply(cat_vars, function(v) as.character(ctx$data[[v]])), sep = "_"))
  M <- length(unique(patt))
  lev <- character(length(y))
  for (pp in unique(patt)) {
    ix <- which(patt == pp)
    lev[ix] <- ifelse(ph[ix] <= stats::median(ph[ix]), "low", "high")
  }
  idx <- split(seq_along(y), paste(patt, lev, sep = "::"))
  os <- vapply(idx, function(I) sum(y[I] == 1), numeric(1))
  of <- vapply(idx, function(I) sum(y[I] == 0), numeric(1))
  es <- vapply(idx, function(I) sum(ph[I]), numeric(1))
  ef <- vapply(idx, function(I) sum(1 - ph[I]), numeric(1))
  keep <- es > 0 & ef > 0
  chisq <- sum((os[keep] - es[keep])^2 / es[keep] + (of[keep] - ef[keep])^2 / ef[keep])
  df <- 2 * M - length(cat_vars) - 2
  if (df <= 0)
    return(list(Statistic = chisq, df = df, p_value = NA, Note = "df <= 0 (too few patterns)"))
  list(Statistic = chisq, df = df, p_value = stats::pchisq(chisq, df, lower.tail = FALSE),
       Note = paste0("split on categorical: ", paste(cat_vars, collapse = ", ")))
}

# le Cessie-van Houwelingen smoothed-residual GOF test (general, multivariate).
# Reference: le Cessie, S. & van Houwelingen, H.C. (1995), Biometrics 51:600-614.
# Adapted (with attribution) from the USGS 'smwrStats' package leCessie.test(),
# which is a work of the US federal government (public domain), with two departures
# from that source: the moment reference is transposed to match the paper (see the
# note inside), and the moments are evaluated in O(n^2 p) rather than O(n^3). It still
# builds an n-by-n kernel matrix, so it stays O(n^2) in memory: a Tier-2 ('slow') test.
gof_lecessie <- function(ctx, opts = list()) {
  if (!ctx$has_model || is.null(ctx$data))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  fits   <- ctx$ph
  resids <- ctx$y - fits
  N      <- length(resids)
  covs   <- ctx$data[, -1, drop = FALSE]                 # model.frame minus response
  covs   <- covs[, !startsWith(names(covs), "("), drop = FALSE]   # no "(offset)" column
  if (ncol(covs) < 1)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no covariates"))

  # per-covariate squared distances (lower-triangle vectors), summed across covariates
  dlist <- lapply(covs, function(x) {
    if (is.numeric(x)) {
      as.numeric(0.5 * (stats::dist(scale(x)))^2)
    } else {
      xx <- as.numeric(as.factor(x)); nc <- length(unique(xx))
      if (nc <= 1) rep(0, N * (N - 1) / 2)
      else as.numeric((stats::dist(xx, method = "manhattan") != 0) * nc / (nc - 1))
    }
  })
  dist.mat <- matrix(0, N, N)
  dist.mat[lower.tri(dist.mat)] <- sqrt(rowSums(as.data.frame(dlist)))
  dist.mat <- dist.mat + t(dist.mat)
  bandwidth <- if (!is.null(opts$bandwidth)) opts$bandwidth else mean(dist.mat)
  R.raw <- pmax(1 - dist.mat / bandwidth, 0)
  Q.raw <- sum(as.numeric(resids %*% R.raw) * resids)

  X   <- ctx$X
  mu2 <- fits * (1 - fits)
  ## Moment reference for Q = r'Rr with r ~ (I-H)(y-mu):  M = (I-H)' R (I-H), where
  ## H = V X (X'VX)^{-1} X'  (le Cessie and van Houwelingen, 1991, Sec. 6).  The
  ## substitution R -> (I-H)'R(I-H) is prescribed in words in le Cessie and van
  ## Houwelingen (1995), Sec. 4.
  ##
  ## FIX (2026-07-29): this previously computed (I-H) R (I-H), the form inherited from
  ## smwrStats::leCessie.test().  That source computes the weighted, NON-symmetric H
  ## correctly and then omits the transpose, so it departs from the paper it cites: in a
  ## weighted fit H' != H and the sandwich is not the moment matrix of the residual.
  ## Measured cost of the omission (n = 200, 4000 replicates): null size 0.0705 vs
  ## 0.0555 at the 5% level, and 0.0762 vs 0.0515 when the fitted probabilities are
  ## extreme.  It is a size bug, not a power bug.  Do not "simplify" the transpose away.
  ##
  ## M is assembled from the rank-p factors of H instead of forming H (N x N) and
  ## multiplying it out, and the variance trace is evaluated elementwise.  Both steps are
  ## exact algebraic identities: the statistic, the degrees of freedom and the p-value are
  ## unchanged to machine precision.  They drop the cost from O(N^3) to O(N^2 p), which is
  ## what makes this test usable past n ~ 1000.
  A   <- mu2 * X                                 # V X          (N x p)
  Bm  <- solve(crossprod(X, A))                  # (X'VX)^{-1}  (p x p)
  ## H = A Bm X' and H' = X Bm' A', so expand M = R - H'R - RH + H'RH through those
  ## factors.  Every product below is O(N^2 p); none is O(N^3).
  RA   <- R.raw %*% A                            # R V X        (N x p)
  AtR  <- t(RA)                                  # (VX)'R  -- R is symmetric
  XB   <- X %*% Bm
  XBt  <- X %*% t(Bm)
  R.cor <- R.raw - XBt %*% AtR - RA %*% t(XB) + XBt %*% (AtR %*% A) %*% t(XB)
  R.cor <- (R.cor + t(R.cor)) / 2                # enforce the symmetry the algebra gives
  dR    <- diag(R.cor)
  E.Q   <- sum(dR * mu2)
  mu4   <- mu2 * (1 - 3 * mu2)
  VarQ1 <- sum(dR^2 * (mu4 - 3 * mu2^2))
  ## For symmetric M, 2 tr(MVMV) = 2 sum_ij M_ij^2 mu2_i mu2_j.  This elementwise form is
  ## le Cessie and van Houwelingen (1995), eq. (A.6), which states the variance this way
  ## before collapsing it to the trace; evaluating it as written is O(N^2), not O(N^3).
  VarQ2 <- 2 * as.numeric(mu2 %*% (R.cor * R.cor) %*% mu2)
  VarQ  <- VarQ1 + VarQ2
  if (!is.finite(VarQ) || VarQ <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: non-positive variance"))
  Test <- Q.raw * 2 * E.Q / VarQ
  df   <- 2 * E.Q^2 / VarQ
  list(Statistic = Test, df = df, p_value = stats::pchisq(Test, df, lower.tail = FALSE), Note = "")
}

# ---- Tier-2 (opt-in, slow) tests ----

# dplyr::ntile equivalent (equal-sized groups, larger groups first).
.gof_ntile <- function(x, n) {
  r <- rank(x, ties.method = "first")
  as.integer(floor((n * (r - 1) / length(x)) + 1))
}

# Observed-vs-expected chi-square over groups (expected from the tested model's ph).
.gof_oe_chisq <- function(y, ph, grp, dfree) {
  idx <- split(seq_along(y), grp)
  o1 <- vapply(idx, function(I) sum(y[I] == 1), numeric(1))
  o0 <- vapply(idx, function(I) sum(y[I] == 0), numeric(1))
  e1 <- vapply(idx, function(I) sum(ph[I]),     numeric(1))
  e0 <- vapply(idx, function(I) sum(1 - ph[I]), numeric(1))
  keep  <- e1 > 0 & e0 > 0
  chisq <- sum((o1[keep] - e1[keep])^2 / e1[keep] + (o0[keep] - e0[keep])^2 / e0[keep])
  if (dfree <= 0)
    return(list(Statistic = chisq, df = dfree, p_value = NA, Note = "df <= 0"))
  list(Statistic = chisq, df = dfree,
       p_value = stats::pchisq(chisq, dfree, lower.tail = FALSE), Note = "")
}

# Fit the overfit GAM (smooth continuous + main categorical + smooth pair
# interactions) and return its fitted probabilities (used for grouping). From GAM.R.
.gof_gam_pi <- function(ctx) {
  if (!requireNamespace("mgcv", quietly = TRUE) || is.null(ctx$data)) return(NULL)
  dat   <- ctx$data
  resp  <- names(dat)[1]
  preds <- names(dat)[-1]
  preds <- preds[!startsWith(preds, "(")]
  if (length(preds) == 0) return(NULL)
  maxlev <- getOption("ebrahim.gof.pr.maxlev", 6)
  is_cat <- vapply(preds, function(v) {
    col <- dat[[v]]
    is.factor(col) || is.character(col) || is.logical(col) ||
      (is.numeric(col) && length(unique(col)) <= maxlev)
  }, logical(1))
  cont <- preds[!is_cat]; cats <- preds[is_cat]
  terms <- character(0)
  if (length(cont) > 0) terms <- c(terms, paste0("s(", cont, ")"))
  if (length(cats) > 0) terms <- c(terms, cats)
  if (length(cont) > 1)
    for (i in 1:(length(cont) - 1)) for (j in (i + 1):length(cont))
      terms <- c(terms, paste0("s(", cont[i], ",", cont[j], ")"))
  if (length(terms) == 0) return(NULL)
  fml <- stats::as.formula(paste(resp, "~", paste(terms, collapse = "+")))
  g <- tryCatch(suppressWarnings(mgcv::gam(fml, family = stats::binomial(), data = dat)),
                error = function(e) NULL)
  if (is.null(g)) return(NULL)
  list(pi = as.numeric(stats::predict(g, type = "response")),
       cont = cont, cats = cats, k = length(preds))
}

gof_gam_hl <- function(ctx, opts = list()) {
  if (!ctx$has_model) return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  gf <- .gof_gam_pi(ctx)
  if (is.null(gf)) return(list(Statistic = NA, df = NA, p_value = NA, Note = "install 'mgcv' / no covariates"))
  .gof_oe_chisq(ctx$y, ctx$ph, .gof_ntile(gf$pi, 10), 10 - 2)
}

gof_gam_pr <- function(ctx, opts = list()) {
  if (!ctx$has_model) return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  gf <- .gof_gam_pi(ctx)
  if (is.null(gf)) return(list(Statistic = NA, df = NA, p_value = NA, Note = "install 'mgcv' / no covariates"))
  if (length(gf$cats) == 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs a categorical covariate"))
  patt <- do.call(paste, c(lapply(gf$cats, function(v) as.character(ctx$data[[v]])), sep = "_"))
  M <- length(unique(patt)); pig <- gf$pi
  lev <- character(length(ctx$y))
  for (pp in unique(patt)) {
    ix <- which(patt == pp); lev[ix] <- ifelse(pig[ix] <= stats::median(pig[ix]), "low", "high")
  }
  .gof_oe_chisq(ctx$y, ctx$ph, paste(patt, lev, sep = "::"), 2 * M - length(gf$cats) - 2)
}

gof_gam_xie <- function(ctx, opts = list()) {
  if (!ctx$has_model) return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  gf <- .gof_gam_pi(ctx)
  if (is.null(gf)) return(list(Statistic = NA, df = NA, p_value = NA, Note = "install 'mgcv' / no covariates"))
  G  <- if (gf$k < 5) 10 else gf$k + 5
  cm <- as.data.frame(ctx$data[c(gf$cont, gf$cats)])
  for (v in gf$cats) cm[[v]] <- as.numeric(as.factor(cm[[v]]))
  cl <- tryCatch(.gof_kmeans(scale(as.matrix(cm)), G, nstart = 10L), error = function(e) NULL)
  if (is.null(cl)) return(list(Statistic = NA, df = NA, p_value = NA, Note = "clustering failed"))
  .gof_oe_chisq(ctx$y, ctx$ph, cl, round(G - gf$k / 2 - 2))
}

# Stute-Zhu cumulative-residual statistic (residuals ordered by linear predictor).
.gof_tsz_stat <- function(resid, eta) {
  n <- length(resid)
  (1 / n^2) * sum(cumsum(resid[order(eta)])^2)
}

# Stute-Zhu GOF test via parametric (model-based) bootstrap. From Stute-Zhu(Bootstrap).R.
# Sequential by default; when run.all.gof(parallel = TRUE) put a PSOCK cluster on
# ctx$cl, the bootstrap replicates run via parallel::parLapply (reproducible via
# clusterSetRNGStream). opts$B sets the number of bootstrap reps (default 200).
gof_stutezhu <- function(ctx, opts = list()) {
  if (!ctx$has_model) return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  B   <- if (is.null(opts$B)) 200L else as.integer(opts$B)
  y   <- ctx$y; ph <- ctx$ph; X <- ctx$X
  eta <- as.numeric(ctx$model$linear.predictors)
  Tobs <- .gof_tsz_stat(y - ph, eta)
  # One bootstrap replicate, self-contained (base-only closure so PSOCK workers
  # need neither this package's namespace nor any serialized surroundings).
  one_rep <- function(i, probs, Xmat) {
    yb <- stats::rbinom(nrow(Xmat), 1, probs)
    fb <- tryCatch(suppressWarnings(stats::glm.fit(Xmat, yb, family = stats::binomial())),
                   error = function(e) NULL)
    if (is.null(fb) || !isTRUE(fb$converged)) return(NA_real_)
    r <- yb - fb$fitted.values
    e <- as.numeric(Xmat %*% fb$coefficients)
    (1 / length(r)^2) * sum(cumsum(r[order(e)])^2)
  }
  environment(one_rep) <- baseenv()
  Tb <- if (!is.null(ctx$cl)) {
    unlist(parallel::parLapply(ctx$cl, seq_len(B), one_rep, probs = ph, Xmat = X))
  } else {
    vapply(seq_len(B), one_rep, numeric(1), probs = ph, Xmat = X)
  }
  # Monte Carlo p-value (1 + #{T_b >= T}) / (1 + B): never 0, valid for any B
  list(Statistic = Tobs, df = NA_real_,
       p_value = (1 + sum(Tb >= Tobs, na.rm = TRUE)) / (1 + sum(!is.na(Tb))),
       Note = paste0(B, " bootstrap reps",
                     if (!is.null(ctx$cl))
                       paste0(" (parallel, ", length(ctx$cl), " workers)") else ""))
}

# --- eHL: e-value Hosmer-Lemeshow (Henzi et al. 2024) ---
# Adapted (base-R reimplementation, with attribution) from the marius-cp/eHL
# repository. Splits the data, fits isotonic recalibration on the training half,
# and accumulates an e-value on the test half; p = min(1, 1 / mean(e-values)).

.gof_ehl_interp <- function(x_min, x_max, q, xout) {
  x_min[which.min(x_min)] <- 0
  x_max[which.max(x_max)] <- 1
  p  <- c(x_max, x_min)
  qq <- c(q, q)
  o  <- order(p)
  stats::approx(x = p[o], y = qq[o], xout = xout, method = "linear", ties = "ordered")$y
}

.gof_ehl <- function(y, P, boot = 10L, s = 0.5) {
  n  <- length(y)
  ev <- numeric(boot)
  for (t in seq_len(boot)) {
    idx <- sample(n, floor(n * s), replace = FALSE)
    o1  <- order(P[idx]);  P1 <- P[idx][o1];  y1 <- y[idx][o1]      # training, sorted
    o2  <- order(P[-idx]); P2 <- P[-idx][o2]; y2 <- y[-idx][o2]     # test, sorted
    iso <- stats::isoreg(P1, y1)
    iK  <- iso$iKnots
    red <- iK[!duplicated(iso$yf[iK], fromLast = TRUE)]
    bin <- rep.int(seq_along(red), times = diff(c(0, red)))
    sp  <- split(seq_along(P1), bin)
    nb  <- length(sp)
    x_min <- x_max <- q <- numeric(nb)
    for (b in seq_len(nb)) {
      I  <- sp[[b]]; sz <- length(I); sy <- sum(y1[I])
      x_min[b] <- min(P1[I]); x_max[b] <- max(P1[I])
      hatpi    <- sy / sz
      q[b]     <- if (hatpi %in% c(0, 1)) (sy + 0.5) / (sz + 1) else hatpi
    }
    qte   <- .gof_ehl_interp(x_min, x_max, q, P2)
    E     <- (qte^y2 * (1 - qte)^(1 - y2)) / (P2^y2 * (1 - P2)^(1 - y2))
    ev[t] <- prod(E)
  }
  mean(ev)
}

gof_ehl <- function(ctx, opts = list()) {
  boot <- if (is.null(opts$boot)) 10L else as.integer(opts$boot)
  HLe  <- tryCatch(.gof_ehl(ctx$y, ctx$ph, boot = boot, s = 0.5), error = function(e) NA_real_)
  if (!is.finite(HLe))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: e-value computation failed"))
  list(Statistic = HLe, df = NA_real_, p_value = min(1, 1 / HLe),
       Note = "e-value test (reported as p = min(1, 1/e))")
}

# BAGofT (binary-adaptive GOF test). The random-forest
# partitioner needs at least two predictors; for a single-predictor model we add
# a constant helper column to the data (not the formula), the workaround
# documented in Kuss (2002) / the thesis, so the test runs instead of erroring.
# Two engines give the same p-value for the same seed: "fast" (bagoft.fast's core,
# needs randomForest, and dcov above five covariates) and "package" (BAGofT itself).
# The default "auto" takes the fast one when its packages are installed.
gof_bagoft <- function(ctx, opts = list()) {
  if (!ctx$has_model || is.null(ctx$data))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs a fitted glm"))
  nsim  <- if (is.null(opts$nsim)) 100L else as.integer(opts$nsim)
  dat   <- ctx$data
  attr(dat, "terms") <- NULL                       # a model.frame's terms attr breaks BAGofT
  preds <- names(dat)[-1]
  added_const <- FALSE
  if (length(preds) < 2) {                          # single predictor: add a constant column
    dat[[".bagoft_const"]] <- 1
    preds <- c(preds, ".bagoft_const")
    added_const <- TRUE
  }
  engine <- if (is.null(opts$engine)) "auto" else match.arg(opts$engine, c("auto", "fast", "package"))
  fast_ok <- requireNamespace("randomForest", quietly = TRUE) &&
    (length(preds) <= 5L || requireNamespace("dcov", quietly = TRUE))
  if (identical(engine, "auto"))
    engine <- if (fast_ok) "fast" else "package"
  extra <- paste0(if (!is.null(opts$Kmax)) paste0(", Kmax=", opts$Kmax) else "",
                  if (!is.null(opts$nsplits)) paste0(", nsplits=", opts$nsplits) else "")
  if (identical(engine, "fast")) {
    if (!fast_ok)
      return(list(Statistic = NA, df = NA, p_value = NA,
                  Note = "Not run: the fast engine needs 'randomForest' (and 'dcov' above five covariates)"))
    if (!is.null(ctx$model$offset) && any(ctx$model$offset != 0))
      return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: offsets are not supported"))
    r <- tryCatch(suppressWarnings(.bagoft_core(
      ctx$X, ctx$y, dat[-1L], link = ctx$model$family$link,
      nsplits = if (is.null(opts$nsplits)) 100L else as.integer(opts$nsplits),
      nsim = nsim, ne = opts$ne, ntree = if (is.null(opts$ntree)) 60 else opts$ntree,
      Kmax = opts$Kmax, nmin = opts$nmin, mtry = opts$mtry, maxnodes = opts$maxnodes)),
      error = function(e) NULL)
    if (is.null(r) || is.null(r$p.value))
      return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: BAGofT computation failed"))
    return(list(Statistic = NA_real_, df = NA_real_, p_value = .bagoft_row_p(r$p.value, nsim),
                Note = paste0("adaptive RF partition; nsim=", nsim, extra,
                              if (added_const) "; constant column added (single predictor)" else "",
                              "; fast engine", .bagoft_row_note(r$p.value, nsim))))
  }
  if (!requireNamespace("BAGofT", quietly = TRUE))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: install the 'BAGofT' package"))
  link  <- ctx$model$family$link
  # Adaptive random-forest partitioner: pass through any tuning the caller sets
  # via control = list(BAGofT = list(...)). parRF knobs: Kmax (max number of
  # partition cells - the granularity of the adaptive split), ntree, nmin, mtry,
  # maxnodes. BAGofT() knobs: nsim (resampling iterations), nsplits, ne (size of
  # the estimation split).
  rf_args <- list(parVar = preds)
  for (a in c("Kmax", "nmin", "ntree", "mtry", "maxnodes"))
    if (!is.null(opts[[a]])) rf_args[[a]] <- opts[[a]]
  parFun <- do.call(BAGofT::parRF, rf_args)
  bag_args <- list(testModel = BAGofT::testGlmBi(formula = stats::formula(ctx$model), link = link),
                   parFun = parFun, data = dat, nsim = nsim)
  if (!is.null(opts$nsplits)) bag_args$nsplits <- as.integer(opts$nsplits)
  if (!is.null(opts$ne))      bag_args$ne      <- opts$ne
  r <- NULL                                          # silence BAGofT's per-sim console chatter
  invisible(suppressMessages(utils::capture.output(
    r <- tryCatch(suppressWarnings(do.call(BAGofT::BAGofT, bag_args)),
                  error = function(e) NULL))))
  if (is.null(r) || is.null(r$p.value))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: BAGofT computation failed"))
  list(Statistic = NA_real_, df = NA_real_, p_value = .bagoft_row_p(r$p.value, nsim),
       Note = paste0("adaptive RF partition; nsim=", nsim, extra,
                     if (added_const) "; constant column added (single predictor)" else "",
                     "; BAGofT package", .bagoft_row_note(r$p.value, nsim)))
}

# BAGofT's p-value is the share of nsim simulated statistics more extreme than the
# observed one, so it is 0 when none is. A Monte Carlo p-value from nsim draws cannot be
# below 1 / (nsim + 1); the row reports that bound and says so.
.bagoft_row_p <- function(p, nsim) max(as.numeric(p), 1 / (nsim + 1))
.bagoft_row_note <- function(p, nsim)
  if (isTRUE(as.numeric(p) < 1 / (nsim + 1))) "; no simulated statistic as extreme, p set to 1/(nsim+1)" else ""

# McCullagh (1985) exact-conditional-moments standardization of the Pearson
# statistic (SAS GOFLOGIT / Kuss 2002 algorithm, ungrouped binary case). The
# Pearson X^2 is standardized by its conditional mean and variance given the
# fitted coefficients, then Z is referred to the normal. Verified to reproduce
# the thesis low-birth-weight result (p = 0.937) to machine precision.
# The macro forms the n x n hat matrix H = X (X'WX)^{-1} X'; only its diagonal (the
# leverages) and its products with vectors are needed, so they are taken from the
# p x p inverse instead: O(n p^2) time and O(n p) memory, the same numbers to rounding.
# Z is referred to the upper tail, as in GOFLOGIT: a Pearson statistic smaller than
# its conditional mean is not counted as misfit.
.gof_mccullagh_moments <- function(X, P, V, k4) {
  inv <- tryCatch(solve(crossprod(X, V * X)), error = function(e) NULL)
  if (is.null(inv)) return(NULL)
  h    <- rowSums((X %*% inv) * X)                        # diag(H)
  u    <- (1 - 2 * P) / V
  Xu   <- crossprod(X, V * u)                             # X'W u
  uhat <- as.numeric(X %*% (inv %*% Xu))                  # H W u
  k3   <- V * (1 - 2 * P)
  list(Esub = 0.5 * sum(k4 * h / V) - 0.5 * sum(uhat * k3 * h),
       RSSu = sum(V * u^2) - sum(Xu * (inv %*% Xu)))      # u'(W - W H W)u
}

gof_mccullagh <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "Not applicable: needs the design matrix X (pass the fitted glm)"))
  X <- ctx$X; ph <- ctx$ph; y <- ctx$y; n <- length(y); p <- ncol(X)
  V  <- ph * (1 - ph)
  X2 <- sum((y - ph)^2 / V)
  k4 <- ph - 7 * ph^2 + 12 * ph^3 - 6 * ph^4
  mo <- .gof_mccullagh_moments(X, ph, V, k4)
  if (is.null(mo))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: singular information matrix"))
  E   <- (n - p) - mo$Esub
  Var <- (1 - p / n) * mo$RSSu             # ungrouped: 2*sum((m-1)/m) = 0
  if (!is.finite(Var) || Var <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: non-positive variance"))
  Z <- (X2 - E) / sqrt(Var)
  list(Statistic = Z, df = NA_real_, p_value = stats::pnorm(Z, lower.tail = FALSE),
       Note = "normal deviate, upper tail")
}

# ---------------------------------------------------------------------------
# GROUPED (covariate-pattern) variants of the global / standardized tests.
#
# Pearson, deviance and McCullagh are classically defined on the J distinct
# COVARIATE PATTERNS, each a Binomial(m_g, P_g). The default rows above use the
# sparse / one-trial form (m_g = 1), which is what the sparse-data simulation
# needs and is identical to the grouped form when every pattern is unique. On
# data with REPEATED covariate patterns (m_g > 1) the two differ, so these
# "(grouped)" variants are reported alongside the ungrouped ones. For fully
# sparse data they coincide with the rows above.
#
# NOTE on Farrington vs EF: the *original* Farrington (1996) test is a GROUPED
# (covariate-pattern, m_g > 1) test. The Ebrahim-Farrington (EF) test in this
# package is deliberately its SPARSE-data counterpart: it does NOT group by
# covariate pattern but forms G data-dependent bins of the predicted risk, so it
# applies directly to fully sparse (one-observation-per-pattern) data. EF and the
# grouped Farrington therefore answer the same question by different partitions.
# ---------------------------------------------------------------------------
.gof_patterns <- function(X, ph, y) {
  key <- apply(round(X, 10), 1L, paste, collapse = "\r")
  ug  <- !duplicated(key); nm <- key[ug]
  list(m  = as.numeric(tapply(rep(1, length(y)), key, sum)[nm]),
       yg = as.numeric(tapply(y,  key, sum)[nm]),
       Pg = as.numeric(tapply(ph, key, `[`, 1L)[nm]),
       Xu = X[ug, , drop = FALSE], J = sum(ug))
}

gof_pearson_grouped <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X"))
  g <- .gof_patterns(ctx$X, ctx$ph, ctx$y); p <- ncol(ctx$X)
  Vg <- g$m * g$Pg * (1 - g$Pg)
  X2 <- sum((g$yg - g$m * g$Pg)^2 / Vg); df <- g$J - p
  list(Statistic = X2, df = df, p_value = stats::pchisq(X2, df, lower.tail = FALSE),
       Note = sprintf("grouped to %d covariate patterns", g$J))
}

gof_deviance_grouped <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X"))
  g <- .gof_patterns(ctx$X, ctx$ph, ctx$y); p <- ncol(ctx$X)
  ys <- pmax(g$yg, 1e-10); ms <- pmax(g$m - g$yg, 1e-10)
  D  <- 2 * sum(g$yg * log(ys / (g$m * g$Pg)) + (g$m - g$yg) * log(ms / (g$m * (1 - g$Pg))))
  df <- g$J - p
  list(Statistic = D, df = df, p_value = stats::pchisq(D, df, lower.tail = FALSE),
       Note = sprintf("grouped to %d covariate patterns", g$J))
}

gof_mccullagh_grouped <- function(ctx, opts = list()) {
  if (is.null(ctx$X))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not applicable: needs the design matrix X"))
  g <- .gof_patterns(ctx$X, ctx$ph, ctx$y); p <- ncol(ctx$X)
  m <- g$m; P <- g$Pg; yg <- g$yg; Xu <- g$Xu; J <- g$J; Ntot <- sum(m)
  V  <- m * P * (1 - P)
  X2 <- sum((yg - m * P)^2 / V)
  k4 <- m * P * (1 - P) * (1 - 6 * P * (1 - P))
  mo <- .gof_mccullagh_moments(Xu, P, V, k4)
  if (is.null(mo))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: singular information matrix"))
  E    <- (J - p) - mo$Esub
  Var  <- 2 * sum((m - 1) / m) + (1 - p / Ntot) * mo$RSSu   # grouped: 2*sum((m-1)/m) is nonzero
  if (!is.finite(Var) || Var <= 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: non-positive variance"))
  Z <- (X2 - E) / sqrt(Var)
  list(Statistic = Z, df = NA_real_, p_value = stats::pnorm(Z, lower.tail = FALSE),
       Note = sprintf("grouped to %d covariate patterns", J))
}

# GiViTI calibration test (Nattino, Finazzi & Bertolini): forward-selects a
# polynomial calibration model and tests it against the identity, using a
# selection-aware null distribution. Wraps givitiR::givitiCalibrationTest, run in
# an isolated callr subprocess so a crash in givitiR's compiled dependencies
# (alabama, rootSolve) returns NA instead of aborting the session. devel defaults
# to "internal" (the model was fit on these data). Verified against the thesis
# low-birth-weight result (internal p = 0.586).
gof_giviti <- function(ctx, opts = list()) {
  if (!requireNamespace("givitiR", quietly = TRUE))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: install the 'givitiR' package"))
  if (!requireNamespace("callr", quietly = TRUE))
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = "Not run: install the 'callr' package (isolates givitiR)"))
  devel <- if (is.null(opts$devel)) "internal" else opts$devel
  o <- as.numeric(ctx$y); e <- as.numeric(ctx$ph)
  ## starting the subprocess draws from the session's random number stream; give it back,
  ## so that the rows after this one do not depend on whether it ran
  rng <- .rng_snapshot(); on.exit(.rng_restore(rng), add = TRUE)
  pv <- tryCatch(
    callr::r(function(o, e, devel) {
      givitiR::givitiCalibrationTest(o = o, e = e, devel = devel)$p.value
    }, args = list(o = o, e = e, devel = devel), timeout = 120),
    error = function(err) NA_real_)
  if (length(pv) != 1 || !is.finite(pv))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: givitiR calibration test failed"))
  list(Statistic = NA_real_, df = NA_real_, p_value = pv,
       Note = paste0("calibration belt; devel=", devel))
}

# Lai & Liu (2018) standardized-power procedure for the Hosmer-Lemeshow test.
# Reference: Lai, X. & Liu, L. (2018), "A simple test procedure in standardizing
# the power of Hosmer-Lemeshow test in large data sets", J. Statist. Comput.
# Simul. It resamples (stratified by outcome, with replacement) to a target size
# n0, refits the model, computes the HL statistic each time, and estimates the
# rejection rate ("standardized power") at n0. The final decision is randomized:
# reject H0 if U(0,1) < power. There is no p-value: the statistic is the power
# and the accept/reject decision is reported in the Note. From Hosmer Bootstrap.R.
gof_lailiu <- function(ctx, opts = list()) {
  if (!ctx$has_model || is.null(ctx$data))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  if (!requireNamespace("ResourceSelection", quietly = TRUE))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "install 'ResourceSelection'"))
  dat <- ctx$data; attr(dat, "terms") <- NULL
  y   <- dat[[1]]
  n   <- length(y)
  n0    <- if (is.null(opts$n0))    n    else as.integer(opts$n0)   # target sample size
  k     <- if (is.null(opts$k))     200L else as.integer(opts$k)    # resamples
  alpha <- if (is.null(opts$alpha)) 0.05 else opts$alpha
  G     <- ctx$G
  fml   <- stats::formula(ctx$model)
  ev_grp <- dat[y == 1, , drop = FALSE]; ne_grp <- dat[y == 0, , drop = FALSE]
  m <- nrow(ev_grp); nn <- nrow(ne_grp)
  if (m < 2 || nn < 2)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "too few events/non-events"))
  m0  <- max(1, min(round(n0 * (m / n)), m))
  n0e <- max(1, min(n0 - m0, nn))
  crit <- stats::qchisq(1 - alpha, df = G - 2)
  # One resample, self-contained (base + ResourceSelection only) so it can run
  # either sequentially or on the PSOCK cluster run.all.gof(parallel = TRUE)
  # puts on ctx$cl (reproducible via clusterSetRNGStream).
  one_rep <- function(i, ev_grp, ne_grp, m, nn, m0, n0e, fml, G) {
    rs <- rbind(ev_grp[sample(seq_len(m),  m0,  replace = TRUE), , drop = FALSE],
                ne_grp[sample(seq_len(nn), n0e, replace = TRUE), , drop = FALSE])
    mb <- tryCatch(suppressWarnings(stats::glm(fml, data = rs, family = stats::binomial())),
                   error = function(e) NULL)
    if (is.null(mb) || !isTRUE(mb$converged)) return(NA_real_)
    tryCatch(suppressWarnings(as.numeric(
      ResourceSelection::hoslem.test(mb$y, stats::fitted(mb), g = G)$statistic)),
      error = function(e) NA_real_)
  }
  environment(one_rep) <- baseenv()
  hl <- if (!is.null(ctx$cl)) {
    unlist(parallel::parLapply(ctx$cl, seq_len(k), one_rep, ev_grp = ev_grp,
                               ne_grp = ne_grp, m = m, nn = nn, m0 = m0,
                               n0e = n0e, fml = fml, G = G))
  } else {
    vapply(seq_len(k), one_rep, numeric(1), ev_grp = ev_grp, ne_grp = ne_grp,
           m = m, nn = nn, m0 = m0, n0e = n0e, fml = fml, G = G)
  }
  valid <- hl[!is.na(hl)]
  if (length(valid) == 0)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "no valid resamples"))
  power    <- mean(valid >= crit)
  decision <- if (stats::runif(1) < power) "REJECT H0 (lack of fit)" else "fail to reject H0"
  list(Statistic = power, df = NA_real_, p_value = NA_real_,
       Note = paste0("standardized power=", round(power, 3),
                     " (n0=", n0, "); decision: ", decision))
}

# Hosmer-Lemeshow F-test (the "modified H&L" of LogisticDx::gof.glm): sort by
# fitted probability, cut into G equal-frequency groups, and ANOVA-F-test the
# deviance residuals against the group factor. Matches the thesis F_test row.
gof_ftest <- function(ctx, opts = list()) {
  y <- ctx$y; ph <- ctx$ph
  grp <- factor(.gof_groups_ef(ph, ctx$G))
  if (nlevels(grp) < 2)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: fewer than 2 groups"))
  eps <- 1e-10
  dr  <- sign(y - ph) * sqrt(pmax(0, -2 * (y * log(pmax(ph, eps)) +
                                           (1 - y) * log(pmax(1 - ph, eps)))))
  av <- tryCatch(stats::anova(stats::lm(dr ~ grp)), error = function(e) NULL)
  if (is.null(av) || !("grp" %in% rownames(av)))
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "Not run: ANOVA failed"))
  list(Statistic = av["grp", "F value"], df = av["grp", "Df"],
       p_value = av["grp", "Pr(>F)"], Note = "ANOVA F on deviance residuals; comparison only, not counted")
}

# Registry of the bundled tests (test wrappers are internal, not exported).
.GOF_REGISTRY <- list(
  "Pearson"             = list(fn = gof_pearson,          family = "Global",       needs_model = FALSE, slow = FALSE),
  "Pearson (grouped)"   = list(fn = gof_pearson_grouped,  family = "Global",       needs_model = TRUE,  slow = FALSE),
  "Deviance"            = list(fn = gof_deviance,         family = "Global",       needs_model = FALSE, slow = FALSE),
  "Deviance (grouped)"  = list(fn = gof_deviance_grouped, family = "Global",       needs_model = TRUE,  slow = FALSE),
  "Osius-Rojek"         = list(fn = gof_osius,            family = "Standardized", needs_model = TRUE,  slow = FALSE),
  "McCullagh"           = list(fn = gof_mccullagh,        family = "Standardized", needs_model = TRUE,  slow = FALSE),
  "McCullagh (grouped)" = list(fn = gof_mccullagh_grouped, family = "Standardized", needs_model = TRUE, slow = FALSE),
  "Copas-RSS"           = list(fn = gof_copas,            family = "Standardized", needs_model = TRUE,  slow = FALSE),
  "Information-Matrix" = list(fn = gof_im,  family = "Global",       needs_model = TRUE,  slow = FALSE),
  "HL"            = list(fn = gof_hl,       family = "Partition",    needs_model = FALSE, slow = FALSE),
  "HL-largeN"     = list(fn = gof_hl_largeN, family = "Partition",   needs_model = FALSE, slow = FALSE),
  "HL-equalwidth" = list(fn = gof_hlw,      family = "Partition",    needs_model = FALSE, slow = FALSE),
  "Pigeon-Heyse"  = list(fn = gof_ph_test,  family = "Partition",    needs_model = FALSE, slow = FALSE),
  "F-test"        = list(fn = gof_ftest,    family = "Partition",    needs_model = FALSE, slow = FALSE),
  "EF"            = list(fn = gof_ef,       family = "Partition",    needs_model = FALSE, slow = FALSE),
  "EF-normal"     = list(fn = gof_ef_normal, family = "Partition",   needs_model = FALSE, slow = FALSE),
  ## EDGE, the cubic basis, at both partitions of the EDGE paper whatever the battery's G
  "EDGE"          = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "poly3", weights = opts$weights, G = "auto")), family = "Directed", needs_model = TRUE, slow = FALSE),
  "EDGE.G10"      = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "poly3", weights = opts$weights, G = 10)), family = "Directed", needs_model = TRUE, slow = FALSE),
  "DEF.poly2"     = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "poly2",  weights = opts$weights, G = opts$G)), family = "Directed", needs_model = TRUE, slow = FALSE),
  "DEF.poly3"     = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "poly3",  weights = opts$weights, G = opts$G)), family = "Directed", needs_model = TRUE, slow = FALSE),
  "DEF.stukel"    = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "stukel", weights = opts$weights, G = opts$G)), family = "Directed", needs_model = TRUE, slow = FALSE),
  "DEF.sym"       = list(fn = function(ctx, opts) gof_def(ctx, list(basis = "sym",    weights = opts$weights, G = opts$G)), family = "Directed", needs_model = TRUE, slow = FALSE),
  "Stukel"        = list(fn = gof_stukel,   family = "Directed",     needs_model = TRUE,  slow = FALSE),
  "Tsiatis"             = list(fn = gof_tsiatis, family = "Covariate-space", needs_model = TRUE, slow = FALSE),
  "Xie"                 = list(fn = gof_xie,     family = "Covariate-space", needs_model = TRUE, slow = FALSE),
  "Pulkstenis-Robinson" = list(fn = gof_pr,      family = "Covariate-space", needs_model = TRUE, slow = FALSE),
  "le-Cessie"           = list(fn = gof_lecessie, family = "Smoothing",      needs_model = TRUE, slow = TRUE),
  "HL-GAM"              = list(fn = gof_gam_hl,    family = "GAM",            needs_model = TRUE, slow = TRUE),
  "PR-GAM"              = list(fn = gof_gam_pr,    family = "GAM",            needs_model = TRUE, slow = TRUE),
  "Xie-GAM"             = list(fn = gof_gam_xie,   family = "GAM",            needs_model = TRUE, slow = TRUE),
  "Stute-Zhu"           = list(fn = gof_stutezhu,  family = "Bootstrap",      needs_model = TRUE,  slow = TRUE),
  "eHL"                 = list(fn = gof_ehl,        family = "Calibration",    needs_model = FALSE, slow = TRUE),
  "GiViTI"              = list(fn = gof_giviti,      family = "Calibration",    needs_model = FALSE, slow = TRUE),
  "GiViTI-external"     = list(fn = function(ctx, opts) gof_giviti(ctx, list(devel = "external")),
                                                       family = "Calibration",  needs_model = FALSE, slow = TRUE,
                                                       on_request = TRUE),
  "Cubic-LR"            = list(fn = gof_cubic_lr,    family = "Calibration",    needs_model = FALSE, slow = FALSE),
  "BAGofT"              = list(fn = gof_bagoft,     family = "Bootstrap",      needs_model = TRUE,  slow = TRUE),
  "Lai-Liu-HL"          = list(fn = gof_lailiu,     family = "Bootstrap",      needs_model = TRUE,  slow = TRUE),
  "Projection"          = list(fn = gof_proj,       family = "Bootstrap",      needs_model = TRUE,  slow = TRUE)
)
