## localize_update.R -- which update does a risk model need? External validation with each part tested against the
## model already updated for the simpler parts, and the size of each part on the risk scale.
##   localize.update(y, p, X)  the verdict: INTERCEPT against the frozen model, SLOPE against the intercept update,
##                             LINK against the logistic recalibration, COV against the flexible recalibration
##   localize.size(y, p, X)    how large each part is, in percentage points of risk, with its chance level
## The numerics are those of the paper's scripts localize_update.R and localize_size.R: for the same seed the two give
## the same p-values and sizes, and test-localize-update.R holds them to that.

#' Which update does a risk model need? Localization with the updated-model reference
#'
#' Checks frozen predictions on validation data and names the parts of the misfit that are
#' present, so that the verdict points to the update the model needs. The four parts and their
#' bases are those of \code{\link{localize.external}}; what changes is the reference. Each part
#' is tested against the model already updated for the simpler parts:
#' \tabular{ll}{
#'   \code{INTERCEPT} \tab against the frozen model \eqn{p} (Monte Carlo) \cr
#'   \code{SLOPE} \tab against the intercept update \eqn{a + \eta} \cr
#'   \code{LINK} \tab against the logistic recalibration \eqn{a + b\eta} \cr
#'   \code{COV} \tab against the flexible recalibration \code{ns(eta, 4)}
#' }
#' with \eqn{\eta = \mathrm{logit}(p)}. Except for \code{INTERCEPT} this is a parametric
#' bootstrap: \code{B} outcome vectors are drawn from the fitted update, the update is refitted
#' to each draw, and each part's bases are projected off that update's own directions with its
#' own fitted weights. The hypothesis of each test is therefore "this part is zero, whatever the
#' simpler parts are": a large error in overall risk or slope is absorbed by the update and
#' cannot inflate the tests of the more demanding parts. This makes \code{localize.update()} the
#' function to use after transport to another population, where
#' \code{\link{localize.external}} can name a rung that is too high.
#'
#' With this reference \code{COV} has a fourth member beside those of
#' \code{\link{localize.external}}: the main effects of the covariates, so that a changed
#' coefficient, of a binary covariate too, has a test of its own. Risks pooled by
#' cross-validation are not fully frozen: their coefficients differ from fold to fold, and the
#' main-effects member responds to that difference even when the form of the model is right.
#' Use \code{main = FALSE} for such risks, and check the error rate by repeating the refits in
#' simulated data.
#'
#' The parts are named by Holm's correction of their p-values (the default) or by Bonferroni's.
#'
#' When the validation data come from several centres, outcomes within a centre are correlated
#' and the plain reference names \code{INTERCEPT} too often. Give \code{cluster} and \code{tau},
#' the standard deviation of a centre random effect on the logit scale, usually estimated by a
#' random-intercept model; every reference draw then carries a new effect for each centre.
#'
#' @inheritSection localize.external Reading a verdict
#'
#' @param y 0/1 outcomes.
#' @param p the frozen predicted probabilities, one per outcome, in \eqn{[0, 1]}; values of
#'   exactly 0 or 1 are moved to \code{1e-10} and \code{1 - 1e-10}.
#' @param X a numeric matrix or data frame of covariates, one row per outcome, without missing
#'   values or constant columns. At least one column; \code{COV} needs two.
#' @param B number of reference draws for each part; each p-value lies on a grid of
#'   \code{1 / (B + 1)}.
#' @param alpha familywise level.
#' @param naming \code{"holm"} (the default) or \code{"bonferroni"}.
#' @param cov_df as in \code{\link{localize.external}}.
#' @param cluster optional vector of centre labels, one per outcome.
#' @param tau the standard deviation of the centre effect on the logit scale; used with
#'   \code{cluster}.
#' @param main if \code{TRUE}, \code{COV} includes the main-effects member; see Details.
#' @param seed optional integer, passed to \code{set.seed()} before the reference draws. The
#'   random number stream of the session is restored on exit.
#' @param plot if \code{TRUE}, draw the verdict with \code{\link{plot.gof_localize}}. The
#'   default draws it in an interactive session and not in scripts, simulations or examples.
#'
#' @return An object of class \code{"gof_localize"}, as described in
#'   \code{\link{localize.external}}. Here \code{intersection} holds the Bonferroni test of each
#'   set of parts, whose closure is Holm's procedure; \code{draws} is the smallest number of
#'   usable draws over the parts (a refit that fails is dropped).
#'
#' @references
#' Ebrahim, E. K., Hussein, O. A. E.-A. and El-Kotory, A. (2026). Which update does a clinical
#' risk model need? Localization of misfit in logistic risk models. Submitted.
#' Reproduction materials: \doi{10.5281/zenodo.24022813}
#'
#' Holm, S. (1979). A simple sequentially rejective multiple test procedure.
#' \emph{Scandinavian Journal of Statistics}, 6(2), 65--70.
#'
#' Vergouwe, Y., Nieboer, D., Oostenbrink, R., Debray, T. P. A., Murray, G. D., Kattan, M. W.,
#' Koffijberg, H., Moons, K. G. M. and Steyerberg, E. W. (2017). A closed testing procedure to
#' select an appropriate method for updating prediction models. \emph{Statistics in Medicine},
#' 36(28), 4529--4539. \doi{10.1002/sim.7179}
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)                     # the published model
#' y <- rbinom(n, 1, plogis(-0.7 + qlogis(p) + 0.8 * X$x1 * X$x2))  # new patients
#' localize.update(y, p, X, B = 49, seed = 1)   # B = 49 to keep the example fast; use 199 or more
#' @seealso \code{\link{localize.size}} for the size of each part;
#'   \code{\link{localize.external}} for the frozen-model reference;
#'   \code{\link{localize.gof}} for a fitted model checked on its own data.
#' @concept external validation
#' @concept calibration
#' @concept model updating
#' @concept familywise error
#' @export
localize.update <- function(y, p, X, B = 199L, alpha = 0.05, naming = c("holm", "bonferroni"), cov_df = 5,
                            cluster = NULL, tau = 0, main = TRUE, seed = NULL, plot = interactive()) {
  dn <- paste(paste(deparse(substitute(y)), collapse = " "), "and", paste(deparse(substitute(p)), collapse = " "))
  .loc_check_flag(plot, "plot"); .loc_check_flag(main, "main")
  naming <- match.arg(naming)
  chk <- .loc_check_yp(y, p, X)
  y <- chk$y; X <- chk$X; p <- chk$p; n <- length(y)
  .loc_check_common(B, alpha, cov_df, "B")
  if (!is.numeric(tau) || length(tau) != 1L || is.na(tau) || tau < 0)
    stop("'tau' must be a single number of at least 0", call. = FALSE)
  if (tau > 0 && is.null(cluster)) stop("'tau' needs 'cluster'", call. = FALSE)
  if (!is.null(cluster) && (length(cluster) != n || anyNA(cluster)))
    stop("'cluster' must hold one label per outcome, without missing values", call. = FALSE)
  if (!is.null(seed)) {             # seed the draws here, then give the caller's stream back
    .seed_restore <- .rng_snapshot()
    on.exit(.rng_restore(.seed_restore), add = TRUE)
    set.seed(seed)
  }

  eta0 <- stats::qlogis(p); one <- matrix(1, n, 1)
  ## reference draws; with clusters, each draw adds a new random effect per cluster, N(0, tau^2) on the logit scale
  rdraw <- function(ph) {
    if (is.null(cluster) || tau <= 0) return(stats::rbinom(n, 1L, ph))
    cl <- as.integer(factor(cluster))
    stats::rbinom(n, 1L, stats::plogis(stats::qlogis(ph) + stats::rnorm(max(cl), 0, tau)[cl]))
  }
  cauchy_p <- function(P) { s <- colMeans(.loc_cauchy(P)); (1 + sum(s[-1] >= s[1])) / length(s) }
  gb <- .loc_group_bases(eta0, X)
  if (main && length(gb$COV)) gb$COV <- c(list(main = X), gb$COV)
  N4 <- .loc_ns(eta0, df = 4)
  ## each part: the update it is tested against (fitted to an outcome vector), and its members at that fit
  parts <- list(
    SLOPE = list(fit = function(yy) stats::fitted(stats::glm(yy ~ 1, offset = eta0, family = stats::binomial())),
                 members = function(ph) list(slope = .loc_member(matrix(eta0), one, ph * (1 - ph)))),
    LINK = list(fit = function(yy) stats::fitted(stats::glm(yy ~ eta0, family = stats::binomial())),
                members = function(ph) lapply(gb$LINK, function(Bm) .loc_member(Bm, cbind(one, eta0), ph * (1 - ph)))),
    COV = list(fit = function(yy) stats::fitted(stats::glm(yy ~ N4, family = stats::binomial())),
               members = function(ph) {
                 el <- stats::qlogis(ph)
                 Z <- cbind(one, el, .loc_ns(el, df = .loc_cov_df(cov_df, n)), N4)
                 lapply(gb$COV, function(Bm) .loc_member(Bm, Z, ph * (1 - ph)))
               }))
  if (!length(gb$COV)) parts$COV <- NULL
  ## the member p-values of one part for one outcome vector (the update refitted to it)
  pv <- function(part, yy) {
    ph <- as.numeric(part$fit(yy)); m <- Filter(Negate(is.null), part$members(ph))
    vapply(m, function(f) f(matrix(yy - ph)), 0)
  }

  w <- p * (1 - p); f0 <- .loc_member(one, NULL, w)
  R <- cbind(y - p, vapply(seq_len(B), function(b) rdraw(p) - p, numeric(n)))
  P0 <- matrix(f0(R), nrow = 1)
  single <- c(INTERCEPT = cauchy_p(P0))
  mem <- c(INTERCEPT.cal = P0[1, 1])
  used <- B
  for (g in names(parts)) {
    pr <- parts[[g]]; ph <- as.numeric(pr$fit(y)); p_obs <- pv(pr, y); k <- length(p_obs)
    if (anyNA(p_obs)) stop("the observed test of ", g, " could not be computed", call. = FALSE)
    draws <- vapply(seq_len(B), function(b) {
      yb <- rdraw(ph)
      tryCatch(suppressWarnings(pv(pr, yb)), error = function(e) rep(NA_real_, k))
    }, numeric(k))
    ## keep a members-by-draws matrix also when the part has a single member (vapply then returns a vector)
    P <- cbind(p_obs, matrix(draws, nrow = k))
    P <- P[, colSums(is.na(P)) == 0, drop = FALSE]
    if (ncol(P) - 1L < 19L) stop("too few usable draws for ", g, " (", ncol(P) - 1L, ")", call. = FALSE)
    used <- min(used, ncol(P) - 1L)
    single[g] <- cauchy_p(P)
    mem <- c(mem, stats::setNames(p_obs, paste0(g, ".", names(p_obs))))
  }

  adj <- stats::p.adjust(single, naming)
  ## the Bonferroni test of each set of parts: the tests whose closure is Holm's procedure
  sets <- unlist(lapply(seq_along(single), function(r) utils::combn(names(single), r, simplify = FALSE)), recursive = FALSE)
  pint <- vapply(sets, function(S) min(1, length(S) * min(single[S])), 0)
  names(pint) <- vapply(sets, paste, "", collapse = "+")
  gn <- sub("\\..*$", "", names(mem))
  groups <- lapply(names(single), function(g) which(gn == g)); names(groups) <- names(single)
  out <- .loc_result(list(named = names(single)[adj <= alpha], intersection = pint, single = single, adjusted = adj,
                          naming = naming),
                     matrix(mem, ncol = 1, dimnames = list(names(mem), NULL)), groups, alpha,
                     setting = "external validation of frozen predictions",
                     calibration = if (tau > 0) "updated model with a centre effect" else "updated model",
                     draws = as.integer(used), n = n, robust = FALSE,
                     method = "Localization of misfit: updated-model reference", data.name = dn)
  if (plot) graphics::plot(out)
  out
}

#' How large is each part of a misfit?
#'
#' Gives the size of the four parts of \code{\link{localize.update}} on the risk scale, for
#' frozen predictions on validation data. Four nested updates are fitted: the intercept update,
#' the logistic recalibration, the flexible recalibration \code{ns(eta, 4)}, and a revision that
#' adds the main effects and the \code{COV} bases as one logistic model. The size of a part is
#' the mean absolute change in predicted risk that its update makes beyond the simpler one, in
#' percentage points.
#'
#' An update fitted to noise also moves the predictions, so each size comes with its chance
#' level: the mean size over \code{R} outcome vectors drawn from the simpler update, under which
#' the part is absent. \code{excess} is the size minus the chance level. The noise interval is the
#' size minus the 97.5\% and 2.5\% points of the chance draws; under no misfit it contains 0
#' with probability of about 0.95.
#'
#' The \code{COV} size is a floor: the revision is fitted with noise, and at moderate sample
#' sizes the excess understates the true size. A size that cannot be told from its chance level
#' does not show that a part is absent; \code{\link{localize.update}} carries that evidence.
#'
#' @param y,p,X as in \code{\link{localize.update}}.
#' @param R number of chance draws for each part.
#' @param seed optional integer, passed to \code{set.seed()} before the chance draws. The
#'   random number stream of the session is restored on exit.
#'
#' @return An object of class \code{"gof_localize_size"}: a list with \code{parts} (a data frame
#'   with \code{part}, \code{size}, \code{chance}, \code{excess} and the noise interval
#'   \code{lo}, \code{hi}, in percentage points), \code{calibration} (the calibration intercept
#'   and slope with Wald 95\% intervals) and \code{ici} (the integrated calibration index of
#'   Austin and Steyerberg, in percentage points).
#'
#' @references
#' Ebrahim, E. K., Hussein, O. A. E.-A. and El-Kotory, A. (2026). Which update does a clinical
#' risk model need? Localization of misfit in logistic risk models. Submitted.
#' Reproduction materials: \doi{10.5281/zenodo.24022813}
#'
#' Austin, P. C. and Steyerberg, E. W. (2019). The Integrated Calibration Index (ICI) and
#' related metrics for quantifying the calibration of logistic regression models.
#' \emph{Statistics in Medicine}, 38(21), 4051--4065. \doi{10.1002/sim.8281}
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)
#' y <- rbinom(n, 1, plogis(-0.7 + qlogis(p)))   # overall risk too high
#' localize.size(y, p, X, R = 20, seed = 1)      # R = 20 to keep the example fast; use 200
#' @seealso \code{\link{localize.update}} for the verdict.
#' @concept external validation
#' @concept calibration
#' @concept model updating
#' @export
localize.size <- function(y, p, X, R = 200L, seed = NULL) {
  chk <- .loc_check_yp(y, p, X)
  y <- chk$y; X <- chk$X; p <- chk$p; n <- length(y)
  if (!is.numeric(R) || length(R) != 1L || is.na(R) || R < 2 || R != round(R))
    stop("'R' must be a whole number of at least 2", call. = FALSE)
  if (!is.null(seed)) {
    .seed_restore <- .rng_snapshot()
    on.exit(.rng_restore(.seed_restore), add = TRUE)
    set.seed(seed)
  }
  eta <- stats::qlogis(p)
  C <- cbind(X, do.call(cbind, .loc_group_bases(eta, X)$COV))
  S4 <- .loc_ns(eta, df = 4)
  ## the update of rung k fitted to an outcome vector; rung 0 is the frozen model
  rung <- function(k, yy) switch(k + 1L,
    p,
    stats::fitted(stats::glm(yy ~ 1, offset = eta, family = stats::binomial())),
    stats::fitted(stats::glm(yy ~ eta, family = stats::binomial())),
    stats::fitted(stats::glm(yy ~ S4, family = stats::binomial())),
    stats::fitted(suppressWarnings(stats::glm(yy ~ S4 + C, family = stats::binomial()))))
  size <- function(k, yy) 100 * mean(abs(rung(k, yy) - rung(k - 1L, yy)))
  obs <- vapply(1:4, size, 0, yy = y)
  ## chance level of part k: outcomes drawn from the simpler update, under which part k is absent
  draws <- vapply(1:4, function(k) {
    base <- rung(k - 1L, y)
    replicate(R, suppressWarnings(size(k, stats::rbinom(n, 1L, base))))
  }, numeric(R))
  chance <- colMeans(draws); q <- apply(draws, 2, stats::quantile, c(.025, .975))
  parts <- data.frame(part = .LOC_ORDER, size = obs, chance = chance, excess = obs - chance,
                      lo = obs - q[2, ], hi = obs - q[1, ])
  wald <- function(f, j) { e <- stats::coef(f)[j]; s <- sqrt(diag(stats::vcov(f)))[j]; c(e, e - 1.96 * s, e + 1.96 * s) }
  cal <- rbind(intercept = wald(stats::glm(y ~ 1, offset = eta, family = stats::binomial()), 1),
               slope = wald(stats::glm(y ~ eta, family = stats::binomial()), 2))
  dimnames(cal) <- list(c("intercept", "slope"), c("estimate", "lower", "upper"))
  lo <- stats::loess(y ~ p, span = 0.75, degree = 2)
  structure(list(parts = parts, calibration = cal, ici = 100 * mean(abs(stats::predict(lo) - p))),
            class = "gof_localize_size")
}

#' @concept calibration
#' @export
print.gof_localize_size <- function(x, digits = 2L, ...) {
  cat("\nSize of each part of the misfit, percentage points of risk\n\n")
  tab <- x$parts[-1]
  tab[] <- lapply(tab, function(v) formatC(v, digits = digits, format = "f"))
  names(tab) <- c("size", "chance", "excess", "noise.lo", "noise.hi")
  rownames(tab) <- x$parts$part
  print(tab, right = TRUE)
  cat(sprintf("\ncalibration intercept %.3f (%.3f to %.3f), slope %.3f (%.3f to %.3f)\n",
              x$calibration[1, 1], x$calibration[1, 2], x$calibration[1, 3],
              x$calibration[2, 1], x$calibration[2, 2], x$calibration[2, 3]))
  cat(sprintf("integrated calibration index %.2f\n\n", x$ici))
  invisible(x)
}

## input checks shared by localize.update() and localize.size()
.loc_check_yp <- function(y, p, X) {
  y <- as.numeric(y); p <- as.numeric(p)
  if (!length(y) || anyNA(y) || !all(y %in% c(0, 1))) stop("'y' must be 0/1, without missing values", call. = FALSE)
  if (length(p) != length(y) || anyNA(p) || any(p < 0 | p > 1))
    stop("'p' must hold one probability in [0, 1] per outcome", call. = FALSE)
  list(y = y, p = pmin(pmax(p, 1e-10), 1 - 1e-10), X = .loc_check_X(X, length(y)))
}
