## localize_gof.R -- error-controlled localization of misfit in a logistic risk model.
##
## Companion software for the localization manuscript. Two user-facing functions:
##   localize.external(y, p, X)  -- frozen predictions on validation data: which of four parts
##                                  of the misfit (INTERCEPT, SLOPE, LINK, COV) is present
##   localize.gof(fit)           -- a fitted logistic glm checked on its own data: LINK and COV
##
## Each group asks one question and has its own Cauchy combination of score tests whose bases are confined to that
## group's part of the misfit (weighted-orthogonal to the earlier parts):
##   INTERCEPT  calibration in the large                      basis 1
##   SLOPE      calibration slope                             basis eta, orthogonal to 1
##   LINK       bends of the score-to-probability map         eta^2, eta^3 | Stukel's two terms | ns(eta, 4),
##                                                            each orthogonal to {1, eta}
##   COV        misfit within covariate patterns              squares and cubes | pairwise products | ns(x_j, 3),
##                                                            each orthogonal to {1, eta, ns(eta, 5)} (functions of eta)
## A group is named when every intersection containing it rejects (closure, Marcus et al. 1976); each intersection
## statistic is the mean Cauchy coordinate of its members' chi-square p-values, calibrated by its rank among the
## reference draws. The numerics are those of the paper's script localize_groups.R: for the same seed the two give
## the same p-values to the last digit, and the tests in test-localize.R hold them to that.

## the order of the groups is the calibration hierarchy; the decision rule acts on the highest one named
.LOC_ORDER <- c("INTERCEPT", "SLOPE", "LINK", "COV")

.LOC_ACTION <- c(
  NONE      = "keep the model; no misfit detected (the verdict is limited by the sample size)",
  INTERCEPT = "update the intercept (overall risk too high or too low)",
  SLOPE     = "logistic recalibration a + b * eta (risks too extreme or too modest)",
  LINK      = "flexible recalibration f(eta) or another link (curved relation between score and risk)",
  COV       = "revise the model; recalibration cannot repair misfit among patients who share a score")

## weighted residualization of the columns of B against Z (w = p(1-p))
.loc_wperp <- function(B, Z, w) {
  if (is.null(Z) || !ncol(Z)) return(B)
  sw <- sqrt(w)
  ## weighted least-squares residuals through a pivoted QR, so redundant columns of Z (e.g. eta inside the span of
  ## an intercept plus a spline in eta) are dropped instead of making the projection singular
  qr.resid(qr(Z * sw, tol = 1e-9), B * sw) / sw
}

## a score-test member: returns a function of the residual matrix R (n x m) giving chi-square p-values (length m).
## The basis is projected with the model weights w and its variance is the model variance sum(w b b').
## Multiplier calibration (r0 = observed residuals given): the variance is robust, sum(e^2 b b'), with e the residuals
## after refitting the projected-off directions Z as an offset model (e = r0 - w Z theta), so that it estimates the
## variance of the empirically projected statistic; R then holds r0 in column 1 and Rademacher signs in the others,
## and the reference draws are sum(b * sign * e). The observed statistic is unchanged, since sum(b w Z) = 0.
.loc_member <- function(B, Z, w, r0 = NULL) {
  Bp <- .loc_wperp(as.matrix(B), Z, w)
  ## drop columns the projection emptied (judged by their norm: the intercept's column is constant but not empty)
  keep <- sqrt(colSums(Bp^2)) > 1e-8 * (1 + sqrt(colSums(as.matrix(B)^2)))
  Bp <- Bp[, keep, drop = FALSE]
  if (!ncol(Bp)) return(NULL)
  ee <- if (is.null(r0)) NULL else if (is.null(Z) || !ncol(Z)) r0 else w * as.numeric(.loc_wperp(matrix(r0 / w), Z, w))
  V  <- if (is.null(ee)) crossprod(Bp, Bp * w) else crossprod(Bp, Bp * ee^2)
  e  <- eigen(V, symmetric = TRUE)
  k  <- sum(e$values > 1e-9 * max(e$values))
  Vi <- e$vectors[, 1:k, drop = FALSE] %*% (t(e$vectors[, 1:k, drop = FALSE]) / e$values[1:k])
  function(R) {
    U <- if (is.null(ee)) crossprod(Bp, R) else cbind(crossprod(Bp, R[, 1]), crossprod(Bp, R[, -1, drop = FALSE] * ee))
    stats::pchisq(colSums(U * (Vi %*% U)), k, lower.tail = FALSE)
  }
}

## the group members on a design: eta and the covariate matrix X (numeric)
.loc_group_bases <- function(eta, X) {
  num <- vapply(seq_len(ncol(X)), function(j) length(unique(X[, j])) > 4, TRUE)
  Xn  <- X[, num, drop = FALSE]
  cs  <- function(v) (v - mean(v)) / stats::sd(v)
  sq  <- if (ncol(Xn)) do.call(cbind, lapply(seq_len(ncol(Xn)), function(j) { z <- cs(Xn[, j]); cbind(z^2, z^3) })) else NULL
  pr  <- if (ncol(X) > 1) do.call(cbind, utils::combn(ncol(X), 2, function(ij) cs(X[, ij[1]]) * cs(X[, ij[2]]), simplify = FALSE)) else NULL
  sp  <- if (ncol(Xn)) do.call(cbind, lapply(seq_len(ncol(Xn)), function(j) splines::ns(Xn[, j], df = 3))) else NULL
  ## with a single covariate every function of it is a function of eta, so there is no covariate-structure part:
  ## the COV group is dropped instead of testing the spline's approximation error
  list(
    LINK = list(poly = cbind(eta^2, eta^3), stukel = cbind(0.5 * eta^2 * (eta >= 0), -0.5 * eta^2 * (eta < 0)),
                spline = splines::ns(eta, df = 4)),
    COV  = if (ncol(X) < 2) list() else Filter(Negate(is.null), list(poly = sq, products = pr, spline = sp)))
}

## the spline dimension used for "functions of eta" when projecting the COV bases: fixed (default 5) or growing with
## n ("auto": the approximation error must vanish faster than n^(-1/2))
.loc_cov_df <- function(cov_df, n) if (identical(cov_df, "auto")) max(5L, as.integer(ceiling(n^(1 / 3)))) else as.integer(cov_df)

## build every member function for one design; Zfit = the model matrix to project off (in-sample) or NULL (external).
## dealias_E (in-sample only, the de-aliased link group): the spline fits on the OBSERVED score of the model columns
## whose conditional mean given the score is curved, frozen, so that the same directions are removed in the data and
## in every bootstrap draw and in-sample LINK reads link shape only
.loc_build_members <- function(eta, X, w, Zfit = NULL, external = TRUE, r0 = NULL, cov_df = 5, dealias_E = NULL) {
  gb <- .loc_group_bases(eta, X)
  one <- matrix(1, length(eta), 1)
  Z2  <- cbind(one, eta, Zfit)
  Z3  <- cbind(Z2, splines::ns(eta, df = .loc_cov_df(cov_df, length(eta))))
  if (!is.null(dealias_E) && ncol(dealias_E)) Z2 <- cbind(Z2, dealias_E)
  mem <- list()
  if (external) {
    mem$INTERCEPT <- list(cal = .loc_member(one, NULL, w, r0))
    mem$SLOPE     <- list(slope = .loc_member(matrix(eta), one, w, r0))
  }
  mem$LINK <- Filter(Negate(is.null), lapply(gb$LINK, function(B) .loc_member(B, Z2, w, r0)))
  mem$COV  <- Filter(Negate(is.null), lapply(gb$COV,  function(B) .loc_member(B, Z3, w, r0)))
  Filter(length, mem)
}

## columns of the model matrix whose weighted mean given the score is curved: weighted F test of a 3-df natural spline
## in eta against a straight line, level .05
.loc_curved_cols <- function(Zfit, eta, w, level = 0.05) {
  n <- length(eta); one <- matrix(1, n, 1)
  L <- cbind(one, eta); S <- cbind(one, splines::ns(eta, df = 3))
  which(vapply(seq_len(ncol(Zfit)), function(j) {
    z <- Zfit[, j]
    r1 <- sum(w * .loc_wperp(matrix(z), L, w)^2); r2 <- sum(w * .loc_wperp(matrix(z), S, w)^2)
    if (r2 <= 0) return(FALSE)
    stats::pf(((r1 - r2) / 2) / (r2 / (n - 4)), 2, n - 4, lower.tail = FALSE) < level
  }, TRUE))
}

## closure over the groups, given member p-values for the observed data (column 1) and the reference draws.
## .cauchy() (legoft.R) is tan((0.5 - p) * pi) with p clipped to [1e-12, 1 - 1e-12], as in the paper's script.
.loc_closure <- function(P, groups, alpha) {
  names_g <- names(groups)
  pint <- list()
  for (r in seq_along(names_g)) for (S in utils::combn(names_g, r, simplify = FALSE)) {
    rows <- unlist(groups[S], use.names = FALSE)
    stat <- colMeans(.cauchy(P[rows, , drop = FALSE]))
    pint[[paste(S, collapse = "+")]] <- (1 + sum(stat[-1] >= stat[1])) / length(stat)
  }
  pint <- unlist(pint)
  has <- function(g) grepl(paste0("(^|\\+)", g, "($|\\+)"), names(pint))
  ## the closed-testing adjusted p-value of a group: the largest p-value of an intersection containing it
  adj <- vapply(names_g, function(g) max(pint[has(g)]), 0)
  named <- vapply(names_g, function(g) all(pint[has(g)] <= alpha), TRUE)
  list(named = names_g[named], intersection = pint, single = pint[names_g], adjusted = adj)
}

## assemble the returned object
.loc_result <- function(cl, P, groups, alpha, setting, calibration, draws, method, data.name, n, dealiased = NULL) {
  top <- if (length(cl$named)) .LOC_ORDER[max(match(cl$named, .LOC_ORDER))] else "NONE"
  structure(list(named = cl$named,
                 highest = if (top == "NONE") NA_character_ else top,
                 action = unname(.LOC_ACTION[[top]]),
                 single = cl$single, adjusted = cl$adjusted, intersection = cl$intersection,
                 members = P[, 1],
                 groups = lapply(groups, function(i) rownames(P)[i]),
                 alpha = alpha, setting = setting, calibration = calibration,
                 draws = draws, n = n, dealiased = dealiased,
                 method = method, data.name = data.name),
            class = "gof_localize")
}

#' Localize the misfit of frozen predictions, with error control
#'
#' Says \emph{which part} of the misfit is present when given probabilities are checked against
#' given 0/1 outcomes: a published risk model on new patients, or any model's predictions on a
#' validation set. Nothing is refitted. The misfit is split into four parts that follow the
#' calibration hierarchy (Van Calster et al. 2016), and each part is named or not with the
#' familywise error rate held at \code{alpha}:
#' \describe{
#'   \item{\code{INTERCEPT}}{calibration in the large: overall risk too high or too low.}
#'   \item{\code{SLOPE}}{the calibration slope: risks too extreme or too modest.}
#'   \item{\code{LINK}}{bends of the map from the score \eqn{\eta = \mathrm{logit}(p)} to risk.}
#'   \item{\code{COV}}{misfit among patients who share a score: a missed interaction or
#'     curved covariate, which no recalibration can repair.}
#' }
#'
#' Each group is a Cauchy combination (Liu and Xie 2020) of score tests whose bases are
#' confined to that group's part, weighted-orthogonal to the parts before it: \code{1};
#' \eqn{\eta}; \eqn{\eta^2, \eta^3}, Stukel's two terms and \code{ns(eta, 4)}; and, for
#' \code{COV}, squares and cubes of the covariates, their pairwise products and
#' \code{ns(x_j, 3)}, each made orthogonal to every function of \eqn{\eta} through
#' \code{ns(eta, cov_df)}. A covariate with four or fewer distinct values enters the products
#' only. With a single covariate every function of it is a function of \eqn{\eta}, so
#' \code{COV} is dropped.
#'
#' A group is named when every intersection of groups containing it rejects: closed testing
#' (Marcus et al. 1976), so the probability of naming any group whose part is absent is at most
#' \code{alpha}. Each intersection is tested by the mean Cauchy coordinate of its members'
#' p-values, referred to its rank among \code{M} reference draws.
#'
#' With \code{calibration = "montecarlo"} (the default) the draws are \eqn{y^* \sim}
#' Bernoulli(\eqn{p}). Nothing is estimated, so under no misfit the p-values are exactly valid
#' at every sample size, and for a part that is absent the error control holds for local
#' departures in the other parts. \code{calibration = "multiplier"} uses a robust variance and
#' sign-flipped residuals, so that an absent part keeps its level for departures of any size in
#' the other parts; that guarantee is asymptotic, and in simulations at \eqn{n \le 4000} this
#' reference was liberal. Use it for study, not yet for decisions.
#'
#' @section Reading a verdict:
#' Act on the highest group named, in the order \code{INTERCEPT < SLOPE < LINK < COV}: it
#' points to the lightest update that repairs the model (Steyerberg et al. 2004).
#' \tabular{ll}{
#'   \code{INTERCEPT} \tab update the intercept \cr
#'   \code{SLOPE} \tab logistic recalibration \eqn{a + b\eta} \cr
#'   \code{LINK} \tab flexible recalibration \eqn{f(\eta)} or another link \cr
#'   \code{COV} \tab revise the model, even if other groups are named too
#' }
#' When nothing is named, no misfit was detected; that verdict is limited by the sample size.
#' When the miscalibration is gross, as after transport to another population, a large error on
#' the logit scale also bends the probability curve, so the named rung can be too high. Then
#' split the validation data at random, reach the verdict and update on one half, and test the
#' updated (frozen) model on the other half with this function.
#'
#' @param y 0/1 outcomes.
#' @param p the predicted probabilities to be checked, one per outcome, in \eqn{[0, 1]}.
#'   Values of exactly 0 or 1 would give an infinite score and are moved to \code{1e-10} and
#'   \code{1 - 1e-10}.
#' @param X a numeric matrix or data frame of covariates, one row per outcome, without missing
#'   values or constant columns. At least one column; \code{COV} needs two.
#' @param M number of reference draws; each p-value lies on a grid of \code{1 / (M + 1)}.
#' @param alpha familywise level.
#' @param calibration \code{"montecarlo"} (the default, exact) or \code{"multiplier"}; see
#'   Details.
#' @param cov_df the degrees of freedom of the spline in \eqn{\eta} that the \code{COV} bases
#'   are made orthogonal to; \code{"auto"} uses \code{max(5, ceiling(n^(1/3)))}.
#' @param seed optional integer, passed to \code{set.seed()} before the reference draws.
#'
#' @return An object of class \code{"gof_localize"}: a list with \code{named} (the groups named,
#'   in hierarchy order), \code{highest} (the highest of them, or \code{NA}), \code{action} (the
#'   update the decision table gives for it), \code{single} (the p-value of each group tested
#'   alone), \code{adjusted} (each group's closed-testing adjusted p-value, the largest p-value of
#'   an intersection containing it: a group is named when it is at most \code{alpha}),
#'   \code{intersection} (the p-value of every intersection, named as in
#'   \code{"SLOPE+COV"}), \code{members} (the chi-square p-value of each member score test),
#'   \code{groups} (the members of each group), \code{alpha}, \code{setting},
#'   \code{calibration}, \code{draws} (the number of reference draws, \code{M}), \code{n},
#'   \code{method} and \code{data.name}.
#'
#' @references
#' Ebrahim, E. K., El-Kotory, A. and Hussein, M. (2026). One goodness-of-fit test is not
#' enough: error-controlled localization of misfit in logistic risk models. Manuscript.
#'
#' Marcus, R., Peritz, E. and Gabriel, K. R. (1976). On closed testing procedures with special
#' reference to ordered analysis of variance. \emph{Biometrika}, 63(3), 655--660.
#' \doi{10.1093/biomet/63.3.655}
#'
#' Liu, Y. and Xie, J. (2020). Cauchy combination test: a powerful test with analytic p-value
#' calculation under arbitrary dependency structures. \emph{Journal of the American Statistical
#' Association}, 115(529), 393--402. \doi{10.1080/01621459.2018.1554485}
#'
#' Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina, M. J. and Steyerberg,
#' E. W. (2016). A calibration hierarchy for risk models was defined: from utopia to
#' empirical data. \emph{Journal of Clinical Epidemiology}, 74, 167--176.
#' \doi{10.1016/j.jclinepi.2015.12.005}
#'
#' Steyerberg, E. W., Borsboom, G. J. J. M., van Houwelingen, H. C., Eijkemans, M. J. C. and
#' Habbema, J. D. F. (2004). Validation and updating of predictive logistic regression models:
#' a study on sample size and shrinkage. \emph{Statistics in Medicine}, 23(16), 2567--2586.
#' \doi{10.1002/sim.1844}
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)              # the published model
#' y <- rbinom(n, 1, plogis(qlogis(p) + 0.8 * X$x1 * X$x2))  # new patients: a missed interaction
#' localize.external(y, p, X, M = 199, seed = 1)   # M = 199 to keep the example fast
#' @seealso \code{\link{localize.gof}} for a fitted model checked on its own data;
#'   \code{\link{deepgof1.external}} and \code{\link{run.all.external}} for one overall test.
#' @concept external validation
#' @concept calibration
#' @concept closed testing
#' @concept familywise error
#' @export
localize.external <- function(y, p, X, M = 999L, alpha = 0.05, calibration = c("montecarlo", "multiplier"),
                              cov_df = 5, seed = NULL) {
  dn <- paste(paste(deparse(substitute(y)), collapse = " "), "and", paste(deparse(substitute(p)), collapse = " "))
  calibration <- match.arg(calibration)
  y <- as.numeric(y); p <- as.numeric(p)
  if (!length(y) || anyNA(y) || !all(y %in% c(0, 1))) stop("'y' must be 0/1, without missing values", call. = FALSE)
  if (length(p) != length(y) || anyNA(p) || any(p < 0 | p > 1))
    stop("'p' must hold one probability in [0, 1] per outcome", call. = FALSE)
  X <- .loc_check_X(X, length(y))
  .loc_check_common(M, alpha, cov_df)
  if (!is.null(seed)) set.seed(seed)

  ## predictions of exactly 0 or 1 would give an infinite score; keep them just inside (0, 1)
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10); eta <- stats::qlogis(p); w <- p * (1 - p)
  r0 <- y - p
  mem <- .loc_build_members(eta, X, w, external = TRUE, r0 = if (calibration == "multiplier") r0 else NULL,
                            cov_df = cov_df)
  R <- if (calibration == "montecarlo") {
    cbind(r0, vapply(seq_len(M), function(m) stats::rbinom(length(p), 1L, p) - p, numeric(length(p))))
  } else {
    ## Rademacher signs, shared by all members so that the draws keep the members' joint dependence
    cbind(r0, matrix(sample(c(-1, 1), length(p) * M, replace = TRUE), length(p), M))
  }
  P <- do.call(rbind, lapply(unlist(mem, recursive = FALSE), function(f) f(R)))
  groups <- lapply(names(mem), function(g) grep(paste0("^", g, "\\."), rownames(P)))
  names(groups) <- names(mem)
  .loc_result(.loc_closure(P, groups, alpha), P, groups, alpha,
              setting = "external validation of frozen predictions",
              calibration = if (calibration == "montecarlo") "Monte Carlo" else "multiplier",
              draws = as.integer(M), n = length(y),
              method = "Localization of misfit: closed testing over orthogonal groups",
              data.name = dn)
}

#' Localize the misfit of a fitted logistic model, with error control
#'
#' The in-sample counterpart of \code{\link{localize.external}}: checks a logistic
#' \code{glm} on the data it was fitted to and names which part of the misfit is present, with
#' the familywise error rate held at \code{alpha}. After a fit with an intercept the
#' \code{INTERCEPT} and \code{SLOPE} parts are identically zero, so two groups remain:
#' \describe{
#'   \item{\code{LINK}}{bends of the map from the linear predictor to risk.}
#'   \item{\code{COV}}{misfit among observations that share a fitted risk: a missed interaction
#'     or curved covariate.}
#' }
#' The bases are those of \code{\link{localize.external}}, also made orthogonal to the model
#' matrix. The reference is the parametric bootstrap: \code{B} outcome vectors are drawn from
#' the fitted model, the model is refitted to each on the same design matrix, and every member
#' p-value is recomputed. A refit that fails is dropped; the number used is returned.
#'
#' In-sample, a model column whose mean given the score is curved can make \code{LINK} read
#' covariate misfit. \code{dealias = TRUE} selects such columns once, on the observed fit (a
#' weighted F test of a 3-df spline in the score against a line, level .05), and makes the
#' \code{LINK} bases orthogonal to their spline fits on the observed score, frozen for every
#' bootstrap draw, so that \code{LINK} reads link shape only. The price is power of \code{LINK}
#' against departures that resemble those columns.
#'
#' Only a logistic fit is served, since the bases and the refits assume the logit link. A
#' model with another link, or any model that is not a \code{glm}, can be checked on
#' validation data with \code{\link{localize.external}}, which needs only its predictions.
#'
#' @section Reading a verdict:
#' Act on the highest group named. \code{COV}: revise the model; the remaining parts are
#' assessed again after the revision. \code{LINK} without \code{COV}: adding a function of the
#' score to the model, or changing the link, is consistent with the data; it is not proof that
#' no revision is needed. Nothing named: no misfit was detected, a verdict limited by the
#' sample size.
#'
#' @param fit a fitted \code{glm} with \code{family = binomial(link = "logit")}, a 0/1 outcome,
#'   one row per observation, no prior weights and no offset.
#' @param X a numeric matrix or data frame of covariates, one row per observation. By default
#'   the numeric columns of the model frame other than the response; a term such as
#'   \code{log(x)} enters as its column \code{log(x)}, while matrix terms (\code{ns()},
#'   \code{poly()}) and factors are left out. \code{COV} needs two columns.
#' @param B number of parametric-bootstrap replicates; each p-value lies on a grid of
#'   \code{1 / (B_used + 1)}.
#' @param alpha familywise level.
#' @param dealias if \code{TRUE}, de-alias the \code{LINK} group from model columns that are
#'   curved in the score; see Details.
#' @param cov_df as in \code{\link{localize.external}}.
#' @param seed optional integer, passed to \code{set.seed()} before the bootstrap.
#'
#' @return An object of class \code{"gof_localize"}, as described in
#'   \code{\link{localize.external}}; \code{draws} is the number of bootstrap replicates used
#'   and \code{dealiased} names the model columns removed from \code{LINK} (if any).
#'
#' @inherit localize.external references
#' @examples
#' set.seed(2)
#' n  <- 500
#' x1 <- rnorm(n); x2 <- rnorm(n)
#' y  <- rbinom(n, 1, plogis(-0.3 + 0.8 * x1 + 0.6 * x2 + 0.8 * x1 * x2))
#' fit <- glm(y ~ x1 + x2, family = binomial())   # the interaction is missed
#' localize.gof(fit, B = 49, seed = 1)   # B = 49 to keep the example fast; use 199 or more
#' @seealso \code{\link{localize.external}}; \code{\link{legoft.localize}} for a two-domain
#'   localization over the classical battery.
#' @concept goodness-of-fit
#' @concept logistic regression
#' @concept closed testing
#' @concept familywise error
#' @export
localize.gof <- function(fit, X = NULL, B = 199L, alpha = 0.05, dealias = FALSE, cov_df = 5, seed = NULL) {
  dn <- paste(deparse(substitute(fit)), collapse = " ")
  if (!inherits(fit, "glm") || !identical(fit$family$family, "binomial") || !identical(fit$family$link, "logit"))
    stop("localize.gof() needs a glm fitted with family = binomial(link = \"logit\"). ",
         "For a model with another link, or any other model, use localize.external() on its ",
         "predicted probabilities for validation data.", call. = FALSE)
  if (any(fit$prior.weights != 1) || !all(fit$y %in% c(0, 1)))
    stop("localize.gof() needs a 0/1 outcome, one row per observation, without prior weights", call. = FALSE)
  if (!is.null(fit$offset) && any(fit$offset != 0))
    stop("localize.gof() does not serve a model with an offset", call. = FALSE)
  if (!is.logical(dealias) || length(dealias) != 1L || is.na(dealias))
    stop("'dealias' must be TRUE or FALSE", call. = FALSE)
  if (is.null(X)) {
    mf <- stats::model.frame(fit)[-1L]
    keep <- vapply(mf, function(v) is.numeric(v) && is.null(dim(v)), TRUE) & !startsWith(names(mf), "(")
    if (!any(keep)) stop("the model has no numeric covariate; give 'X'", call. = FALSE)
    X <- mf[keep]
  }
  X <- .loc_check_X(X, length(fit$y))
  .loc_check_common(B, alpha, cov_df, "B")
  if (!is.null(seed)) set.seed(seed)

  mm <- stats::model.matrix(fit)
  ph0 <- as.numeric(stats::fitted(fit))
  dcols <- if (dealias) .loc_curved_cols(mm[, -1, drop = FALSE], stats::qlogis(ph0), ph0 * (1 - ph0)) else NULL
  dE <- NULL
  if (length(dcols)) {   # fitted once, on the observed score, and frozen
    e0 <- stats::qlogis(ph0); w0 <- ph0 * (1 - ph0)
    Zc <- mm[, -1, drop = FALSE][, dcols, drop = FALSE]
    dE <- Zc - .loc_wperp(Zc, cbind(1, splines::ns(e0, df = 3)), w0)
  }
  members_at <- function(f) {
    ph <- as.numeric(stats::fitted(f)); eta <- stats::qlogis(ph)
    .loc_build_members(eta, X, ph * (1 - ph), Zfit = mm[, -1, drop = FALSE], external = FALSE, cov_df = cov_df,
                       dealias_E = dE)
  }
  pv_of <- function(f, y) { m <- members_at(f); r <- matrix(y - as.numeric(stats::fitted(f)))
                            vapply(unlist(m, recursive = FALSE), function(g) g(r), 0) }
  y0 <- as.numeric(fit$y); ph <- as.numeric(stats::fitted(fit))
  P <- matrix(pv_of(fit, y0), ncol = 1)
  rn <- names(pv_of(fit, y0))
  for (b in seq_len(B)) {
    yb <- stats::rbinom(length(ph), 1L, ph)
    fb <- tryCatch(suppressWarnings(stats::glm.fit(mm, yb, family = stats::binomial())), error = function(e) NULL)
    if (is.null(fb)) next
    class(fb) <- c("glm", "lm"); fb$y <- yb
    P <- cbind(P, tryCatch(pv_of(fb, yb), error = function(e) rep(NA_real_, length(rn))))
  }
  P <- P[, colSums(is.na(P)) == 0, drop = FALSE]
  rownames(P) <- rn
  if (ncol(P) - 1L < 19L) stop("too few usable bootstrap replicates (", ncol(P) - 1L, ")", call. = FALSE)
  gn <- Filter(function(g) any(grepl(paste0("^", g, "\\."), rn)), c("LINK", "COV"))   # COV is absent with one covariate
  groups <- lapply(gn, function(g) grep(paste0("^", g, "\\."), rn)); names(groups) <- gn
  .loc_result(.loc_closure(P, groups, alpha), P, groups, alpha,
              setting = "in-sample checking of a fitted logistic model",
              calibration = "parametric bootstrap", draws = ncol(P) - 1L, n = length(y0),
              method = "Localization of misfit: closed testing over orthogonal groups",
              data.name = dn,
              dealiased = colnames(mm)[-1][dcols])
}

## input checks shared by both functions
.loc_check_X <- function(X, n) {
  if (is.null(X)) stop("'X' is required", call. = FALSE)
  if (is.data.frame(X) && !all(vapply(X, is.numeric, TRUE)))
    stop("'X' must be numeric; code a factor as numeric columns first", call. = FALSE)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("'X' must be numeric", call. = FALSE)
  if (nrow(X) != n) stop("'X' must have one row per outcome (", nrow(X), " rows for ", n, " outcomes)", call. = FALSE)
  if (ncol(X) < 1L) stop("'X' needs at least one column", call. = FALSE)
  if (anyNA(X) || any(!is.finite(X))) stop("'X' has missing or infinite values", call. = FALSE)
  if (any(apply(X, 2, function(v) length(unique(v)) < 2L)))
    stop("'X' has a constant column; drop it", call. = FALSE)
  X
}

.loc_check_common <- function(M, alpha, cov_df, what = "M") {
  if (!is.numeric(M) || length(M) != 1L || is.na(M) || M < 19 || M != round(M))
    stop("'", what, "' must be a whole number of at least 19", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must lie strictly between 0 and 1", call. = FALSE)
  if (!identical(cov_df, "auto") && (!is.numeric(cov_df) || length(cov_df) != 1L || is.na(cov_df) || cov_df < 1))
    stop("'cov_df' must be a positive whole number or \"auto\"", call. = FALSE)
  invisible(TRUE)
}

#' @concept closed testing
#' @export
print.gof_localize <- function(x, digits = 4L, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("data:     ", x$data.name, " (n = ", x$n, ")\n", sep = "")
  cat("setting:  ", x$setting, "\n", sep = "")
  cat(sprintf("reference: %s, %d draws\n", x$calibration, x$draws))
  if (length(x$dealiased))
    cat("LINK de-aliased from: ", paste(x$dealiased, collapse = ", "), "\n", sep = "")
  tab <- data.frame(alone = formatC(x$single, digits = digits, format = "f"),
                    adjusted = formatC(x$adjusted, digits = digits, format = "f"),
                    named = ifelse(names(x$single) %in% x$named, "*", ""),
                    row.names = names(x$single))
  cat("\np-values by group (alone, and closed-testing adjusted):\n")
  print(tab, right = TRUE)
  v <- if (length(x$named)) paste(x$named, collapse = " + ") else "none"
  cat(sprintf("\ngroups named at FWER %s: %s\n", format(x$alpha), v))
  cat("recommended update: ", x$action, "\n\n", sep = "")
  invisible(x)
}
