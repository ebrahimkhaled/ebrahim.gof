# Ridge precision matrices and ridge networks: the functions of the package
# ridgegraph 0.1.1 by Mohammad Arashi, added here with his agreement. They go with
# Arashi and Ebrahim (2026, submitted), "What does a ridge network estimate?
# Effective penalty and leakage in ridge precision matrix estimation".
#
# The sample covariance is taken as S = X'X / n_eff, where X is the data
# matrix after column centring and n_eff = nrow(X) - 1. After centring, S is
# distributed as the sample covariance of n_eff independent Gaussian vectors,
# so n_eff is the sample size that enters the effective-penalty formulas.

.rg_check_x <- function(x) {
  if (!is.matrix(x) || !is.numeric(x))
    stop("'x' must be a numeric matrix")
  if (anyNA(x)) stop("'x' must not contain missing values")
  if (nrow(x) < 3L || ncol(x) < 2L)
    stop("'x' must have at least 3 rows and 2 columns")
  invisible(TRUE)
}

.rg_spec <- function(x) {
  xc <- scale(x, center = TRUE, scale = FALSE)
  neff <- nrow(x) - 1L
  e <- eigen(crossprod(xc) / neff, symmetric = TRUE)
  list(vec = e$vectors, val = pmax(e$values, 0), neff = neff, xc = xc)
}

.rg_q <- function(sp, lambda) sp$vec %*% (t(sp$vec) / (sp$val + lambda))

.rg_pcor <- function(a) {
  d <- 1 / sqrt(diag(a))
  r <- -a * outer(d, d)
  diag(r) <- 1
  r
}

# 1 - rhat written as (n - r) / n + (lambda / n) sum_{s_j > 0} 1 / (s_j + lambda),
# with r the number of positive eigenvalues, to avoid cancellation for small
# lambda
.rg_shat <- function(sp, lambda, tol = 1e-10) {
  pos <- sp$val[sp$val > tol * max(sp$val)]
  (sp$neff - length(pos)) / sp$neff + lambda * sum(1 / (pos + lambda)) / sp$neff
}

.rg_kappa <- function(sp, lambda) lambda / .rg_shat(sp, lambda)

#' Ridge precision matrix and ridge partial correlations
#'
#' Computes the ridge precision estimator \eqn{(S + \lambda I)^{-1}} from a data
#' matrix and the corresponding partial correlations
#' \eqn{-q_{ij} / (q_{ii} q_{jj})^{1/2}}.
#'
#' The partial correlations do not depend on the scale of the estimator. For
#' Gaussian data they estimate the partial correlations of
#' \eqn{(\Sigma + \kappa I)^{-1}}, where \eqn{\kappa} is the effective penalty
#' returned by \code{\link{effective_penalty}}.
#'
#' @param x numeric data matrix with observations in rows. Columns are centred
#'   internally and \eqn{S = X^\top X / (n - 1)}.
#' @param lambda a single positive penalty.
#' @return A symmetric \eqn{p \times p} matrix. For \code{ridge_pcor} the
#'   diagonal is set to one.
#' @references Arashi, M. and Ebrahim, E. K. (2026). What does a ridge network
#'   estimate? Effective penalty and leakage in ridge precision matrix
#'   estimation. Submitted.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(40 * 10), 40, 10)
#' r <- ridge_pcor(x, lambda = 0.5)
#' range(r[upper.tri(r)])
#' @export
ridge_precision <- function(x, lambda) {
  .rg_check_x(x)
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0)
    stop("'lambda' must be a single positive number")
  .rg_q(.rg_spec(x), lambda)
}

#' @rdname ridge_precision
#' @export
ridge_pcor <- function(x, lambda) {
  .rg_pcor(ridge_precision(x, lambda))
}

#' Effective penalty of the ridge precision estimator
#'
#' Estimates the effective penalty \eqn{\kappa = \lambda (1 + \delta)} of the
#' ridge precision estimator through the degrees-of-freedom ratio
#' \eqn{\hat r = \mathrm{tr}\{S (S + \lambda I)^{-1}\} / (n - 1)}, as
#' \eqn{\hat\kappa = \lambda / (1 - \hat r)}; returns the limit of
#' \eqn{\hat\kappa} as \eqn{\lambda \to 0}; and finds the penalty that attains a
#' prescribed effective penalty. The columns of \code{x} are centred and
#' \eqn{S = X^\top X / (n - 1)}, so the effective sample size in all formulas
#' is \eqn{n - 1}, where \eqn{n} is the number of rows.
#'
#' The estimated effective penalty is continuous and strictly increasing in
#' \eqn{\lambda}, so \code{penalty_for_kappa} has a unique solution.
#'
#' @param x numeric data matrix with observations in rows.
#' @param lambda vector of positive penalties.
#' @param tol relative tolerance below which an eigenvalue of \eqn{S} is treated
#'   as zero.
#' @param kappa target effective penalty, larger than
#'   \code{effective_penalty_floor(x)}.
#' @return \code{effective_penalty} returns a data frame with columns
#'   \code{lambda}, \code{df_ratio} (\eqn{\hat r}), \code{scale}
#'   (\eqn{1 - \hat r}, an estimate of \eqn{1 / (1 + \delta)}) and \code{kappa}.
#'   \code{effective_penalty_floor} returns the harmonic mean of the positive
#'   eigenvalues of \eqn{S} when their number equals the effective sample size,
#'   and zero otherwise. \code{penalty_for_kappa} returns a penalty.
#' @references Arashi, M. and Ebrahim, E. K. (2026). What does a ridge network
#'   estimate? Effective penalty and leakage in ridge precision matrix
#'   estimation. Submitted.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(30 * 60), 30, 60)
#' effective_penalty(x, c(0.01, 0.1, 1))
#' effective_penalty_floor(x)
#' @export
effective_penalty <- function(x, lambda) {
  .rg_check_x(x)
  if (!is.numeric(lambda) || any(lambda <= 0))
    stop("'lambda' must contain positive numbers")
  sp <- .rg_spec(x)
  sc <- vapply(lambda, function(l) .rg_shat(sp, l), numeric(1))
  data.frame(lambda = lambda, df_ratio = 1 - sc, scale = sc,
             kappa = lambda / sc)
}

#' @rdname effective_penalty
#' @export
effective_penalty_floor <- function(x, tol = 1e-10) {
  .rg_check_x(x)
  sp <- .rg_spec(x)
  pos <- sp$val[sp$val > tol * max(sp$val)]
  if (length(pos) < sp$neff) 0 else sp$neff / sum(1 / pos)
}

#' @rdname effective_penalty
#' @importFrom stats uniroot
#' @export
penalty_for_kappa <- function(x, kappa) {
  .rg_check_x(x)
  sp <- .rg_spec(x)
  fl <- effective_penalty_floor(x)
  if (kappa <= fl)
    stop("'kappa' must exceed the limit of the effective penalty, ",
         format(fl, digits = 4))
  uniroot(function(l) .rg_kappa(sp, l) - kappa, c(1e-12, 1e8),
          tol = 1e-12)$root
}

#' Extrapolation of ridge precision estimates across penalties
#'
#' Combines ridge estimates at two effective penalties \eqn{\kappa_1} and
#' \eqn{\kappa_2 = r \kappa_1} with weights \eqn{\kappa_2 / (\kappa_2 -
#' \kappa_1)} and \eqn{-\kappa_1 / (\kappa_2 - \kappa_1)}, after rescaling each
#' estimate by its estimated scale. The population version of this combination
#' differs from the precision matrix by a term of second order in the penalties.
#'
#' The sum of the absolute weights is \eqn{(r + 1) / (r - 1)}, which multiplies
#' the sampling variability.
#'
#' @param x numeric data matrix with observations in rows.
#' @param lambda first (smaller) penalty.
#' @param ratio ratio \eqn{\kappa_2 / \kappa_1 > 1} of the effective penalties.
#' @return A list with the combined matrix \code{precision}, its partial
#'   correlations \code{pcor}, the two penalties \code{lambda}, the two estimated
#'   effective penalties \code{kappa} and the \code{weights}.
#' @references Arashi, M. and Ebrahim, E. K. (2026). What does a ridge network
#'   estimate? Effective penalty and leakage in ridge precision matrix
#'   estimation. Submitted.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(50 * 20), 50, 20)
#' fit <- extrapolated_ridge(x, lambda = 0.5)
#' fit$kappa
#' @importFrom stats uniroot
#' @export
extrapolated_ridge <- function(x, lambda, ratio = 2) {
  .rg_check_x(x)
  if (ratio <= 1) stop("'ratio' must exceed 1")
  sp <- .rg_spec(x)
  k1 <- .rg_kappa(sp, lambda)
  l2 <- uniroot(function(l) .rg_kappa(sp, l) - ratio * k1,
                c(lambda, 1e8), tol = 1e-12)$root
  k2 <- .rg_kappa(sp, l2)
  s1 <- .rg_shat(sp, lambda)
  s2 <- .rg_shat(sp, l2)
  w1 <- k2 / (k2 - k1)
  w2 <- -k1 / (k2 - k1)
  est <- w1 * s1 * .rg_q(sp, lambda) + w2 * s2 * .rg_q(sp, l2)
  structure(list(precision = est, pcor = .rg_pcor(est),
                 lambda = c(lambda, l2), kappa = c(k1, k2),
                 weights = c(w1, w2)), class = "extrapolated_ridge")
}

#' Multi-penalty permutation selection of edges
#'
#' Selects the pairs whose ridge partial correlation exceeds a permutation
#' threshold at some penalty in a grid. For each penalty the threshold is an
#' order statistic of the maximal absolute partial correlation over the observed
#' data and \code{B} data sets in which each column is permuted independently, at
#' level \code{alpha} divided by the number of penalties.
#'
#' If the rows of \code{x} are independent and identically distributed and its
#' columns are mutually independent, the probability that any edge is selected is
#' at most \code{alpha}, for every sample size and dimension. No edge can be
#' selected unless \code{B + 1 >= length(mult) / alpha}.
#'
#' @param x numeric data matrix with observations in rows.
#' @param mult penalty grid in units of \eqn{\mathrm{tr}(S) / p}.
#' @param B number of permutations.
#' @param alpha level in \eqn{(0, 1)}.
#' @return A list with the logical \code{adjacency} matrix, the penalties
#'   \code{lambda} and the thresholds \code{threshold}.
#' @references Arashi, M. and Ebrahim, E. K. (2026). What does a ridge network
#'   estimate? Effective penalty and leakage in ridge precision matrix
#'   estimation. Submitted.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(60 * 8), 60, 8)
#' x[, 2] <- x[, 1] + rnorm(60, sd = 0.5)
#' sel <- mps_select(x, B = 49)
#' which(sel$adjacency & upper.tri(sel$adjacency), arr.ind = TRUE)
#' @export
mps_select <- function(x, mult = c(0.3, 1, 3, 10), B = 99L, alpha = 0.1) {
  .rg_check_x(x)
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must lie in (0, 1)")
  p <- ncol(x)
  ut <- upper.tri(diag(p))
  sp <- .rg_spec(x)
  lams <- mult * sum(sp$val) / p
  robs <- lapply(lams, function(l) .rg_pcor(.rg_q(sp, l)))
  tobs <- vapply(robs, function(r) max(abs(r[ut])), numeric(1))
  tb <- matrix(0, B, length(lams))
  for (b in seq_len(B)) {
    spb <- .rg_spec(apply(x, 2L, sample))
    tb[b, ] <- vapply(lams, function(l) max(abs(.rg_pcor(.rg_q(spb, l))[ut])),
                      numeric(1))
  }
  if (B + 1 < length(lams) / alpha)
    warning("no edge can be selected: B + 1 must be at least ",
            "length(mult) / alpha = ", ceiling(length(lams) / alpha))
  k <- ceiling((1 - alpha / length(lams)) * (B + 1))
  tau <- vapply(seq_along(lams), function(g) sort(c(tobs[g], tb[, g]))[k],
                numeric(1))
  adj <- matrix(FALSE, p, p)
  for (g in seq_along(lams)) adj <- adj | (abs(robs[[g]]) > tau[g])
  diag(adj) <- FALSE
  list(adjacency = adj, lambda = lams, threshold = tau)
}
