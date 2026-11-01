## VENDORED into ebrahim.gof 3.0.0 on 2026-11-01 from the cavity paper source
## (PDFs/paper_seriesB_closedform/R/cavity_engine.R and R/cluster_correction.R). Only the shipped design
## is vendored, so the paper and the package compute the same numbers for it. If the paper changes in
## review, re-vendor rather than editing this copy.
# =====================================================================
# cavity.gof() -- the leave-one-out (cavity) Hosmer-Lemeshow test for a penalised logistic fit
#
#   1. fit the ridge, lasso or elastic-net model at the given penalty;
#   2. compute the cavity index xi_i = eta_i - q_i (y_i - pi_i) / (1 - w_i q_i), the one-step
#      leave-one-out linear predictor (active set only for the lasso and the elastic net);
#   3. recalibrate y on xi: logistic (n <= 5000) or logistic-normal (larger n);
#   4. group on the quantiles of xi and form the standardised group residuals e_g;
#   5. refer S = sum e_g^2 to chi-squared on G - 2 degrees of freedom (and to the exact weighted law).
#
# The reference is valid for predictors that behave like independent draws. When many rows of X
# repeat (twins), the observations share the error of their index and the closed-form reference
# rejects too often; the function reports the share of observations with a twin and can refer S to
# its distribution under the fitted model instead (reference = "simulate").
#
# lambda is on the THEORY scale (objective: log-likelihood minus lambda / 2 times the squared norm of
# the slopes), lambda = n * lambda_glmnet for the ridge. Pass lambda_scale = "glmnet" to give glmnet's
# value instead.
# =====================================================================

.cav_expit <- function(z) 1 / (1 + exp(-z))

.cav_ridge <- function(X1, y, lambda, D, beta = rep(0, ncol(X1)), tol = 1e-10, maxit = 100) {
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta)
    pi  <- .cav_expit(eta)
    w   <- pmax(pi * (1 - pi), 1e-8)
    bn  <- drop(solve(crossprod(X1, w * X1) + lambda * D, crossprod(X1, w * (eta + (y - pi) / w))))
    if (max(abs(bn - beta)) < tol) return(bn)
    beta <- bn
  }
  beta
}

# Gauss-Hermite nodes and weights for a standard normal (weights sum to one)
.cav_gh <- function(n) {
  i <- seq_len(n - 1); J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- sqrt(i); J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}

.cav_pvalue <- function(S, Sigma) {
  lam <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  lam <- lam[lam > 1e-9]
  if (!length(lam)) return(NA_real_)
  d <- CompQuadForm::davies(S, lambda = lam, lim = 1e5, acc = 1e-7)
  p <- d$Qq
  if (d$ifault != 0 || is.na(p) || p < 0 || p > 1) p <- CompQuadForm::imhof(S, lambda = lam)$Qq
  min(max(p, 0), 1)
}

# ---- one fit: coefficients, linear predictor, the leave-one-out quadratic form q_i -------------
.cav_fit <- function(X, y, lambda, alpha, beta0 = NULL) {
  n <- nrow(X); X1 <- cbind(1, X)
  if (alpha == 0) {
    D <- diag(c(0, rep(1, ncol(X))))
    beta <- .cav_ridge(X1, y, lambda, D, if (is.null(beta0)) rep(0, ncol(X1)) else beta0)
    E <- seq_len(ncol(X)); XE <- X1; DE <- D; ridge <- lambda
  } else {
    lam_g <- lambda / n
    fit <- glmnet::glmnet(X, y, family = "binomial", alpha = alpha, lambda = lam_g,
                          standardize = FALSE, thresh = 1e-12, maxit = 1e6)
    beta <- as.numeric(stats::coef(fit, s = lam_g))
    E <- which(beta[-1] != 0)
    XE <- cbind(1, X[, E, drop = FALSE]); DE <- diag(c(0, rep(1, length(E))), length(E) + 1)
    ridge <- lambda * (1 - alpha)
  }
  eta <- drop(X1 %*% beta)
  pi  <- pmin(pmax(.cav_expit(eta), 1e-9), 1 - 1e-9)
  w   <- pi * (1 - pi)
  Minv <- solve(crossprod(XE, w * XE) + ridge * DE)
  q <- rowSums((XE %*% Minv) * XE)
  list(beta = beta, eta = eta, pi = pi, w = w, q = q, XE = XE, Minv = Minv, nnz = length(E))
}

.cav_index <- function(f, y) f$eta - f$q / (1 - f$w * f$q) * (y - f$pi)

# ---- recalibration, grouping and the statistic, for a given index ---------------------------------
.cav_statistic <- function(xi, y, G, recal, smooth.degree = 0) {
  n <- length(y)
  if (stats::sd(xi) < 1e-8)
    stop("the cavity index is constant: the fit is degenerate (for example a lasso with no active predictor)")
  rc <- suppressWarnings(stats::glm.fit(cbind(1, xi), y, family = stats::binomial()))
  if (!rc$converged || any(!is.finite(rc$coefficients)) || abs(rc$coefficients[2]) > 50)
    stop("the recalibration of y on the cavity index did not converge or separates the data")
  th <- unname(rc$coefficients); tau <- 0
  if (recal == "logistic") {
    m <- pmin(pmax(.cav_expit(th[1] + th[2] * xi), 1e-9), 1 - 1e-9); v <- m * (1 - m)
    Dm <- v * cbind(1, xi)
  } else {                                     # logistic-normal: E expit(a + b xi + tau Z), Z ~ N(0, 1)
    gh <- .cav_gh(24)
    curve <- function(p) {
      P <- 0
      for (k in seq_along(gh$x)) P <- P + gh$w[k] * .cav_expit(p[1] + p[2] * xi + exp(p[3]) * gh$x[k])
      pmin(pmax(P, 1e-9), 1 - 1e-9)
    }
    nll <- function(p) { m <- curve(p); -sum(y * log(m) + (1 - y) * log(1 - m)) }
    op <- stats::optim(c(th, log(0.5)), nll, method = "L-BFGS-B",
                       lower = c(-20, -20, log(1e-3)), upper = c(20, 20, log(20)))
    th <- op$par[1:2]; tau <- exp(op$par[3])
    m <- curve(op$par); v <- m * (1 - m)
    Dm <- matrix(0, n, 3)
    for (k in seq_along(gh$x)) {
      pk <- .cav_expit(th[1] + th[2] * xi + tau * gh$x[k]); dk <- gh$w[k] * pk * (1 - pk)
      Dm <- Dm + dk * cbind(1, xi, gh$x[k] * tau)
    }
    if (tau <= 1.001e-3) Dm <- Dm[, 1:2, drop = FALSE]   # tau on its bound: the logistic fit
  }
  groups <- pmin(ceiling(rank(xi, ties.method = "random") / (n / G)), G)
  idx <- split(seq_len(n), groups); G <- length(idx)
  Vg  <- vapply(idx, function(I) sum(v[I]), 0.0)
  obs <- vapply(idx, function(I) sum(y[I]), 0.0)
  ex  <- vapply(idx, function(I) sum(m[I]), 0.0)
  e   <- (obs - ex) / sqrt(Vg)
  U   <- t(vapply(idx, function(I) colSums(Dm[I, , drop = FALSE]), numeric(ncol(Dm)))) / sqrt(Vg)
  Om  <- diag(G) - U %*% solve(crossprod(Dm, Dm / v), t(U))
  S   <- sum(e^2)
  out <- list(S = S, e = e, Omega = Om, groups = groups, m = m, theta = c(th, tau = tau),
              table = data.frame(group = seq_len(G), n = lengths(idx), observed = obs, expected = ex,
                                 mean.index = vapply(idx, function(I) mean(xi[I]), 0.0), residual = e))
  if (smooth.degree > 0) {
    xbar <- out$table$mean.index
    Pz <- as.matrix(stats::poly(xbar, smooth.degree))
    zr <- crossprod(Pz, e); out$S.smooth <- as.numeric(t(zr) %*% solve(crossprod(Pz)) %*% zr)
    es <- eigen(Om, symmetric = TRUE); Sh <- es$vectors %*% (sqrt(pmax(es$values, 0)) * t(es$vectors))
    out$p.smooth <- .cav_pvalue(out$S.smooth, Sh %*% (Pz %*% solve(crossprod(Pz), t(Pz))) %*% Sh)
  }
  out
}

# ---- share of observations with a twin: a cosine x_k' M^{-1} x_j / sqrt(q_k q_j) above cut ------
# computed row by row, so no n x n matrix is stored
.cav_twin_share <- function(f, cut) {
  L <- f$XE %*% f$Minv; q <- f$q; n <- length(q); twin <- logical(n)
  for (k in seq_len(n - 1)) {
    j <- (k + 1):n
    hit <- j[drop(f$XE[j, , drop = FALSE] %*% L[k, ]) / sqrt(q[k] * q[j]) > cut]
    if (length(hit)) { twin[k] <- TRUE; twin[hit] <- TRUE }
  }
  mean(twin)
}

#' Leave-one-out (cavity) Hosmer-Lemeshow test for penalised logistic regression
#'
#' Tests the fit of a ridge, lasso or elastic-net logistic regression in high dimensions, including
#' \eqn{p > n}. The textbook Hosmer-Lemeshow test groups on fitted probabilities that are shrunk and
#' have seen their own responses, and after a penalised fit it rejects most correct models. This test
#' groups, recalibrates and standardises on the leave-one-out linear predictor instead, computed from
#' the one fit by a single Newton step, and refers the statistic to the familiar
#' \eqn{\chi^2_{G-2}} law.
#'
#' @param X numeric matrix of predictors (no intercept column). Standardise the predictors before the
#'   call if wanted, not inside a procedure that uses all \eqn{n} observations to fit each one; the fit
#'   itself does not rescale them.
#' @param y binary response (0/1).
#' @param lambda penalty. On the theory scale (log-likelihood minus \eqn{\lambda/2} times the squared
#'   norm of the slopes; \eqn{\lambda = n \lambda_{glmnet}}) unless \code{lambda_scale = "glmnet"}.
#' @param alpha elastic-net mixing parameter: 0 (default) is the ridge, fitted internally; values in
#'   (0, 1] use \pkg{glmnet} (lasso at 1), and the leave-one-out index is computed on the active set.
#' @param G number of groups (default 10).
#' @param lambda_scale \code{"theory"} (default) or \code{"glmnet"}.
#' @param recalibration \code{"auto"} (logistic for \eqn{n \le 5000}, logistic-normal above),
#'   \code{"logistic"} or \code{"logistic-normal"}.
#' @param smooth.degree degree of the smooth version (orthogonal polynomials in the group means of the
#'   index); 0 skips it.
#' @param split if \code{TRUE}, also compute the sample-split test (index from a fit on the other half,
#'   so no approximation and no dependence between indices) and its cross-fitted version, which
#'   combines the two halves by the Cauchy combination.
#' @param reference \code{"closed"} (default) or \code{"simulate"}: refer \eqn{S} also to its
#'   distribution under the fitted model, by simulating responses from the fitted probabilities,
#'   refitting and recomputing \eqn{S}, \code{B} times. Use it when \code{twin.share} is not near zero.
#' @param B number of simulated data sets for \code{reference = "simulate"}.
#' @param check.twins if \code{TRUE} (default), compute the share of observations with a twin. It costs
#'   of order \eqn{n^2} times the size of the active set.
#' @param twin.cut cosine above which two observations count as twins (default 0.999).
#'
#' @return An object of class \code{"cavity.gof"}: a list with \code{statistic} (\eqn{S}), \code{df},
#'   \code{p.value} (\eqn{\chi^2_{G-2}}), \code{p.exact} (the weighted chi-squared law of the grouped
#'   residuals), \code{p.smooth}, \code{p.split}, \code{p.crossfit}, \code{p.simulated} (when
#'   requested), \code{twin.share}, \code{table} (observed and expected counts and standardised
#'   residuals by group), \code{recalibration} (fitted intercept, slope and noise), \code{index} (the
#'   cavity index), \code{beta}, \code{lambda}, \code{alpha}.
#'
#' @details The reference law is derived for predictors that behave like independent draws. When many
#'   rows of \code{X} repeat, observations with the same row share the error of their index and the
#'   closed-form reference rejects too often; in a data set where half of the rows repeated, it rejected
#'   57 per cent of correct models. The function then warns, and \code{reference = "simulate"} gives a
#'   reference conditional on the design.
#'
#' @references Ebrahim, E. K. and Arashi, M. (2026). A leave-one-out goodness-of-fit test for penalised
#'   logistic regression in high dimensions. Submitted. Code archive:
#'   \doi{10.5281/zenodo.23861454}.
#'
#' @examples
#' set.seed(1)
#' n <- 300; p <- 60
#' X <- matrix(rnorm(n * p), n)
#' y <- rbinom(n, 1, 1 / (1 + exp(-drop(X %*% rnorm(p, sd = 0.3)))))
#' res <- cavity.gof(X, y, lambda = 50)
#' res
#' @export
cavity.gof <- function(X, y, lambda, alpha = 0, G = 10, lambda_scale = c("theory", "glmnet"),
                       recalibration = c("auto", "logistic", "logistic-normal"), smooth.degree = 0,
                       split = FALSE, reference = c("closed", "simulate"), B = 199,
                       check.twins = TRUE, twin.cut = 0.999) {
  lambda_scale <- match.arg(lambda_scale); recalibration <- match.arg(recalibration)
  reference <- match.arg(reference)
  X <- as.matrix(X); y <- as.numeric(y); n <- length(y)
  if (!is.numeric(X) || nrow(X) != n) stop("X must be a numeric matrix with one row per response")
  if (anyNA(X) || anyNA(y)) stop("X and y must not contain missing values")
  if (!all(y %in% c(0, 1))) stop("y must be binary (0/1)")
  if (length(lambda) != 1 || !is.finite(lambda) || lambda <= 0) stop("lambda must be one positive number")
  if (alpha < 0 || alpha > 1) stop("alpha must lie in [0, 1]")
  if (alpha > 0 && !requireNamespace("glmnet", quietly = TRUE))
    stop("alpha > 0 needs the package 'glmnet'")
  if (lambda_scale == "glmnet") lambda <- n * lambda
  recal <- if (recalibration == "auto") { if (n <= 5000) "logistic" else "logistic-normal" } else recalibration

  f  <- .cav_fit(X, y, lambda, alpha)
  xi <- .cav_index(f, y)
  st <- .cav_statistic(xi, y, G, recal, smooth.degree)
  G  <- nrow(st$table)
  out <- list(statistic = st$S, df = G - 2, p.value = stats::pchisq(st$S, G - 2, lower.tail = FALSE),
              p.exact = .cav_pvalue(st$S, st$Omega), p.smooth = st$p.smooth, S.smooth = st$S.smooth,
              table = st$table, recalibration = st$theta, recal.family = recal, index = xi,
              beta = f$beta, nnz = f$nnz, lambda = lambda, alpha = alpha, n = n,
              twin.share = NA_real_, p.split = NULL, p.crossfit = NULL, p.simulated = NULL)

  if (check.twins) {
    out$twin.share <- .cav_twin_share(f, twin.cut)
    if (out$twin.share > 0.01 && reference == "closed")
      warning(sprintf(paste0("%.0f%% of the observations have a twin (repeated or nearly repeated rows); ",
                             "the closed-form reference can reject too often. Use reference = \"simulate\"."),
                      100 * out$twin.share), call. = FALSE)
  }

  if (split) {
    A <- sort(sample.int(n, floor(n / 2))); Bh <- setdiff(seq_len(n), A)
    half <- function(train, test) {
      ft <- .cav_fit(X[train, , drop = FALSE], y[train], lambda * length(train) / n, alpha)
      xs <- drop(cbind(1, X[test, , drop = FALSE]) %*% ft$beta)
      stats::pchisq(.cav_statistic(xs, y[test], G, recal)$S, G - 2, lower.tail = FALSE)
    }
    p1 <- half(A, Bh); p2 <- half(Bh, A)
    out$p.split <- p1
    out$p.crossfit <- 0.5 - atan(mean(tan((0.5 - c(p1, p2)) * pi))) / pi
  }

  if (reference == "simulate") {
    Sb <- vapply(seq_len(B), function(b) {
      yb <- stats::rbinom(n, 1, f$pi)
      fb <- .cav_fit(X, yb, lambda, alpha, beta0 = if (alpha == 0) f$beta else NULL)
      tryCatch(.cav_statistic(.cav_index(fb, yb), yb, G, recal)$S, error = function(e) NA_real_)
    }, 0.0)
    out$p.simulated <- (1 + sum(Sb >= st$S, na.rm = TRUE)) / (1 + sum(!is.na(Sb)))
    out$B <- sum(!is.na(Sb))
  }
  class(out) <- "cavity.gof"
  out
}

#' @export
print.cavity.gof <- function(x, ...) {
  cat("\nLeave-one-out (cavity) Hosmer-Lemeshow test\n")
  cat(sprintf("  fit: %s, lambda = %.4g (theory scale), n = %d, active predictors = %d\n",
              if (x$alpha == 0) "ridge" else sprintf("elastic net (alpha = %g)", x$alpha), x$lambda, x$n, x$nnz))
  cat(sprintf("  recalibration: %s, slope = %.3f\n", x$recal.family, x$recalibration[2]))
  cat(sprintf("  S = %.3f on %d df:  P = %.4g  (exact weighted law: P = %.4g)\n",
              x$statistic, x$df, x$p.value, x$p.exact))
  if (!is.null(x$p.smooth)) cat(sprintf("  smooth version: S = %.3f, P = %.4g\n", x$S.smooth, x$p.smooth))
  if (!is.null(x$p.split)) cat(sprintf("  split test: P = %.4g; cross-fitted: P = %.4g\n", x$p.split, x$p.crossfit))
  if (!is.null(x$p.simulated)) cat(sprintf("  simulated reference (B = %d): P = %.4g\n", x$B, x$p.simulated))
  if (!is.na(x$twin.share)) cat(sprintf("  share of observations with a twin: %.3f%s\n", x$twin.share,
                                        if (x$twin.share > 0.01) "  (closed-form reference unreliable)" else ""))
  invisible(x)
}
