## VENDORED into ebrahim.gof 2.7.0 on 2026-09-05 from the Series B paper source
## (PDFs/paper_seriesB_closedform/R/closedform_engine.R). The research engine there carries
## diagnostic paths -- alternative references, alternative groupings, alternative inflations --
## that exist to justify the shipped design; only the shipped design is vendored here, so the
## paper and the package compute the same numbers for it. If the paper changes in review,
## re-vendor rather than editing this copy.
# =====================================================================
# calm.gof() -- the closed-form (no-bootstrap) reference for the
#               shrinkage-corrected grouped goodness-of-fit statistics
#
#   SC.HL             the decile (Hosmer-Lemeshow) basis
#   SC.EDGE           the EDGE basis of degree three
#   SC.EDGE.adaptive  the EDGE basis whose degree is chosen from the data
#
# Companion to shrink.gof(), which refers the same statistics to a bootstrap.
# CALM = Calibration Assessment under Lambda-shrunk Models.
#
# The procedure:
#   1. fit the ridge (or other smooth-penalty) model by penalized IRLS;
#   2. form the grouped standardized residuals and subtract the estimated
#      shrinkage displacement, giving the corrected residual v;
#   3. build the null covariance of v from the DE-NOISED Bernoulli variance,
#      estimated from the fit alone through the observable adjustments of
#      Bellec (2025);
#   4. inflate that covariance by a factor of order p/n for the selection
#      effect of grouping on a fitted index, and read the exact tail of the
#      resulting weighted chi-squared law (Davies 1980).
#
# lambda is on the THEORY scale: lambda = n * lambda_glmnet. Pass
# lambda_scale = "glmnet" to give glmnet's value instead; both are printed.
# =====================================================================

.calm_expit <- function(z) 1 / (1 + exp(-z))

.calm_ridge <- function(X1, y, lambda, D, tol = 1e-10, maxit = 100) {
  beta <- rep(0, ncol(X1))
  for (it in seq_len(maxit)) {
    eta <- drop(X1 %*% beta)
    pi  <- .calm_expit(eta)
    w   <- pmax(pi * (1 - pi), 1e-8)
    z   <- eta + (y - pi) / w
    bn  <- drop(solve(crossprod(X1, w * X1) + lambda * D, crossprod(X1, w * z)))
    if (max(abs(bn - beta)) < tol) { beta <- bn; break }
    beta <- bn
  }
  beta
}

.calm_qf <- function(u, Z) {
  zr <- crossprod(Z, u)
  as.numeric(t(zr) %*% solve(crossprod(Z)) %*% zr)
}

# ---- everything the statistics need, from one fit ---------------------
.calm_pieces <- function(X, y, lambda, G) {
  X <- as.matrix(X); y <- as.numeric(y); n <- length(y)
  X1 <- cbind(1, X); D <- diag(c(0, rep(1, ncol(X))))
  beta <- .calm_ridge(X1, y, lambda, D)
  eta  <- drop(X1 %*% beta)
  pi   <- pmin(pmax(.calm_expit(eta), 1e-6), 1 - 1e-6)
  w    <- pmax(pi * (1 - pi), 1e-8)
  groups <- pmin(ceiling(rank(pi, ties.method = "first") / (n / G)), G)
  idx <- split(seq_len(n), groups); G <- length(idx)
  Vg  <- vapply(idx, function(I) sum(w[I]), 0.0)
  r   <- (vapply(idx, function(I) sum(y[I]), 0.0) -
          vapply(idx, function(I) sum(pi[I]), 0.0)) / sqrt(Vg)
  U   <- t(vapply(idx, function(I) colSums(w[I] * X1[I, , drop = FALSE]),
                  numeric(ncol(X1)))) / sqrt(Vg)
  Fm <- crossprod(X1, w * X1); K <- lambda * D; M <- Fm + K
  Fi <- solve(Fm); Mi <- solve(M)
  bt <- beta + drop(Fi %*% (K %*% beta))      # one-step debiased estimate
  mu <- drop(U %*% Mi %*% (K %*% bt))         # the shrinkage displacement
  v  <- r - mu                                # the corrected residual
  pbar <- vapply(idx, function(I) mean(pi[I]), 0.0)
  Z <- tryCatch(as.matrix(stats::poly(pbar, 3)), error = function(e) NULL)
  XW2X <- crossprod(X1, w^2 * X1)
  df_M <- sum(diag(Mi %*% Fm)); trV_M <- sum(w) - sum(diag(Mi %*% XW2X))
  list(v = v, r = r, mu = mu, Z = Z, groups = groups, pi = pi, beta = beta,
       S_dec = sum(v^2), S_edge = if (is.null(Z)) NA_real_ else .calm_qf(v, Z),
       Om_MLE = diag(G) - U %*% Fi %*% t(U),
       df_M = df_M, trV_M = trV_M, gamma_hat = df_M / trV_M, v_hat = trV_M / n)
}

# ---- exact tail probability of a Gaussian quadratic form --------------
.calm_pvalue <- function(S, Sigma, Z = NULL) {
  if (is.na(S)) return(NA_real_)
  if (is.null(Z)) {
    lam <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  } else {
    e  <- eigen(Sigma, symmetric = TRUE)
    Sh <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
    P  <- Z %*% solve(crossprod(Z), t(Z))
    lam <- eigen(Sh %*% P %*% Sh, symmetric = TRUE, only.values = TRUE)$values
  }
  lam <- lam[lam > 1e-9]
  if (!length(lam)) return(NA_real_)
  d <- CompQuadForm::davies(S, lambda = lam, lim = 1e5, acc = 1e-7)
  p <- d$Qq
  if (d$ifault != 0 || is.na(p) || p < 0 || p > 1)
    p <- CompQuadForm::imhof(S, lambda = lam)$Qq
  min(max(p, 0), 1)
}

# ---- observable adjustments (Bellec 2025, eq. 3.20) -------------------
# Bellec's model has no intercept, so the design and the fitted index are centred first;
# Sigma_x is unknown and replaced by the sample covariance of the centred design.
.calm_bellec <- function(X, y, pc) {
  X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  Xc  <- scale(X, center = TRUE, scale = FALSE)
  eta <- drop(cbind(1, X) %*% pc$beta); psi <- y - pc$pi
  Xb  <- eta - mean(eta)
  r2  <- sum(psi^2) / n
  S   <- crossprod(Xc) / n
  eS  <- eigen(S, symmetric = TRUE)
  Sih <- eS$vectors %*% (1 / sqrt(pmax(eS$values, 1e-8)) * t(eS$vectors))
  t2 <- sum((Sih %*% crossprod(Xc, psi))^2) / n^2 +
        (2 * pc$v_hat / n) * sum(psi * Xb) +
        (pc$v_hat^2 / n) * sum((Xb - pc$gamma_hat * psi)^2) - (p / n) * r2
  a2 <- ((pc$v_hat / n) * sum((Xb - pc$gamma_hat * psi)^2) + sum(psi * Xb) / n -
         pc$gamma_hat * r2)^2 / max(t2, 1e-8)
  s2 <- max(sum((Xb - pc$gamma_hat * psi)^2) / n - a2, 1e-8)
  list(m = Xb - pc$gamma_hat * psi, a2 = a2, s2 = s2, mean_eta = mean(eta))
}

# Gauss-Hermite nodes and weights for a standard normal, by Golub-Welsch (no dependency)
.calm_gh <- function(n) {
  i <- seq_len(n - 1); J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- sqrt(i); J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}

# ---- the de-noised null Bernoulli variance ---------------------------
.calm_w0 <- function(X, y, pc, nodes = 20) {
  b   <- .calm_bellec(X, y, pc)
  k   <- b$a2 / (b$a2 + b$s2)
  mu  <- b$m * k
  tau <- sqrt(b$a2 * b$s2 / (b$a2 + b$s2))
  gh  <- .calm_gh(nodes)
  nll <- function(par) {
    alpha <- par[1]; cc <- exp(par[2]); ll <- 0
    for (q in seq_len(nodes)) {
      z  <- mu + tau * gh$x[q]
      pr <- pmin(pmax(.calm_expit(alpha + cc * z), 1e-9), 1 - 1e-9)
      ll <- ll + gh$w[q] * ifelse(y == 1, pr, 1 - pr)
    }
    -sum(log(pmax(ll, 1e-300)))
  }
  # log c is boxed: an unconstrained search diverges on the occasional sample towards a
  # step-function link, which is not a scale the logistic model can identify.
  op <- stats::optim(c(b$mean_eta, 0), nll, method = "L-BFGS-B",
                     lower = c(-10, log(0.05)), upper = c(10, log(20)))
  alpha <- op$par[1]; cc <- exp(op$par[2])
  w0 <- numeric(length(y))
  for (q in seq_len(nodes)) {
    z <- mu + tau * gh$x[q]; pr <- .calm_expit(alpha + cc * z)
    w0 <- w0 + gh$w[q] * pr * (1 - pr)
  }
  list(w0 = pmax(w0, 1e-8), c = cc, alpha = alpha,
       est_signal_sd = cc * sqrt(b$a2), a2 = b$a2, s2 = b$s2)
}

.calm_reference <- function(X, y, lambda, pc) {
  X1 <- cbind(1, as.matrix(X)); n <- nrow(X1); q <- ncol(X1)
  idx <- split(seq_len(n), pc$groups); G <- length(idx)
  w_h <- pmax(pc$pi * (1 - pc$pi), 1e-8)
  Vg  <- vapply(idx, function(I) sum(w_h[I]), 0.0)
  Cmat <- matrix(0, G, n); for (g in seq_len(G)) Cmat[g, idx[[g]]] <- 1
  U  <- t(vapply(idx, function(I) colSums(w_h[I] * X1[I, , drop = FALSE]),
                 numeric(q))) / sqrt(Vg)
  Fm <- crossprod(X1, w_h * X1)
  Astar <- Cmat / sqrt(Vg) - U %*% solve(Fm, t(X1))
  dn <- .calm_w0(X, y, pc)
  list(Sigma = Astar %*% (dn$w0 * t(Astar)), c = dn$c,
       est_signal_sd = dn$est_signal_sd, a2 = dn$a2, s2 = dn$s2)
}

.calm_edge_basis <- function(pc, degree) {
  idx <- split(seq_along(pc$pi), pc$groups)
  pbar <- vapply(idx, function(I) mean(pc$pi[I]), 0.0)
  tryCatch(as.matrix(stats::poly(pbar, degree)), error = function(e) NULL)
}

#' Closed-form goodness-of-fit test for penalized logistic regression (CALM)
#'
#' Refers the shrinkage-corrected grouped goodness-of-fit statistics to an
#' analytic reference distribution, so that no bootstrap is needed. It is the
#' companion of \code{\link{shrink.gof}}, which refers the same statistics to a
#' prepivoting bootstrap: the statistics are identical, only the reference
#' differs. CALM stands for Calibration Assessment under Lambda-shrunk Models.
#'
#' Under a ridge penalty the grouped standardized residuals are displaced by
#' shrinkage. Subtracting an estimate of that displacement restores the
#' maximum likelihood covariance to first order, but the resulting reference is
#' conservative, because it standardizes by the Bernoulli variance at the
#' \emph{fitted} probabilities, which shrinkage inflates towards one quarter.
#' CALM replaces it by a de-noised estimate of the null Bernoulli variance,
#' obtained from the fit alone through the observable adjustments of Bellec
#' (2025), and inflates the result for the effect of grouping on an index that
#' depends on the response.
#'
#' Two bases are returned. The EDGE basis projects the corrected residual onto
#' orthogonal polynomials in the group-mean fitted probability, and
#' \code{SC.EDGE.adaptive} chooses the degree from the data: it keeps the cubic
#' direction when the aspect ratio is small, or when the observable index
#' correlation \eqn{\hat\rho} satisfies \eqn{\hat\rho^{6} \ge \tau}, and uses
#' degree two otherwise. That basis is the one to prefer. The decile basis
#' \code{SC.HL} is the shrinkage-corrected Hosmer-Lemeshow statistic and is
#' reported for continuity with that tradition; because it spends every group
#' direction it cannot avoid the direction the selection effect occupies, and
#' its reference relies on a constant calibrated by simulation, which does not
#' transfer to every design. See the reference for the designs in which it
#' fails.
#'
#' @param X numeric matrix or data frame of predictors, without an intercept
#'   column. Standardize the columns as you would before any ridge fit.
#' @param y numeric or integer vector of 0/1 responses.
#' @param lambda the ridge penalty. On the theory scale by default, that is on
#'   the scale of the log-likelihood; pass \code{lambda_scale = "glmnet"} to
#'   give \code{glmnet}'s value instead, which is \code{lambda / n}.
#' @param G number of equal-frequency groups. Default 10.
#' @param basis which statistics to compute: any of \code{"decile"},
#'   \code{"adaptive"} and \code{"edge"}. Default: all three.
#' @param lambda_scale \code{"theory"} (default) or \code{"glmnet"}.
#' @param tau threshold of the degree rule, in (0, 1). Default 0.2. Values
#'   below 0.2 keep the cubic direction too often in the proportional regime.
#' @param inflate whether to inflate the reference for the selection effect,
#'   \code{"kappa"} (default) or \code{"none"}. On the adaptive basis the
#'   inflation changes the size by at most 0.008 in the designs examined.
#'
#' @return An object of class \code{"calm.gof"}: a list with components
#'   \code{SC.HL}, \code{SC.EDGE} and \code{SC.EDGE.adaptive} (each a list with
#'   \code{statistic} and \code{p.value}, and for the adaptive basis also the
#'   chosen \code{degree} and the observable \code{rho_hat}), together with the
#'   penalty on both scales and the observables used to build the reference.
#'
#' @section Scope:
#'   The reference is validated for aspect ratios \eqn{p/n} up to 0.4 and for
#'   penalties that shrink towards zero; the lasso is not covered, because the
#'   displacement requires a differentiable penalty. The test asks whether the
#'   logistic form is correct along the fitted index. It does not ask whether
#'   the shrunk probabilities are calibrated, which under a penalty they are
#'   not, by an amount the analyst chose when selecting \code{lambda}.
#'
#'   The reference is validated for outcomes that are not strongly unbalanced.
#'   Once \eqn{p/n} is an appreciable fraction the level is lost as the
#'   prevalence falls: at \eqn{p/n = 0.25} the smooth basis rejects 0.140,
#'   0.574 and 0.884 of correctly specified models at prevalences 0.30, 0.15
#'   and 0.08, and the decile basis 0.060, 0.204 and 0.492. At fixed dimension
#'   the decile basis is unaffected and the smooth one degrades far more
#'   slowly, to 0.130 at prevalence 0.08. What fails there is the
#'   response-dependent grouping inherited from the Hosmer-Lemeshow
#'   construction rather than the reference itself. With few events and
#'   \eqn{p/n} an appreciable fraction, neither statistic is validated and
#'   \code{\link{shrink.gof}} is the less badly behaved of the two.
#'
#' @references
#' Bellec, P. C. (2025). Observable adjustments in single-index models for
#' regularized M-estimators with bounded p/n. \emph{The Annals of Statistics},
#' \strong{53}(2), 531--560.
#'
#' Davies, R. B. (1980). Algorithm AS 155: the distribution of a linear
#' combination of chi-squared random variables. \emph{Journal of the Royal
#' Statistical Society, Series C}, \strong{29}(3), 323--333.
#'
#' Ebrahim, E. K. (2026). A closed-form reference distribution for
#' goodness-of-fit testing under penalized logistic regression in the
#' proportional regime.
#'
#' @seealso \code{\link{shrink.gof}} for the bootstrap reference for the same
#'   statistics, and \code{\link{edge.gof}} for the EDGE test on an unpenalized fit.
#'
#' @examples
#' set.seed(1)
#' n <- 300; p <- 20
#' X <- matrix(rnorm(n * p), n)
#' y <- rbinom(n, 1, 1 / (1 + exp(-(X[, 1] - 0.5 * X[, 2]))))
#' calm.gof(X, y, lambda = 100)
#'
#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept CALM
#' @concept penalized regression
#' @concept ridge regression
#' @concept high-dimensional
#' @concept Hosmer-Lemeshow
#' @export
calm.gof <- function(X, y, lambda, G = 10,
                     basis = c("decile", "adaptive", "edge"),
                     lambda_scale = c("theory", "glmnet"),
                     tau = 0.2,
                     inflate = c("kappa", "none")) {
  if (!requireNamespace("CompQuadForm", quietly = TRUE))
    stop("calm.gof() needs the 'CompQuadForm' package for the exact tail probability; ",
         "install it with install.packages(\"CompQuadForm\").", call. = FALSE)
  basis <- match.arg(basis, c("decile", "adaptive", "edge"), several.ok = TRUE)
  lambda_scale <- match.arg(lambda_scale); inflate <- match.arg(inflate)
  stopifnot(is.numeric(tau), length(tau) == 1, tau > 0, tau < 1)
  X <- as.matrix(X); y <- as.numeric(y); n <- length(y)
  if (nrow(X) != n) stop("X and y have different numbers of observations.", call. = FALSE)
  if (!all(y %in% c(0, 1))) stop("y must contain only 0 and 1.", call. = FALSE)
  prev <- mean(y)
  bal  <- min(prev, 1 - prev)
  if (bal < 0.35 && ncol(X) / n >= 0.05) {
    warning(sprintf(
      "calm.gof(): outcome prevalence is %.3f at p/n = %.3f. The reference is validated for
outcomes that are not strongly unbalanced; in this regime the level is lost as the prevalence
falls (see the Scope section of ?calm.gof). Treat the p-value as indicative only.",
      prev, ncol(X) / n), call. = FALSE)
  } else if (bal < 0.20) {
    warning(sprintf(
      "calm.gof(): outcome prevalence is %.3f. At fixed dimension the decile basis is unaffected,
but the smooth basis degrades with strong imbalance (see the Scope section of ?calm.gof).",
      prev), call. = FALSE)
  }
  if (ncol(X) >= n) stop("calm.gof() needs p < n; the reference is built from the sample covariance of X.", call. = FALSE)
  lambda_glmnet <- if (lambda_scale == "glmnet") lambda else lambda / n
  if (lambda_scale == "glmnet") lambda <- n * lambda

  pc  <- .calm_pieces(X, y, lambda, G)
  ref <- .calm_reference(X, y, lambda, pc)
  Sigma0 <- ref$Sigma
  kappa <- ncol(X) / n
  rho_hat <- sqrt(ref$a2 / (ref$a2 + ref$s2))
  # basis-specific inflation for the selection effect: the term lives almost entirely in the
  # cubic direction, so it is large for a basis that keeps that direction and small for one that does not
  gof <- function(which) if (inflate == "none") 0 else
    kappa * c(dec = 0.24, e2 = 0.11, e3 = 1.15)[[which]]

  out <- list(prevalence = prev, lambda = lambda, lambda_glmnet = lambda_glmnet, G = G, kappa = kappa,
              rho_hat = rho_hat, tau = tau, inflate = inflate,
              est_signal_sd = ref$est_signal_sd, scale_c = ref$c)
  if ("decile" %in% basis)
    out$SC.HL <- list(statistic = pc$S_dec,
                      p.value = .calm_pvalue(pc$S_dec, (1 + gof("dec")) * Sigma0),
                      g_hat = gof("dec"))
  if ("edge" %in% basis)
    out$SC.EDGE <- list(statistic = pc$S_edge,
                        p.value = .calm_pvalue(pc$S_edge, (1 + gof("e3")) * Sigma0, pc$Z),
                        degree = 3L, g_hat = gof("e3"))
  if ("adaptive" %in% basis) {
    # keep the cubic when there is no selection effect to avoid (small aspect ratio), or when
    # the index is accurate enough to carry a cubic signal; never fall below degree two, since
    # rho_hat is a null quantity and falls when a departure is present
    k <- if (kappa < 0.05 || rho_hat^6 >= tau) 3L else 2L
    Zk <- .calm_edge_basis(pc, k)
    Sk <- if (is.null(Zk)) NA_real_ else .calm_qf(pc$v, Zk)
    gk <- gof(if (k == 3L) "e3" else "e2")
    out$SC.EDGE.adaptive <- list(statistic = Sk,
                                 p.value = .calm_pvalue(Sk, (1 + gk) * Sigma0, Zk),
                                 degree = k, rho_hat = rho_hat, g_hat = gk)
  }
  class(out) <- "calm.gof"
  out
}

#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept CALM
#' @concept penalized regression
#' @concept ridge regression
#' @concept high-dimensional
#' @concept Hosmer-Lemeshow
#' @export
print.calm.gof <- function(x, ...) {
  cat("\nCALM: closed-form goodness of fit for penalized logistic regression\n")
  cat(sprintf("lambda = %.4g (theory scale) = %.4g on glmnet's scale | G = %d | p/n = %.3f\n",
              x$lambda, x$lambda_glmnet, x$G, x$kappa))
  for (nm in c("SC.HL", "SC.EDGE", "SC.EDGE.adaptive")) {
    if (is.null(x[[nm]])) next
    cat(sprintf("  %-17s statistic %8.3f   p = %.4f", nm, x[[nm]]$statistic, x[[nm]]$p.value))
    if (nm == "SC.EDGE.adaptive")
      cat(sprintf("   [degree %d, rho_hat %.2f]", x[[nm]]$degree, x[[nm]]$rho_hat))
    cat("\n")
  }
  cat("\nSC.EDGE.adaptive is the statistic to prefer; see ?calm.gof for the scope of the reference.\n\n")
  invisible(x)
}
