# The projection goodness-of-fit test of Escanciano (2006), in the form studied for
# logistic regression by Liu et al. (2024, Sec. 2.1), with their model-based bootstrap.
#
#   T = n^-2 sum_i sum_j sum_l e_i e_j A_ijl,   e_i = y_i - p-hat_i,
#   A_ijl = integral over the unit sphere of I(X_i'w <= X_l'w) I(X_j'w <= X_l'w) dw,
#
# with dw the uniform probability measure. Put u = X_i - X_l and v = X_j - X_l. Then
#   u, v both nonzero:  A_ijl = (pi - angle(u, v)) / (2 pi)
#   exactly one zero:   A_ijl = 1/2
#   both zero:          A_ijl = 1
# which is Escanciano's (2006, Appendix) closed form up to the constant 2 pi. The printed
# formula in Liu et al. reads arccos(...) without "pi -"; the integral fixes the orientation.
# There is no length factor: A_ijl depends on directions only. X is the model matrix without
# its intercept, on the raw scale (the weight is not invariant to rescaling a column).

# The n x n weight A_ij = sum_l A_ijl. Pure R: one n x n cross-product and one n x n acos
# per l, O(n^3) time and O(n^2) memory. A difference shorter than eps counts as zero.
.proj_weight <- function(X, tol = 1e-10) {
  X <- as.matrix(X)
  n <- nrow(X)
  eps <- tol * max(1, sqrt(max(rowSums(X * X))))
  if (ncol(X) == 1L) {
    # one covariate: w is +1 or -1 with probability 1/2 each, so A_ijl is
    # (I(u <= 0) I(v <= 0) + I(u >= 0) I(v >= 0)) / 2, with |u| <= eps read as u = 0
    d  <- outer(X[, 1], X[, 1], "-")          # d[i, l] = x_i - x_l
    Lo <- (d <= eps) + 0
    Hi <- (d >= -eps) + 0
    return(0.5 * (tcrossprod(Lo) + tcrossprod(Hi)))
  }
  S  <- matrix(0, n, n)                        # sum over l of the angles
  Zm <- matrix(0, n, n)                        # Zm[i, l] = 1 when X_i = X_l
  dg <- seq.int(1L, n * n, by = n + 1L)
  for (l in seq_len(n)) {
    D  <- X - rep(X[l, ], each = n)
    nm <- sqrt(rowSums(D * D))
    z  <- nm <= eps
    nm[z] <- 1
    C  <- tcrossprod(D / nm)                   # cosines between unit differences
    C[dg] <- 1                                 # i = j: angle 0 exactly
    if (any(z)) {                              # a zero difference: angle 0, i.e. A_ijl = 1/2
      C[z, ] <- 1
      C[, z] <- 1
      Zm[z, l] <- 1
    }
    S <- S + acos(pmin.int(pmax.int(C, -1), 1))
  }
  # (pi - angle) / (2 pi) summed over l, plus 1/2 more where both differences are zero
  n / 2 - S / (2 * pi) + 0.5 * tcrossprod(Zm)
}

.proj_stat <- function(e, A) sum(e * (A %*% e)) / length(e)^2

#' Projection Goodness-of-Fit Test for Binary Regression
#'
#' @description
#' \code{projection.gof()} computes the projection test of Escanciano (2006) as defined for
#' logistic regression by Liu et al. (2024, Sec. 2.1), and refers it to their
#' model-based bootstrap. The residual-marked empirical process is taken along every
#' direction of the covariate space and its Cramer-von Mises norm is integrated over the
#' unit sphere, so the test is consistent against any departure of the mean function,
#' including departures that are invisible along the fitted linear predictor (where the
#' Stute-Zhu test looks).
#'
#' @details
#' With \eqn{e_i = y_i - \hat p_i} and \eqn{X_i} the covariate vector of observation
#' \eqn{i} (the model matrix without its intercept), the statistic is
#' \deqn{T = n^{-2} \sum_i \sum_j \sum_l e_i e_j A_{ijl},}
#' \deqn{A_{ijl} = \int_{S^p} I(X_i^T w \le X_l^T w) I(X_j^T w \le X_l^T w)\, dw,}
#' with \eqn{dw} the uniform probability measure on the unit sphere. For
#' \eqn{u = X_i - X_l} and \eqn{v = X_j - X_l} both nonzero,
#' \eqn{A_{ijl} = \{\pi - \angle(u, v)\} / (2\pi)}; it is \eqn{1/2} when exactly one of
#' them is zero and \eqn{1} when both are (Escanciano 2006, Appendix). The constant does
#' not affect the bootstrap p-value. Liu et al. print the weight as a multiple of
#' \eqn{\arccos(\cdot)}; the integral fixes the orientation as \eqn{\pi - \arccos(\cdot)},
#' and a Monte Carlo over random directions agrees with the form used here (see the
#' package tests).
#'
#' The weight uses the covariates on their own scale, as Liu et al. do. The angle is not
#' invariant to rescaling one column, so \code{scale = TRUE}, which standardizes each
#' column first, gives a different test; it is offered for covariates in unrelated units.
#'
#' \strong{Reference distribution.} The p-value comes from the model-based bootstrap of
#' Liu et al. (2024, Sec. 2.1; Dikta et al. 2006): \code{B} responses are drawn from the
#' fitted probabilities, the model is refitted to each with the same design, link and
#' offset, and \eqn{p = (1 + \#\{T^* \ge T\}) / (B + 1)}. A refit that fails scores
#' \eqn{T^* = -\infty}, which can only make the test conservative; the count is returned.
#' The default \code{B = 1000} is the number of bootstrap samples Liu et al. use in both
#' of their data examples.
#'
#' \strong{Cost.} The weight matrix is computed once, in \eqn{O(n^3 p)} time and
#' \eqn{O(n^2)} memory, and each bootstrap replicate then costs one refit and one
#' quadratic form. In pure R the weight takes about 0.2 s at \eqn{n = 200} and 3 to 5 s
#' at \eqn{n = 500}, and the whole test with \code{B = 1000} about 1 s and 6 s; at \eqn{n = 2000} it takes minutes and several n-by-n matrices of
#' memory. With a single covariate the weight has an exact rank form and is much faster.
#'
#' @param object A fitted binary \code{\link[stats]{glm}} (\code{family = binomial},
#'   any link) with unit prior weights and at least one covariate.
#' @param B Number of bootstrap replicates (default 1000, as in Liu et al.'s examples).
#' @param scale Logical; standardize each covariate column before computing the angles.
#'   Default \code{FALSE}, the statistic of Liu et al.
#' @param tol A difference vector shorter than \code{tol} times the largest covariate
#'   norm (at least 1) is treated as the zero vector.
#'
#' @return An object of class \code{"htest"} with elements \code{statistic} (\eqn{T}),
#'   \code{parameter} (\code{B}), \code{p.value}, \code{method}, \code{data.name}, and
#'   additionally \code{boot} (the \code{B} bootstrap statistics), \code{n_failed}
#'   (refits that failed) and \code{n_nonconverged} (refits that did not converge; their
#'   fitted values are still used).
#'
#' @references
#' Escanciano JC (2006). "A consistent diagnostic test for regression models using
#' projections." \emph{Econometric Theory}, 22(6), 1030--1051.
#' \doi{10.1017/S0266466606060506}
#'
#' Liu H, Li X, Chen F, Haerdle W, Liang H (2024). "A comprehensive comparison of
#' goodness-of-fit tests for logistic regression models." \emph{Statistics and
#' Computing}, 34, 175. \doi{10.1007/s11222-024-10487-5}
#'
#' Dikta G, Kvesic M, Schmidt C (2006). "Bootstrap approximations in model checks for
#' binary data." \emph{Journal of the American Statistical Association}, 101(474),
#' 521--530. \doi{10.1198/016214505000001032}
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' n  <- 150
#' x1 <- rnorm(n); x2 <- rnorm(n)
#' y  <- rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2))
#' fit <- glm(y ~ x1 + x2, family = binomial())
#' projection.gof(fit, B = 99)
#'
#' \donttest{
#' ## an omitted interaction: invisible to a test that looks only along the
#' ## fitted linear predictor, visible along other directions
#' y2  <- rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2 + 1.5 * x1 * x2))
#' bad <- glm(y2 ~ x1 + x2, family = binomial())
#' projection.gof(bad)          # B = 1000
#' }
#'
#' @seealso \code{\link{run.all.gof}}, where the test is the row \code{"Projection"}.
#' @concept goodness-of-fit
#' @concept logistic regression
#' @concept bootstrap
#' @export
projection.gof <- function(object, B = 1000, scale = FALSE, tol = 1e-10) {
  if (!inherits(object, "glm") || object$family$family != "binomial")
    stop("projection.gof: 'object' must be a fitted binomial glm.")
  B <- as.integer(B)
  if (length(B) != 1L || is.na(B) || B < 1L)
    stop("projection.gof: 'B' must be a positive integer.")
  y <- as.numeric(object$y)
  n <- length(y)
  if (!all(y %in% c(0, 1)) || any(object$prior.weights != 1))
    stop("projection.gof: needs binary (0/1) data with unit prior weights.")
  X <- stats::model.matrix(object)
  Z <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(Z) == 0L)
    stop("projection.gof: the model has no covariates.")
  if (isTRUE(scale)) {
    Z <- Z[, apply(Z, 2, stats::sd) > 0, drop = FALSE]
    Z <- base::scale(Z)
  }
  A    <- .proj_weight(Z, tol = tol)
  mu   <- as.numeric(stats::fitted(object))
  Tobs <- .proj_stat(y - mu, A)

  off <- object$offset
  if (is.null(off)) off <- rep(0, n)
  fam <- object$family
  Tb  <- numeric(B)
  n_failed <- 0L
  n_nc     <- 0L
  for (b in seq_len(B)) {
    yb <- stats::rbinom(n, 1, mu)
    fb <- tryCatch(suppressWarnings(stats::glm.fit(X, yb, family = fam, offset = off)),
                   error = function(e) NULL)
    if (is.null(fb) || any(!is.finite(fb$fitted.values))) {
      Tb[b] <- -Inf
      n_failed <- n_failed + 1L
      next
    }
    if (!isTRUE(fb$converged)) n_nc <- n_nc + 1L
    Tb[b] <- .proj_stat(yb - fb$fitted.values, A)
  }
  structure(list(
    statistic = c(T = Tobs),
    parameter = c(B = B),
    p.value   = (1 + sum(Tb >= Tobs)) / (B + 1),
    method    = paste0("Projection goodness-of-fit test (Escanciano 2006; Liu et al. 2024), ",
                       "model-based bootstrap", if (isTRUE(scale)) ", standardized covariates"),
    data.name = paste(deparse(stats::formula(object)), collapse = " "),
    boot = Tb, n_failed = n_failed, n_nonconverged = n_nc),
    class = "htest")
}

# Battery wrapper. The weight is O(n^3) in time and O(n^2) in memory, so the row is
# skipped above max_n (default 3000) unless the caller raises it.
gof_proj <- function(ctx, opts = list()) {
  if (!ctx$has_model)
    return(list(Statistic = NA, df = NA, p_value = NA, Note = "needs a glm model"))
  B     <- if (is.null(opts$B)) 1000L else as.integer(opts$B)
  max_n <- if (is.null(opts$max_n)) 3000 else opts$max_n
  if (ctx$n > max_n)
    return(list(Statistic = NA, df = NA, p_value = NA,
                Note = sprintf(paste0("Not run: n = %d > %d (O(n^3) weight); raise it with ",
                                      "control = list(Projection = list(max_n = ...))"),
                               ctx$n, max_n)))
  r <- projection.gof(ctx$model, B = B, scale = isTRUE(opts$scale))
  list(Statistic = unname(r$statistic), df = NA_real_, p_value = r$p.value,
       Note = paste0(B, " model-based bootstrap refits",
                     if (isTRUE(opts$scale)) "; standardized covariates" else "",
                     if (r$n_failed > 0) sprintf("; %d failed refits", r$n_failed) else ""))
}
