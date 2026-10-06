## deepgof_certificate.R -- can DeepGOF-1 see a departure BEFORE you spend a bootstrap on it?
##
## The scored statistic S is built from affine maps, rectifiers, max-pool and componentwise
## maxima, so it is piecewise affine and its recession function
##      Sinf(u) = lim_{t->Inf} S(z0 + t u) / t
## exists, is independent of z0, and IS THE SAME NETWORK WITH EVERY OFFSET SET TO ZERO. Under a
## fixed alternative the residual map grows like sqrt(n) times a deterministic direction while the
## bootstrap scores stay bounded, so the power of the deployed test tends to one exactly when
## Sinf(direction) > 0. One forward pass on the frozen weights therefore decides, in advance,
## whether this test can detect a named departure.
##
## THE OFFSETS. Zeroing them is not an approximation: on the shipped weights the closed form
## matches slopes measured far along each ray to 1.3e-15. Note that the direction is scaled by the
## training sd ONLY -- the centring mu cancels in a difference of maps, and subtracting it here
## would silently add a bias term to a quantity that is supposed to be homogeneous.
##
## THE SCALE. Sinf is positively homogeneous, so doubling a direction doubles the slope and the
## SIGN -- the part that certifies consistency -- never moves. The number is still reported at
## unit Euclidean norm, as the paper's certificate table is, so that margins from different
## departures can be read against each other and against the worst direction on the sphere
## (-0.281 for these weights). Without that, a "margin" would say more about how loudly the
## departure was written down than about the test.

## the frozen network with its biases removed: the positively homogeneous part of S
.dg_weights_nobias <- function(M) {
  W <- M$W
  for (nm in names(W)) if (grepl("bias", nm)) W[[nm]] <- W[[nm]] * 0
  M$W <- W
  M
}

## The cell-wise direction of a departure, on the SAME grid the test would use: the map of the
## default reading of deepgof1() with the observed residual y - pihat replaced by its expectation
## under the alternative, pistar - pihat. Same covariates, same axis rule, same rank bins, same
## standardization -- if this drifted from deepgof1() the certificate would describe a test nobody
## runs. Ties among covariate values go by row order here; deepgof1() breaks them at random, which
## moves single records between neighbouring cells and not the limit direction.
.dg_direction <- function(fit, pistar, K = 6L) {
  mm <- stats::model.matrix(fit)
  X  <- .dg_covariates(fit)
  if (is.null(X)) {
    v <- colnames(mm)[colnames(mm) != "(Intercept)"]
    X <- mm[, v[!is.na(stats::coef(fit)[v])], drop = FALSE]
    tcols <- stats::setNames(lapply(colnames(X), function(v) match(v, colnames(mm))), colnames(X))
  } else {
    tcols <- .dg_term_cols(fit, mm, colnames(X))
  }
  if (ncol(X) < 1L) stop("deepgof1() needs at least one covariate", call. = FALSE)
  ph <- as.numeric(fit$fitted.values)
  n  <- length(ph)
  if (ncol(X) == 1L) {
    ## one covariate: the 36 quantile cells along it, as deepgof1() reads it
    idx <- .dg_cells(X, K)$idx[[1]]
    ax  <- colnames(X)
  } else {
    ax  <- .dg_pick_axes(mm, tcols, stats::coef(fit))
    bn  <- function(v) pmin(K, 1L + floor(K * (rank(v, ties.method = "first") - 1) / n))
    idx <- as.integer((bn(X[, ax[1]]) - 1L) * K + bn(X[, ax[2]]))
  }
  list(map = as.numeric(.dg_cellmap(idx, pistar - ph, ph, K)), axes = ax)
}

#' Can DeepGOF-1 detect this departure? A certificate from the frozen weights
#'
#' Decides, without running the test, whether \code{\link{deepgof1}} is consistent against a
#' departure you name: its power tends to one precisely when the certificate below is positive.
#' The answer costs one forward pass on the shipped weights -- no bootstrap, no simulation, no
#' data beyond the fitted model.
#'
#' Give the departure in whichever form you have it. \code{pi_star} is the vector of true
#' success probabilities under the alternative, one per observation of \code{fit}.
#' \code{eta_shift} is the term the fitted linear predictor is missing, so that the truth is
#' \code{plogis(eta_hat + eta_shift)} -- an omitted quadratic is \code{eta_shift = c * x^2}.
#' \code{direction} skips the data entirely and takes a \code{K} by \code{K} map direction, which
#' is what a theoretical question ("a checkerboard at the cell scale") usually looks like.
#'
#' A non-positive certificate does not mean the test never rejects; it means the guarantee of
#' power tending to one is unavailable for that direction, which for this network happens on a
#' thin cone of sub-resolution oscillation -- the failure mode every binned statistic shares.
#' Because the grid reads two covariates, a departure that integrates to zero within every cell
#' has direction zero and is certified blind, which is the representation limit of the method
#' stated as a computation rather than a caveat.
#'
#' With a fitted model the direction is built on the map of the default reading of
#' \code{\link{deepgof1}} (\code{reading = "axes"}): the two covariates of its axis rule, or the
#' 36 cells along the only covariate. For the all-pairs reading, name the pair's map yourself
#' through \code{direction}.
#'
#' @param fit a fitted \code{glm} with \code{family = binomial()}: the null model whose
#'   adequacy is in question. Not needed if \code{direction} is given.
#' @param pi_star true success probabilities under the alternative, one per observation.
#' @param eta_shift the omitted term on the linear-predictor scale, one per observation; used
#'   only when \code{pi_star} is not given.
#' @param direction a \code{K} by \code{K} matrix, or a vector of length \code{K * K} in row-major
#'   order, giving the map direction directly.
#' @param K grid resolution. Leave at 6: the shipped weights are valid at no other resolution.
#'
#' @return An object of class \code{"gof_certificate"}: \code{certificate} (the recession slope at
#'   unit map norm; positive certifies consistency), \code{verdict}, \code{mirror} (the slope of
#'   the same departure with its sign flipped, since a test may see a departure one way and not
#'   the other), \code{direction} (normalized), \code{axes} and \code{method}.
#'
#' @references
#' Ebrahim EK, Hussein OAE-A, El-Kotory A (2026). "Where Does a Logistic Risk Model Fail?
#' An Audited Neural Goodness-of-Fit Test for Model Development and External Validation."
#' arXiv:2609.29575 [stat.ME]. \doi{10.48550/arXiv.2609.29575}
#'
#' @examples
#' set.seed(1)
#' n  <- 200
#' x1 <- runif(n, -3, 3); x2 <- rnorm(n)
#' y  <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2))
#' fit <- glm(y ~ x1 + x2, family = binomial())
#'
#' # an omitted quadratic in x1: is this test consistent against it?
#' gof_certificate(fit, eta_shift = 0.9 * (x1^2 - mean(x1^2)))
#'
#' # a map direction named directly: a threshold departure along the first axis
#' gof_certificate(direction = matrix(rep(c(-1, -1, -1, 1, 1, 1), each = 6), 6, 6, byrow = TRUE))
#'
#' # the blind cone is thin but not empty: multi-start minimization over the sphere reaches
#' # -0.281 for these weights, at a direction that oscillates cell by cell. Directions like
#' # that are the ones this function is for.
#'
#' @seealso \code{\link{deepgof1}}
#' @concept goodness-of-fit
#' @concept consistency
#' @concept DeepGOF-1
#' @export
gof_certificate <- function(fit, pi_star = NULL, eta_shift = NULL, direction = NULL, K = 6L) {
  M <- deepgof1_weights
  if (K != M$K)
    stop("the shipped DeepGOF-1 weights are valid only at K = ", M$K, call. = FALSE)
  M0 <- .dg_weights_nobias(M)
  slope <- function(v) .dg_score(v / M$sd, M0)      # sd only: see the note at the top of the file

  axes <- NULL
  if (!is.null(direction)) {
    v <- if (is.matrix(direction)) as.numeric(t(direction)) else as.numeric(direction)
    if (length(v) != K * K)
      stop("direction must have ", K * K, " entries (a ", K, " by ", K, " map)", call. = FALSE)
  } else {
    if (missing(fit) || !inherits(fit, "glm") || fit$family$family != "binomial")
      stop("give either direction, or a glm fitted with family = binomial()", call. = FALSE)
    if (is.null(pi_star)) {
      if (is.null(eta_shift))
        stop("name the departure: pi_star, eta_shift, or direction", call. = FALSE)
      eta_shift <- as.numeric(eta_shift)
      if (length(eta_shift) != length(fit$fitted.values))
        stop("eta_shift must have one entry per observation", call. = FALSE)
      pi_star <- stats::plogis(stats::qlogis(as.numeric(fit$fitted.values)) + eta_shift)
    }
    pi_star <- as.numeric(pi_star)
    if (length(pi_star) != length(fit$fitted.values))
      stop("pi_star must have one entry per observation", call. = FALSE)
    if (any(pi_star < 0 | pi_star > 1)) stop("pi_star must lie in [0, 1]", call. = FALSE)
    dd <- .dg_direction(fit, pi_star, K)
    v <- dd$map; axes <- dd$axes
  }

  nv <- sqrt(sum(v^2))
  if (nv > 0) v <- v / nv                          # report the slope per unit map norm
  cert <- slope(v)
  ## A departure that integrates to zero in every cell leaves NO direction at all -- the map does
  ## not grow, so the question "does the slope have the right sign" is the wrong question. Say so
  ## rather than returning a slope of zero that reads like a near miss.
  flat <- max(abs(v)) < 1e-10
  structure(list(certificate = cert,
                 mirror      = slope(-v),
                 verdict     = if (flat) "no direction: the departure cancels within every cell"
                               else if (cert > 0) "detectable: power tends to one"
                               else "not certified: this direction lies in the blind cone",
                 direction   = matrix(v, K, K, byrow = TRUE),
                 flat = flat, K = K, axes = axes,
                 method = "DeepGOF-1 consistency certificate (recession slope of the frozen network)"),
            class = "gof_certificate")
}

#' @concept goodness-of-fit
#' @concept DeepGOF-1
#' @export
print.gof_certificate <- function(x, ...) {
  cat("\n\t", x$method, "\n\n", sep = "")
  cat(sprintf("certificate = %.4f (mirror %.4f)\n", x$certificate, x$mirror))
  if (!is.null(x$axes))
    cat(sprintf("grid: %d x %d over ranks of %s and %s\n", x$K, x$K, x$axes[1], x$axes[2]))
  cat(x$verdict, "\n\n", sep = "")
  invisible(x)
}
