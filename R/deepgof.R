# =====================================================================
# deepgof1() -- DeepGOF-1: a pretrained goodness-of-fit test for logistic regression
#
# The test statistic is a small convolutional network (18,273 frozen parameters) that
# reads the fitted model's RESIDUAL MAP -- a 6x6 grid of standardized residual sums over
# the ranks of the two strongest covariates -- and the p-value is the rank of the observed
# score inside the analyst's own parametric bootstrap. The level is therefore a property
# of the calibration, not of what the network learned.
#
# The network was trained ONCE, offline, on simulated departures and ships frozen in
# R/sysdata.rda. Nothing is trained here: the forward pass is a handful of small matrix
# multiplies in base R, so the package needs no Python, no GPU and no extra dependency.
#
# Reference implementation for
#   Ebrahim, E.K. (2026). DeepGOF-1: a pretrained convolutional goodness-of-fit test for
#   logistic regression with a computable consistency certificate.
# The training corpus, training code, benchmark harness and every per-replicate p-value
# behind that paper are archived separately (see the paper's data-availability statement);
# only the deployable test lives in this package.
# =====================================================================

## ---- the forward pass -------------------------------------------------------------------
## im2col: (C,H,W) -> ((C*9) x (H*W)) with zero padding of 1, so a convolution is ONE
## matrix multiply and 200 forward passes stay fast enough in base R.
.dg_im2col3 <- function(x) {
  C <- dim(x)[1]; H <- dim(x)[2]; Wd <- dim(x)[3]
  xp <- array(0, c(C, H + 2L, Wd + 2L))
  xp[, 2:(H + 1L), 2:(Wd + 1L)] <- x
  out <- matrix(0, C * 9L, H * Wd)
  k <- 0L
  for (di in 0:2) for (dj in 0:2) {
    k <- k + 1L
    patch <- xp[, (1L + di):(H + di), (1L + dj):(Wd + dj), drop = FALSE]
    out[((k - 1L) * C + 1L):(k * C), ] <- matrix(patch, C, H * Wd)
  }
  out
}

.dg_conv3 <- function(x, Wt, b) {
  C <- dim(x)[1]; H <- dim(x)[2]; Wd <- dim(x)[3]; O <- dim(Wt)[1]
  Wm <- matrix(0, O, C * 9L)
  k <- 0L
  for (di in 1:3) for (dj in 1:3) {
    k <- k + 1L
    Wm[, ((k - 1L) * C + 1L):(k * C)] <- Wt[, , di, dj, drop = FALSE][, , 1, 1]
  }
  array(Wm %*% .dg_im2col3(x) + b, c(O, H, Wd))
}

.dg_maxpool2 <- function(x) {
  C <- dim(x)[1]; H <- dim(x)[2]; Wd <- dim(x)[3]
  o <- array(0, c(C, H %/% 2L, Wd %/% 2L))
  for (i in seq_len(H %/% 2L)) for (j in seq_len(Wd %/% 2L))
    o[, i, j] <- apply(x[, (2*i-1):(2*i), (2*j-1):(2*j), drop = FALSE], 1, max)
  o
}

## Score one STANDARDIZED map (length K*K).
## THE PORT TRAP: PyTorch's view(-1,1,K,K) fills ROW-major and R's array() fills
## COLUMN-major, so a naive array(v, c(1,K,K)) silently TRANSPOSES the map. byrow = TRUE
## restores PyTorch's ordering, which is also the order the map builder below emits.
## An earlier port that got this wrong agreed to 7.5e-3 -- close enough to look right.
.dg_score <- function(z, M = deepgof1_weights) {
  K <- M$K
  x <- array(matrix(as.numeric(z), K, K, byrow = TRUE), c(1L, K, K))
  x <- pmax(.dg_conv3(x, M$W[["conv.0.weight"]], M$W[["conv.0.bias"]]), 0)
  x <- pmax(.dg_conv3(x, M$W[["conv.2.weight"]], M$W[["conv.2.bias"]]), 0)
  x <- .dg_maxpool2(x)
  x <- pmax(.dg_conv3(x, M$W[["conv.5.weight"]], M$W[["conv.5.bias"]]), 0)
  e <- c(apply(x, 1, max), apply(x, 1, mean))          # GlobalMax then GlobalMean -> 64
  h <- pmax(as.vector(M$W[["head.0.weight"]] %*% e + M$W[["head.0.bias"]]), 0)
  as.numeric(M$W[["head.2.weight"]] %*% h + M$W[["head.2.bias"]])
}

## ---- the residual map -------------------------------------------------------------------
## Identical to the training-time builder; changing it invalidates the shipped weights.
.dg_map <- function(fit, K = 6L) {
  d  <- stats::model.frame(fit)
  mm <- stats::model.matrix(fit)
  vars <- colnames(mm)[-1]
  if (length(vars) < 2L) stop("deepgof1() needs at least two covariates", call. = FALSE)
  b  <- stats::coef(fit)[vars]
  sc <- abs(b) * apply(mm[, vars, drop = FALSE], 2, stats::sd)
  ax <- vars[sort(order(sc, decreasing = TRUE)[1:2])]
  ph <- as.numeric(stats::fitted(fit))
  r  <- stats::model.response(d) - ph
  n  <- length(r)
  bn <- function(v) pmin(K, 1L + floor(K * (rank(v, ties.method = "first") - 1) / n))
  idx <- factor((bn(mm[, ax[1]]) - 1L) * K + bn(mm[, ax[2]]), levels = 1:(K * K))
  s  <- tapply(r, idx, sum); vv <- tapply(ph * (1 - ph), idx, sum)
  s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
  list(map = as.numeric(s) / sqrt(pmax(as.numeric(vv), 1e-8)), axes = ax)
}

#' DeepGOF-1: a pretrained goodness-of-fit test for logistic regression
#'
#' Tests whether a fitted binomial \code{glm} is correctly specified, using a
#' convolutional network that was trained once, offline, on simulated departures and is
#' shipped frozen with this package. The analyst never trains anything: the network reads
#' the model's residual map and the p-value is the rank of the observed score within the
#' analyst's own parametric bootstrap, so the level does not depend on what the network
#' learned.
#'
#' The residual map is a \code{K} by \code{K} grid over the empirical ranks of the two
#' covariates with the largest \eqn{|\hat\beta_j| \hat\sigma_j}; each cell holds a
#' standardized residual sum, approximately standard normal under a correct model. Misfit
#' therefore has a location on the map -- an omitted quadratic paints a stripe, an omitted
#' interaction a saddle -- which is what the convolutional statistic reads.
#'
#' The test is a small-sample instrument. Against the classical partition tests it gains
#' most at \eqn{n} of 50 to 200 and the gain decays as \eqn{n} grows; because the grid uses
#' only two covariates, misfit that lives off those axes is harder for it to see than for
#' covariate-space or smoothing tests. See \code{\link{run.all.gof}} to run it alongside
#' the classical battery.
#'
#' @param fit a fitted \code{glm} with \code{family = binomial()} and at least two
#'   covariates.
#' @param B number of parametric-bootstrap replicates. The p-value lies on a grid of
#'   \code{1/(B+1)}, so \code{B = 199} makes the nominal .05 attainable exactly.
#' @param K grid resolution. Leave at 6: the shipped weights were trained at \code{K = 6}
#'   and are not valid at any other resolution.
#'
#' @return An object of class \code{"deepgof1"}: a list with \code{statistic} (the observed
#'   score), \code{p.value}, \code{B}, \code{K}, \code{axes} (the two covariates the grid
#'   was built on), \code{boot} (the \code{B} bootstrap scores) and \code{method}.
#'
#' @section Reproducibility:
#' A bootstrap refit that fails to converge is scored \code{+Inf}, so it counts against
#' rejection -- the conservative direction. Set a seed before calling for a reproducible
#' p-value.
#'
#' @references
#' Ebrahim, E.K. (2026). DeepGOF-1: a pretrained convolutional goodness-of-fit test for
#' logistic regression with a computable consistency certificate.
#'
#' Besag, J. and Clifford, P. (1989). Generalized Monte Carlo significance tests.
#' \emph{Biometrika} \strong{76}, 633--642. \doi{10.1093/biomet/76.4.633}
#'
#' @examples
#' set.seed(1)
#' n  <- 150
#' x1 <- runif(n, -3, 3); x2 <- rnorm(n)
#' # a model with an omitted quadratic term
#' y  <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.9 * (x1^2 - mean(x1^2))))
#' fit <- glm(y ~ x1 + x2, family = binomial())
#' deepgof1(fit, B = 49)   # B = 49 to keep the example fast; use the default in practice
#'
#' @seealso \code{\link{run.all.gof}}, \code{\link{ef.gof}}, \code{\link{legoft}}
#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept DeepGOF-1
#' @concept convolutional neural network
#' @concept pretrained
#' @export
deepgof1 <- function(fit, B = 199L, K = 6L) {
  if (!inherits(fit, "glm") || fit$family$family != "binomial")
    stop("deepgof1() expects a glm fitted with family = binomial()", call. = FALSE)
  if (K != deepgof1_weights$K)
    stop("the shipped DeepGOF-1 weights are valid only at K = ", deepgof1_weights$K,
         call. = FALSE)
  M <- deepgof1_weights
  score <- function(m) .dg_score((m - M$mu) / M$sd, M)

  obs  <- .dg_map(fit, K)
  Sobs <- score(obs$map)

  ph  <- as.numeric(stats::fitted(fit))
  dat <- stats::model.frame(fit)
  yn  <- names(dat)[1]
  frm <- stats::formula(fit)
  Sb  <- numeric(B)
  for (b in seq_len(B)) {
    dat[[yn]] <- stats::rbinom(length(ph), 1L, ph)
    fb <- tryCatch(suppressWarnings(stats::glm(frm, data = dat, family = stats::binomial())),
                   error = function(e) NULL)
    ## a failed refit counts AGAINST rejection (conservative); -Inf would deflate p
    Sb[b] <- if (is.null(fb)) Inf else score(.dg_map(fb, K)$map)
  }
  structure(list(statistic = Sobs,
                 p.value   = (1 + sum(Sb >= Sobs)) / (B + 1),
                 B = B, K = K, axes = obs$axes, boot = Sb,
                 method    = "DeepGOF-1: pretrained residual-map goodness-of-fit test",
                 data.name = deparse(substitute(fit))),
            class = "deepgof1")
}

#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept DeepGOF-1
#' @concept convolutional neural network
#' @concept pretrained
#' @export
print.deepgof1 <- function(x, ...) {
  cat("\n\t", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  cat(sprintf("S = %.4f, B = %d, p-value = %.4f\n", x$statistic, x$B, x$p.value))
  cat(sprintf("grid: %d x %d over ranks of %s and %s\n", x$K, x$K, x$axes[1], x$axes[2]))
  cat("alternative hypothesis: the logistic model is misspecified\n\n")
  invisible(x)
}
