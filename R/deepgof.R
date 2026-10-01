# =====================================================================
# deepgof1() -- DeepGOF-1: a pretrained goodness-of-fit test for logistic regression
#
# The test statistic is a small convolutional network (18,273 frozen parameters) that
# reads the fitted model's RESIDUAL MAP -- a 6x6 grid of standardized residual sums over
# the ranks of two covariates (the two strongest by default, or every pair in turn with
# reading = "allpairs") -- and the p-value is the rank of the observed score inside the
# analyst's own parametric bootstrap. The level is therefore a property of the
# calibration, not of what the network learned.
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
## `tb` breaks ties among covariate values: a random permutation of the rows, drawn once by
## deepgof1() and passed unchanged to the observed map and to every bootstrap map. With
## tb = NULL, ties go by row order, as in 2.7.0 and the training corpus (whose covariates
## were continuous, so it had no ties). Row order must not be used on real data: a file
## sorted by the outcome then puts tied rows into cells by their outcome, the bootstrap
## outcomes are not sorted, and a correct model can be rejected.
.dg_map <- function(fit, K = 6L, tb = NULL)
  .dg_map_core(stats::model.matrix(fit), fit$y, as.numeric(stats::fitted(fit)),
               stats::coef(fit), K, tb)

## the same map from a design matrix, outcome, fitted risks and coefficients, so a bootstrap
## refit by glm.fit() needs no model frame
.dg_map_core <- function(mm, y, ph, coefs, K = 6L, tb = NULL) {
  vars <- colnames(mm)[-1]
  if (length(vars) < 2L) stop("deepgof1() needs at least two covariates", call. = FALSE)
  b  <- coefs[vars]
  sc <- abs(b) * apply(mm[, vars, drop = FALSE], 2, stats::sd)
  ax <- vars[sort(order(sc, decreasing = TRUE)[1:2])]
  r  <- y - ph
  n  <- length(r)
  rk <- function(v) {
    if (is.null(tb)) return(rank(v, ties.method = "first"))
    o <- order(v, tb); rr <- integer(n); rr[o] <- seq_len(n); rr
  }
  bn <- function(v) pmin(K, 1L + floor(K * (rk(v) - 1) / n))
  idx <- factor((bn(mm[, ax[1]]) - 1L) * K + bn(mm[, ax[2]]), levels = 1:(K * K))
  s  <- tapply(r, idx, sum); vv <- tapply(ph * (1 - ph), idx, sum)
  s[is.na(s)] <- 0; vv[is.na(vv)] <- 0
  list(map = as.numeric(s) / sqrt(pmax(as.numeric(vv), 1e-8)), axes = ax)
}

## ---- the all-pairs reading ------------------------------------------------------------------
## The maps are drawn over the model's COVARIATES, not its model-matrix columns: a covariate that
## enters as ns(x, 3), poly(x, 2) or I(x^2) is ranked by x itself, and an interaction adds no axis
## of its own. Ranking by a spline basis column or by x^2 would fold the axis. Covariates that the
## model uses untransformed are read from the model frame; transformed ones are fetched from the
## data the model was fitted to, on the rows the fit kept. Factors and logicals enter by their
## codes. Returns NULL when some covariate cannot be recovered; deepgof1() then falls back to the
## model-matrix columns.
.dg_covariates <- function(fit) {
  mf <- stats::model.frame(fit)
  tl <- attr(stats::terms(mf), "term.labels")
  vars <- unique(unlist(lapply(tl, function(l) all.vars(str2lang(l)))))
  if (!length(vars)) return(matrix(numeric(0), nrow(mf), 0L))
  src <- fit$data
  get1 <- function(v) {
    if (v %in% names(mf) && is.null(dim(mf[[v]]))) return(mf[[v]])
    if (is.data.frame(src) && v %in% names(src)) {
      if (nrow(src) == nrow(mf) && is.null(fit$na.action)) return(src[[v]])
      ids <- rownames(mf)
      if (!anyDuplicated(rownames(src)) && all(ids %in% rownames(src)))
        return(src[ids, v, drop = TRUE])
    }
    ## a model fitted without 'data': look where the formula was written
    if (is.null(fit$na.action)) {
      e <- tryCatch(get(v, envir = environment(stats::formula(fit))), error = function(e) NULL)
      if (is.atomic(e) && is.null(dim(e)) && length(e) == nrow(mf)) return(e)
    }
    NULL
  }
  cols <- lapply(vars, get1)
  if (any(vapply(cols, is.null, TRUE))) return(NULL)
  X <- vapply(cols, function(v) {
    if (is.character(v)) v <- factor(v)
    as.numeric(if (is.factor(v)) as.integer(v) else v)
  }, numeric(nrow(mf)))
  X <- matrix(X, nrow(mf), dimnames = list(NULL, vars))
  ## a covariate constant on the fitted rows has no ranks to read
  X[, apply(X, 2, function(v) length(unique(v)) > 1L), drop = FALSE]
}

## The cells depend only on the covariate ranks, which the bootstrap does not change, so they are
## found once per call: one integer cell index per row for every pair of columns (or, with a
## single column, 36 quantile cells along it). Each replicate then needs only the cell sums of its
## own residuals. The cell order is the one .dg_map emits, so the network sees the same layout.
.dg_cells <- function(X, K = 6L, tb = NULL) {
  n <- nrow(X)
  rk <- function(v) {
    if (is.null(tb)) return(rank(v, ties.method = "first"))
    o <- order(v, tb); rr <- integer(n); rr[o] <- seq_len(n); rr
  }
  if (ncol(X) == 1L) {
    ## one covariate: the map is 36 quantile cells along its ranks, read row by row
    i <- pmin(K * K, 1L + floor(K * K * (rk(X[, 1]) - 1) / n))
    return(list(idx = list(as.integer(i)), pairs = matrix(colnames(X), 1L, 2L)))
  }
  bins <- apply(X, 2, function(v) pmin(K, 1L + floor(K * (rk(v) - 1) / n)))
  prs <- utils::combn(colnames(X), 2L)
  list(idx = lapply(seq_len(ncol(prs)), function(j)
         as.integer((bins[, prs[1L, j]] - 1L) * K + bins[, prs[2L, j]])),
       pairs = t(prs))
}

## ---- one engine for every reading ---------------------------------------------------------------
## For each covariate, the model-matrix columns of its main-effect terms: the columns of every term whose
## only variable is that covariate (x, ns(x, 3), poly(x, 2), I(x^2), a factor's dummies). Interaction
## terms belong to no single covariate and are left out.
.dg_term_cols <- function(fit, mm, vars) {
  tl <- attr(stats::terms(fit), "term.labels")
  tv <- lapply(tl, function(l) all.vars(str2lang(l)))
  asg <- attr(mm, "assign")
  stats::setNames(lapply(vars, function(v)
    which(asg %in% which(vapply(tv, function(z) length(z) == 1L && z == v, TRUE)))), vars)
}

## The axis rule over covariates: each covariate is scored by the standard deviation of its terms'
## total contribution to the linear predictor, and the two largest give the axes, in model order. For a
## covariate that enters as one untransformed column this is |b| sd(x), the rule of versions 2.7.0 and
## 2.8.0, so for such models the axes, the map and the p-value are unchanged.
.dg_pick_axes <- function(mm, tcols, coefs) {
  sc <- vapply(tcols, function(cl) {
    if (!length(cl)) return(0)
    b <- coefs[cl]; b[is.na(b)] <- 0
    stats::sd(as.numeric(mm[, cl, drop = FALSE] %*% b))
  }, 0)
  names(tcols)[sort(order(sc, decreasing = TRUE)[1:2])]
}

## Every statistic of one call from ONE set of bootstrap refits:
##   axes     the network score of the map over the two axes of the axis rule (chosen again per refit)
##   allpairs the largest network score over the maps of every pair of covariates
##   ss, max  the sum of squared cells and the largest |cell| of the axis-rule map (comparison
##            statistics; not offered by deepgof1())
## Returns the observed values, the B x 4 matrix of bootstrap values, and the observed maps.
.dg_engine <- function(fit, B, K, X, mm, tb, score, need = c("axes", "allpairs"), tcols = NULL) {
  ph <- as.numeric(stats::fitted(fit)); y <- as.numeric(fit$y); n <- length(y)
  rk <- function(v) {
    if (is.null(tb)) return(rank(v, ties.method = "first"))
    o <- order(v, tb); rr <- integer(n); rr[o] <- seq_len(n); rr
  }
  one_d <- ncol(X) == 1L
  bins <- if (one_d) NULL else apply(X, 2, function(v) pmin(K, 1L + floor(K * (rk(v) - 1) / n)))
  cells <- if ("allpairs" %in% need || one_d) .dg_cells(X, K, tb) else NULL
  if (is.null(tcols))            # X holds model-matrix columns: each is its own term
    tcols <- stats::setNames(lapply(colnames(X), function(v) match(v, colnames(mm))), colnames(X))
  stats_of <- function(r, p, coefs) {
    out <- c(axes = NA_real_, allpairs = NA_real_, ss = NA_real_, max = NA_real_)
    if (one_d) {
      m <- .dg_cellmap(cells$idx[[1]], r, p, K); s <- score(m)
      return(list(v = c(axes = s, allpairs = s, ss = sum(m^2), max = max(abs(m))), map = m,
                  axes = colnames(X), top = 1L))
    }
    ax <- NULL; mA <- NULL
    if (any(c("axes", "ss", "max") %in% need)) {
      ax <- .dg_pick_axes(mm, tcols, coefs)
      mA <- .dg_cellmap(as.integer((bins[, ax[1]] - 1L) * K + bins[, ax[2]]), r, p, K)
      out[c("axes", "ss", "max")] <- c(score(mA), sum(mA^2), max(abs(mA)))
    }
    top <- NA_integer_
    if ("allpairs" %in% need) {
      so <- vapply(cells$idx, function(i) score(.dg_cellmap(i, r, p, K)), 0)
      out["allpairs"] <- max(so); top <- which.max(so)
      attr(out, "pairscores") <- so
    }
    list(v = out, map = mA, axes = ax, top = top)
  }
  obs <- stats_of(y - ph, ph, stats::coef(fit))
  itc <- attr(stats::terms(fit), "intercept") > 0L
  boot <- matrix(NA_real_, B, 4L, dimnames = list(NULL, c("axes", "allpairs", "ss", "max")))
  for (b in seq_len(B)) {
    ys <- stats::rbinom(n, 1L, ph)
    fb <- tryCatch(suppressWarnings(stats::glm.fit(mm, ys, weights = fit$prior.weights,
                                                   offset = fit$offset, family = fit$family,
                                                   control = fit$control, intercept = itc)),
                   error = function(e) NULL)
    ## a failed refit counts AGAINST rejection (conservative); -Inf would deflate p
    boot[b, ] <- if (is.null(fb)) Inf else {
      pb <- as.numeric(fb$fitted.values)
      stats_of(ys - pb, pb, fb$coefficients)$v[c("axes", "allpairs", "ss", "max")]
    }
  }
  list(obs = obs, boot = boot, cells = cells)
}

## rank p-value of each statistic, and the p-value of the smaller of the axes and all-pairs p-values,
## calibrated exactly: the observed data and the B replicates are exchangeable under the null, so each
## of the B + 1 is given its own rank p-value within the set, and the observed minimum is ranked
## among the B + 1 minima
.dg_pvalues <- function(obs, boot) {
  B <- nrow(boot)
  p1 <- function(k) (1 + sum(boot[, k] >= obs[[k]])) / (B + 1)
  out <- vapply(colnames(boot), p1, 0)
  all_a <- c(obs[["axes"]], boot[, "axes"]); all_p <- c(obs[["allpairs"]], boot[, "allpairs"])
  within <- function(v) vapply(seq_along(v), function(i) sum(v >= v[i]) / length(v), 0)
  mins <- pmin(within(all_a), within(all_p))
  c(out, combined = sum(mins <= mins[1]) / (B + 1))
}

## the standardized map of one pair, from its precomputed cell index. The sums go through
## sum(), as in .dg_map's tapply(), so the two builders agree to the last bit; rowsum()
## accumulates differently and moves the score by about 1e-15.
.dg_cellmap <- function(i, r, ph, K = 6L) {
  f <- factor(i, levels = seq_len(K * K))
  s  <- vapply(split(r, f), sum, 0)
  vv <- vapply(split(ph * (1 - ph), f), sum, 0)
  s / sqrt(pmax(vv, 1e-8))
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
#' A residual map is a \code{K} by \code{K} grid over the empirical ranks of two
#' covariates; each cell holds a standardized residual sum, approximately standard normal
#' under a correct model. Misfit therefore has a location on the map -- an omitted
#' quadratic paints a stripe, an omitted interaction a saddle -- which is what the
#' convolutional statistic reads.
#'
#' Four readings of the map are offered. The default, \code{reading = "axes"}, draws one map,
#' over the two covariates whose terms contribute most to the linear predictor (the standard
#' deviation of each covariate's total contribution; for a covariate that enters as one
#' untransformed column this is \eqn{|\hat\beta_j| \hat\sigma_j}), and chooses them again in
#' every bootstrap replicate. A covariate that enters as \code{ns(x, 3)}, \code{poly(x, 2)},
#' \code{I(x^2)} or a factor is one covariate, and the map is drawn over the ranks of \code{x}
#' itself. When the misfit lies on the covariates with the strongest effects, among others that
#' carry little signal, this rule finds them almost every time and has the most power. It reads
#' a covariate by its effect in the fitted model, so it can pass over a covariate whose effect is
#' a pure U-shape with no linear slope.
#'
#' \code{reading = "allpairs"} scores the map of every pair of covariates and takes the
#' largest score as the statistic. The same maximum is taken in every bootstrap replicate, so
#' the p-value needs no correction for the choice of pair. It does not depend on the fitted
#' effects, so it finds a U-shaped covariate, and with few covariates (three or so) it has
#' more power than the default; with many covariates the maximum over \eqn{p(p-1)/2} pairs
#' costs power when one pair carries the misfit. The pair that reaches the maximum is returned
#' as \code{axes}, and its map as \code{map}, so the result also says where the misfit lies;
#' \code{covariates} restricts the pairs to a chosen set.
#'
#' \code{reading = "combined"} runs both from the same bootstrap refits and reports the
#' smaller of their two p-values, calibrated exactly: the observed data and the \code{B}
#' replicates are exchangeable under the null, so the observed minimum is ranked among the
#' \code{B + 1} minima. It costs no more than the all-pairs reading. \code{reading = "columns"}
#' is the rule of versions 2.7.0 and 2.8.0, over two model-matrix columns, kept to reproduce
#' earlier results; for models whose covariates all enter as one untransformed column it gives
#' the same p-value as the default.
#'
#' Transformed covariates are read from the data the model was fitted to; factors enter by
#' their level codes.
#'
#' With one covariate there is no pair to choose, and every reading is the same map of 36
#' quantile cells along its ranks. The shipped network was trained on two-covariate maps
#' only: on such models the bootstrap still gives it its level, but it has less power than
#' a network trained on this map would.
#'
#' Ties among covariate values, as with binary, categorical or rounded covariates, are
#' broken at random. The random order is drawn once per call and used for the observed map
#' and for every bootstrap map, so the p-value does not depend on the order of the rows in
#' the data. With heavily tied axes the p-value can vary noticeably from one seed to the
#' next; report the seed.
#'
#' The test is a small-sample instrument. Against the classical partition tests it gains
#' most at \eqn{n} of 50 to 200 and the gain decays as \eqn{n} grows; because each map uses
#' two covariates at a time, misfit that depends on three or more covariates jointly is
#' harder for it to see than for covariate-space or smoothing tests. It is not one of the tests
#' \code{\link{run.all.gof}} selects: call it directly on the same fitted model and read its
#' p-value beside the panel. See \code{\link{run.all.gof}} for
#' the classical battery.
#'
#' @param fit a fitted \code{glm} with \code{family = binomial()}, a 0/1 outcome and at
#'   least one covariate.
#' @param B number of parametric-bootstrap replicates. The p-value lies on a grid of
#'   \code{1/(B+1)}, so \code{B = 199} makes the nominal .05 attainable exactly.
#' @param K grid resolution. Leave at 6: the shipped weights were trained at \code{K = 6}
#'   and are not valid at any other resolution.
#' @param reading \code{"axes"} (the default), \code{"allpairs"}, \code{"combined"} or
#'   \code{"columns"} (the rule of versions 2.7.0 and 2.8.0). See Details for when to use which.
#' @param covariates optional character vector naming the covariates to form the pairs
#'   from, for the all-pairs and combined readings. By default every covariate of the model that is
#'   not constant is used.
#'
#' @return An object of class \code{"deepgof1"}: a list with \code{statistic} (the observed
#'   score), \code{p.value}, \code{B}, \code{K}, \code{reading}, \code{axes} (the
#'   covariates of the map that gave the statistic), \code{map} (that map, a \code{K} by \code{K}
#'   matrix of standardized residual sums whose rows follow the first axis), \code{pairs}
#'   (for the all-pairs and combined readings, the observed score of every pair), \code{components}
#'   (for \code{"combined"}, the p-values of the axis rule and the all-pairs reading), \code{boot} (the \code{B}
#'   bootstrap statistics of the reported statistic) and \code{method}.
#'
#' @section Reproducibility:
#' A bootstrap refit that fails to converge is scored \code{+Inf}, so it counts against
#' rejection -- the conservative direction. Set a seed before calling for a reproducible
#' p-value. When some covariate column has tied values, the random tie-breaking uses the
#' same seed, so from version 2.8.0 such calls give a different p-value for a given seed
#' than earlier versions did; calls without ties give the same p-value as before.
#' Version 2.9.0 refits each bootstrap sample on the fitted model's design matrix. Earlier
#' versions refitted the formula on the model frame, which fails for every term that
#' transforms a covariate (\code{log(x)}, \code{ns(x, 3)}, \code{poly(x, 2)}): each
#' replicate was then scored \code{+Inf} and the p-value was 1. For models without such
#' terms the default reading gives the same p-value as in 2.8.0.
#'
#' @references
#' Ebrahim EK (2026). "DeepGOF-1: A Pretrained Convolutional Goodness-of-Fit Test for
#' Logistic Regression with a Computable Consistency Certificate." Manuscript under review.
#' Reproduction materials and frozen weights: \doi{10.5281/zenodo.22113220}
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
#' # every pair of covariates, with the pair where the misfit is largest
#' deepgof1(fit, B = 49, reading = "allpairs")$axes
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
deepgof1 <- function(fit, B = 199L, K = 6L, reading = c("axes", "allpairs", "combined", "columns"),
                     covariates = NULL) {
  if (!inherits(fit, "glm") || fit$family$family != "binomial")
    stop("deepgof1() expects a glm fitted with family = binomial()", call. = FALSE)
  if (K != deepgof1_weights$K)
    stop("the shipped DeepGOF-1 weights are valid only at K = ", deepgof1_weights$K,
         call. = FALSE)
  ## the bootstrap draws one Bernoulli outcome per row, so grouped or weighted fits are not served
  if (any(fit$prior.weights != 1) || !all(fit$y %in% c(0, 1)))
    stop("deepgof1() needs a 0/1 outcome, one row per observation, without prior weights",
         call. = FALSE)
  reading <- match.arg(reading)
  if (!is.null(covariates) && reading %in% c("axes", "columns"))
    stop("'covariates' applies to the all-pairs and combined readings only", call. = FALSE)
  M <- deepgof1_weights
  score <- function(m) .dg_score((m - M$mu) / M$sd, M)

  mm <- stats::model.matrix(fit)
  tcols <- NULL
  if (reading != "columns") {
    X <- .dg_covariates(fit)
    if (is.null(X)) {
      ## the raw covariates are not recoverable (e.g. the data are gone): use the columns
      warning("deepgof1(): could not recover the covariates from the data the model was ",
              "fitted to; the maps are drawn over model-matrix columns instead", call. = FALSE)
      v <- colnames(mm)[colnames(mm) != "(Intercept)"]
      X <- mm[, v[!is.na(stats::coef(fit)[v])], drop = FALSE]
    } else {
      tcols <- .dg_term_cols(fit, mm, colnames(X))
    }
    if (!is.null(covariates)) {
      bad <- setdiff(covariates, colnames(X))
      if (length(bad))
        stop("'covariates' must name covariates of 'fit': ", paste(bad, collapse = ", "),
             call. = FALSE)
      X <- X[, covariates, drop = FALSE]
      if (!is.null(tcols)) tcols <- tcols[covariates]
    }
    if (ncol(X) < 1L) stop("deepgof1() needs at least one covariate", call. = FALSE)
  } else {
    v <- colnames(mm)[colnames(mm) != "(Intercept)"]
    if (length(v) < 1L) stop("deepgof1() needs at least one covariate", call. = FALSE)
    X <- mm[, v, drop = FALSE]
  }
  ## with one covariate there is no pair to choose: every reading is the 36-cell map along it
  pairs_path <- reading == "allpairs" || ncol(X) == 1L

  ## Ties among covariate values are broken at random, once per call: the same permutation
  ## serves the observed map and every bootstrap map, so the bootstrap reproduces the tie
  ## structure and the row order of the data cannot enter the statistic. Nothing is drawn
  ## when no covariate has a tie, so results for continuous covariates are unchanged.
  tied <- any(apply(X, 2, anyDuplicated) > 0L)
  tb <- if (tied) sample.int(nrow(mm)) else NULL

  ph  <- as.numeric(stats::fitted(fit))
  y   <- as.numeric(fit$y)                       # 0/1, also for a factor or logical response
  itc <- attr(stats::terms(fit), "intercept") > 0L

  if (reading == "columns" && !pairs_path) {
    ## The rule of versions 2.7.0 and 2.8.0, kept as it was for reproducing their results: one map
    ## over the two model-matrix columns with the largest |b| sd, chosen again in every replicate.
    obs  <- .dg_map(fit, K, tb)
    Sobs <- score(obs$map)
    axes <- obs$axes
    map  <- obs$map
    pairs <- NULL
    ## The refits use the fitted model's own design matrix. Refitting the formula on the model
    ## frame, as versions up to 2.8.0 did, fails for every term that transforms a covariate
    ## (log(x), ns(x, 3), poly(x, 2)), because the model frame holds the transformed column and
    ## not x; every replicate was then scored +Inf and the p-value was 1. The design matrix also
    ## keeps a spline basis fixed, which is what a parametric bootstrap under the fit requires.
    Sb <- numeric(B)
    for (b in seq_len(B)) {
      ys <- stats::rbinom(length(ph), 1L, ph)
      fb <- tryCatch(suppressWarnings(stats::glm.fit(mm, ys, weights = fit$prior.weights,
                                                     offset = fit$offset, family = fit$family,
                                                     control = fit$control, intercept = itc)),
                     error = function(e) NULL)
      ## a failed refit counts AGAINST rejection (conservative); -Inf would deflate p
      Sb[b] <- if (is.null(fb)) Inf else
        score(.dg_map_core(mm, ys, as.numeric(fb$fitted.values), fb$coefficients, K, tb)$map)
    }
    pval <- (1 + sum(Sb >= Sobs)) / (B + 1)
  } else {
    need <- switch(reading, axes = "axes", allpairs = "allpairs", combined = c("axes", "allpairs"),
                   columns = "allpairs")
    E <- .dg_engine(fit, B, K, X, mm, tb, score, need = need, tcols = tcols)
    P <- .dg_pvalues(E$obs$v, E$boot)
    key <- if (reading == "combined") "combined" else if (pairs_path) "allpairs" else "axes"
    pval <- P[[key]]
    Sobs <- if (reading == "combined") E$obs$v[["axes"]] else E$obs$v[[if (pairs_path) "allpairs" else "axes"]]
    Sb   <- E$boot[, if (key == "combined") "axes" else key]
    if (pairs_path) {
      top <- E$obs$top
      axes <- if (ncol(X) == 1L) colnames(X) else E$cells$pairs[top, ]
      map  <- .dg_cellmap(E$cells$idx[[top]], y - ph, ph, K)
    } else {
      axes <- E$obs$axes
      map  <- E$obs$map
    }
    pairs <- NULL
    if ("allpairs" %in% need && ncol(X) > 1L) {
      so <- attr(E$obs$v, "pairscores")
      pairs <- data.frame(axis1 = E$cells$pairs[, 1L], axis2 = E$cells$pairs[, 2L], score = so,
                          stringsAsFactors = FALSE)
    }
    if (reading == "combined")
      attr(pval, "components") <- c(axes = P[["axes"]], allpairs = P[["allpairs"]])
  }
  structure(list(statistic = Sobs,
                 p.value   = as.numeric(pval),
                 components = attr(pval, "components"),
                 B = B, K = K, reading = reading, axes = axes,
                 map = matrix(map, K, K, byrow = TRUE), pairs = pairs, boot = Sb,
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
  if (length(x$axes) == 1L) {
    cat(sprintf("map: %d quantile cells over ranks of %s\n", x$K * x$K, x$axes))
  } else if (identical(x$reading, "allpairs")) {
    cat(sprintf("map: %d x %d, maximum over %d pair%s of covariates, reached at %s and %s\n",
                x$K, x$K, nrow(x$pairs), if (nrow(x$pairs) == 1L) "" else "s",
                x$axes[1], x$axes[2]))
  } else {
    cat(sprintf("grid: %d x %d over ranks of %s and %s\n", x$K, x$K, x$axes[1], x$axes[2]))
  }
  if (identical(x$reading, "combined") && !is.null(x$components))
    cat(sprintf("combined reading: axis rule p = %.4f, all-pairs p = %.4f\n",
                x$components[["axes"]], x$components[["allpairs"]]))
  alt <- if (is.null(x$alternative)) "the logistic model is misspecified" else x$alternative
  cat("alternative hypothesis: ", alt, "\n\n", sep = "")
  invisible(x)
}

#' DeepGOF-1 for frozen predictions: external validation of a risk model
#'
#' Tests whether given probabilities are calibrated for given 0/1 outcomes within subgroups of
#' the covariates, for a model that was fitted elsewhere and is not refitted: a published risk
#' score checked on new patients, or the predictions of any model, logistic or not, on a
#' validation set. The null hypothesis is that each \eqn{y_i} is Bernoulli(\eqn{p_i}) with
#' \eqn{p_i} as given.
#'
#' The residual map of \code{\link{deepgof1}} is drawn over the ranks of the covariates and
#' standardized by the given probabilities, and the shipped network scores it. The reference
#' distribution is built by drawing \eqn{y^* \sim} Bernoulli(\eqn{p}) and scoring again; with
#' nothing estimated, the law of the score is the same under the whole null, so the rank
#' p-value is exactly valid at every sample size (Besag and Clifford 1989).
#'
#' Because the map lays the residuals out over the covariates, the test checks calibration
#' within patient subgroups (strong calibration in the hierarchy of Van Calster et al. 2016).
#' A wrong overall rate or a wrong calibration slope is the same for every patient and is
#' found better by tests along the predicted risk, such as the calibration belt
#' (\code{\link{run.all.gof}} with GiViTI); a miscalibration that differs between patients
#' with the same predicted risk, a missed U-shape, threshold or interaction, is found better
#' by this test, which also shows where it lies.
#'
#' The axis rule needs a ranking of the covariates by their effect. Without coefficients, it
#' regresses the logit of \code{p} on the covariates by least squares and scores each
#' covariate by \eqn{|b_j| \hat\sigma_j}; for a published logistic model in these covariates
#' this recovers its own coefficients exactly. Name the axes with \code{axes} to fix them in
#' advance.
#'
#' @param y 0/1 outcomes.
#' @param p the predicted probabilities to be checked, strictly between 0 and 1, one per
#'   outcome.
#' @param X a numeric matrix or data frame of covariates, one row per outcome, with column
#'   names. At least one column.
#' @param B number of Monte Carlo draws; the p-value lies on a grid of \code{1 / (B + 1)}.
#' @param K grid resolution; leave at 6, the resolution the shipped weights were trained at.
#' @param reading \code{"combined"} (the default), \code{"axes"} or \code{"allpairs"}, as in
#'   \code{\link{deepgof1}}. The combined reading is the exact minimum of the other two,
#'   calibrated over the same draws.
#' @param axes optional names of two columns of \code{X} for the axis-rule map, fixing it in
#'   advance.
#' @return An object of class \code{"deepgof1"}, as returned by \code{\link{deepgof1}}.
#' @references
#' Besag, J. and Clifford, P. (1989). Generalized Monte Carlo significance tests.
#' \emph{Biometrika}, 76(4), 633--642.
#'
#' Van Calster, B., Nieboer, D., Vergouwe, Y., De Cock, B., Pencina, M. J. and Steyerberg,
#' E. W. (2016). A calibration hierarchy for risk models was defined: from utopia to
#' empirical data. \emph{Journal of Clinical Epidemiology}, 74, 167--176.
#' @examples
#' set.seed(1)
#' n <- 400
#' X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
#' p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)          # the published model
#' y <- rbinom(n, 1, plogis(qlogis(p) + 0.8 * X$x1 * X$x2)) # the new patients
#' deepgof1.external(y, p, X, B = 99)
#' @seealso \code{\link{deepgof1}}
#' @concept external validation
#' @concept calibration
#' @export
deepgof1.external <- function(y, p, X, B = 199L, K = 6L, reading = c("combined", "axes", "allpairs"),
                              axes = NULL) {
  reading <- match.arg(reading)
  if (K != deepgof1_weights$K)
    stop("the shipped DeepGOF-1 weights are valid only at K = ", deepgof1_weights$K, call. = FALSE)
  y <- as.numeric(y); p <- as.numeric(p)
  if (!all(y %in% c(0, 1))) stop("'y' must be 0/1", call. = FALSE)
  if (length(p) != length(y) || anyNA(p) || any(p <= 0 | p >= 1))
    stop("'p' must hold one probability strictly between 0 and 1 per outcome", call. = FALSE)
  X <- as.data.frame(X)
  if (nrow(X) != length(y)) stop("'X' must have one row per outcome", call. = FALSE)
  if (!all(vapply(X, is.numeric, TRUE)))
    stop("'X' must be numeric; code a factor as numeric columns first", call. = FALSE)
  X <- as.matrix(X)
  if (ncol(X) < 1L || is.null(colnames(X)) || anyDuplicated(colnames(X)))
    stop("'X' needs at least one column, with distinct names", call. = FALSE)
  if (anyNA(X)) stop("'X' has missing values", call. = FALSE)
  if (!is.null(axes) && (length(axes) != 2L || !all(axes %in% colnames(X)) || axes[1] == axes[2]))
    stop("'axes' must name two different columns of 'X'", call. = FALSE)
  M <- deepgof1_weights
  score <- function(m) .dg_score((m - M$mu) / M$sd, M)
  n <- length(y)

  ## ties broken at random once, as in deepgof1(): the same order serves every draw
  tb <- if (any(apply(X, 2, anyDuplicated) > 0L)) sample.int(n) else NULL
  cells <- .dg_cells(X, K, tb)
  one <- ncol(X) == 1L
  if (one) {
    kax <- 1L
  } else {
    if (is.null(axes)) {
      ## the axis rule without coefficients: least squares of logit(p) on the covariates
      b <- stats::coef(stats::lm.fit(cbind(1, X), stats::qlogis(p)))[-1L]
      b[is.na(b)] <- 0
      sc <- abs(b) * apply(X, 2, stats::sd)
      axes <- colnames(X)[sort(order(sc, decreasing = TRUE)[1:2])]
    }
    ax <- colnames(X)[sort(match(axes, colnames(X)))]
    kax <- which(cells$pairs[, 1L] == ax[1] & cells$pairs[, 2L] == ax[2])
  }
  need_all <- reading != "axes" && !one
  stats_of <- function(yy) {
    e <- yy - p
    if (!need_all) {
      s <- score(.dg_cellmap(cells$idx[[kax]], e, p, K))
      return(c(axes = s, allpairs = s))
    }
    s <- vapply(cells$idx, function(i) score(.dg_cellmap(i, e, p, K)), 0)
    structure(c(axes = s[[kax]], allpairs = max(s)), pairscores = s)
  }
  obs <- stats_of(y)
  boot <- t(vapply(seq_len(B), function(b) as.numeric(stats_of(stats::rbinom(n, 1L, p))), c(0, 0)))
  colnames(boot) <- c("axes", "allpairs")
  p1 <- function(k) (1 + sum(boot[, k] >= obs[[k]])) / (B + 1)
  within <- function(v) vapply(seq_along(v), function(i) sum(v >= v[i]) / length(v), 0)
  mins <- pmin(within(c(obs[["axes"]], boot[, "axes"])), within(c(obs[["allpairs"]], boot[, "allpairs"])))
  pv <- c(axes = p1("axes"), allpairs = p1("allpairs"), combined = sum(mins <= mins[1]) / (B + 1))

  key <- if (one) "axes" else reading
  top <- if (need_all) which.max(attr(obs, "pairscores")) else kax
  show <- if (reading == "allpairs") top else kax
  pairs <- if (need_all) data.frame(axis1 = cells$pairs[, 1L], axis2 = cells$pairs[, 2L],
                                    score = attr(obs, "pairscores"), stringsAsFactors = FALSE) else NULL
  structure(list(statistic = obs[[if (reading == "allpairs") "allpairs" else "axes"]],
                 p.value   = as.numeric(pv[[key]]),
                 components = if (reading == "combined" && !one) pv[c("axes", "allpairs")] else NULL,
                 B = B, K = K, reading = if (one) "axes" else reading,
                 axes = if (one) colnames(X) else cells$pairs[show, ],
                 map = matrix(.dg_cellmap(cells$idx[[show]], y - p, p, K), K, K, byrow = TRUE),
                 pairs = pairs, boot = boot[, if (reading == "allpairs") "allpairs" else "axes"],
                 method = "DeepGOF-1 for frozen predictions: exact Monte Carlo calibration test",
                 data.name = paste(deparse(substitute(y)), "and", deparse(substitute(p))),
                 alternative = "the probabilities are miscalibrated within covariate subgroups"),
            class = "deepgof1")
}
