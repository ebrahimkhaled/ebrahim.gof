#' EDGE-Belt: the grouped calibration picture behind the EDGE test
#'
#' @description
#' Returns, and draws, the group-level quantities that \code{\link{edge.gof}}
#' (\code{\link{def.gof}}) computes and then discards: the standardized grouped
#' residuals, their estimation-adjusted null standard deviations, the shape the
#' test projects onto, and the bands that say which groups are out of line.
#' The p-value reported with the belt is the one \code{def.gof()} itself returns,
#' so the picture and the test cannot disagree.
#'
#' @details
#' The observations are sorted by predicted risk and split into \code{G}
#' equal-frequency groups, exactly as in \code{def.gof()}. For group \eqn{g} with
#' \eqn{O_g} events, \eqn{E_g = \sum \hat p_i} and
#' \eqn{V_g = \sum \hat p_i(1-\hat p_i)}, the standardized residual is
#' \eqn{r_g = (O_g - E_g)/\sqrt{V_g}}.
#'
#' Because \eqn{\beta} was estimated from the same data, \eqn{r} is not a vector of
#' independent standard normals: under the fitted null its covariance is
#' \eqn{\Omega = I - U(X'WX)^{-1}U'}, whose diagonal is \strong{below one}. The belt
#' divides each residual by its own \eqn{\sqrt{\Omega_{gg}}}, so a band drawn at
#' \eqn{\pm 1.96} is the correct pointwise band and not the too-wide naive one.
#'
#' \strong{Two bands, and why the wider one matters.} With \code{G} groups, a
#' pointwise 95\% band is expected to exclude about \code{0.05 * G} groups when the
#' model is perfect: at \code{G = 40} that is two groups flagged by chance alone.
#' The belt therefore also carries a \emph{simultaneous} band, the level-\code{level}
#' quantile of \eqn{\max_g |r_g|/\sqrt{\Omega_{gg}}} under \eqn{r \sim N(0, \Omega)},
#' obtained by simulation from the known \eqn{\Omega}. Read the simultaneous band
#' when asking "is this model miscalibrated anywhere"; read the pointwise band only
#' for a group singled out in advance.
#'
#' \strong{The curve.} The orange curve is \eqn{B(B'B)^{-1}B'r}, the component of the
#' residual vector that the EDGE statistic squares (\eqn{B = Z} for the unit form,
#' \eqn{B = Z\sqrt{V_g}} for the score form). It is the shape the test looked for,
#' drawn on the same axes as the residuals it was extracted from.
#'
#' \strong{What this is not.} The belt is \code{G} discrete groups, not a continuous
#' curve over the risk range: it is coarser than the calibration belt of Nattino et
#' al. (2014), implemented in \pkg{givitiR}, and it inherits whatever the
#' equal-frequency grouping hides. It is a diagnostic companion to the test, not a
#' replacement for a calibration curve.
#'
#' @param object A fitted binary logistic \code{\link[stats]{glm}}, or a binary
#'   (0/1) response vector (then supply \code{predicted_probs}).
#' @param predicted_probs Numeric predicted probabilities; required when
#'   \code{object} is not a glm.
#' @param X Optional design matrix, used for the estimation adjustment when
#'   \code{object} is not a glm. Without it the belt falls back to
#'   \eqn{\Omega = I}, which gives conservative (too wide) bands.
#' @param external \code{TRUE} when the predicted risks were frozen before these
#'   outcomes were seen, as in external validation of a model fitted elsewhere.
#'   Nothing was estimated from these data, so \eqn{\Omega = I} is then the
#'   \emph{exact} reference rather than a conservative fallback, and the belt says
#'   so instead of warning. The test is then \code{def.gof(..., external = TRUE)}:
#'   a column of ones joins the basis, so the curve also carries the overall level.
#'   Only \code{weights = "unit"} is available. Default \code{FALSE}.
#' @param G Number of equal-frequency groups, or \code{"auto"} for
#'   \code{max(10, ceiling(n/25))}. Default \code{"auto"}.
#' @param basis Shape basis: \code{"poly3"} (default), \code{"poly2"},
#'   \code{"stukel"} or \code{"sym"}.
#' @param weights \code{"unit"} (default, the statistic as published) or
#'   \code{"score"}.
#' @param method Reference for the p-value, passed to \code{def.gof()}.
#' @param level Band level. Default \code{0.95}.
#' @param nsim Draws used for the simultaneous band. Default \code{20000}.
#' @param seed Optional seed for that simulation, for a reproducible band.
#'
#' @return An object of class \code{"edge_belt"}: a list with \code{test} (the
#'   \code{def.gof()} result), \code{groups} (one row per group: size, mean
#'   predicted risk, observed and expected events, variance, residual, null SD,
#'   standardized residual, the fitted direction, and the two flags),
#'   \code{critical} (the two band half-widths and the level), and the settings.
#'
#' @references
#' Nattino, G., Finazzi, S. and Bertolini, G. (2014). A new calibration test and a
#' reappraisal of the calibration belt for the assessment of prediction models based
#' on dichotomous outcomes. \emph{Statistics in Medicine} \strong{33}, 2390-2407.
#'
#' @examples
#' set.seed(1)
#' n  <- 800
#' x  <- runif(n, -3, 3)
#' y  <- rbinom(n, 1, 1 - exp(-exp(0.8 * x)))   # log-log truth, fitted as logit
#' fit <- glm(y ~ x, family = binomial())
#' b <- edge.belt(fit)
#' b
#' plot(b)
#'
#' @seealso \code{\link{edge.gof}}, \code{\link{def.gof}}.
#' @importFrom stats fitted predict model.matrix qlogis poly quantile rnorm qnorm
#' @concept calibration
#' @concept model diagnostics
#' @concept EDGE
#' @export
edge.belt <- function(object, predicted_probs = NULL, X = NULL, G = "auto",
                      basis   = c("poly3", "poly2", "stukel", "sym"),
                      weights = c("unit", "score"),
                      method  = c("satterthwaite", "imhof"),
                      external = FALSE, level = 0.95, nsim = 20000L, seed = NULL) {

  basis   <- match.arg(basis)
  weights <- match.arg(weights)
  method  <- match.arg(method)
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1)
    stop("'level' must be a single number in (0, 1).")
  if (external) {
    if (inherits(object, "glm"))
      stop("external = TRUE needs the frozen predictions, not the fit that produced them: ",
           "pass the validation outcomes as 'object' and their predicted risks as 'predicted_probs'.")
    X <- NULL                       # nothing was estimated here, so there is nothing to adjust for
  }

  ## the official numbers come from the tested function itself, never recomputed here
  tst <- def.gof(object, predicted_probs = predicted_probs, X = X, G = G,
                 basis = basis, method = method, weights = weights, external = external)

  p <- suppressWarnings(.edge_belt_pieces(object, predicted_probs, X, G, basis, external))

  sd_g <- sqrt(pmax(diag(p$Omega), 0))
  ok   <- sd_g > 1e-8
  z    <- rep(NA_real_, p$Gn)
  z[ok] <- p$r[ok] / sd_g[ok]

  ## the shape the statistic squares, in the form the user asked for
  B   <- if (weights == "score") p$Z * sqrt(p$Vg) else p$Z
  fit <- drop(B %*% solve(crossprod(B), crossprod(B, p$r)))

  crit_point <- stats::qnorm(1 - (1 - level) / 2)
  crit_sim   <- .edge_belt_simband(p$Omega, sd_g, ok, level, nsim, seed)

  ng  <- vapply(p$idx, length, integer(1))
  sdo <- sqrt(p$Vg * pmax(diag(p$Omega), 0)) / ng      # the same band, on the risk scale

  groups <- data.frame(
    group      = seq_len(p$Gn),
    size       = as.integer(ng),
    risk_mean  = p$pbar,
    observed   = p$og,
    expected   = p$eg,
    variance   = p$Vg,
    residual   = p$r,
    sd_null    = sd_g,
    z          = z,
    direction  = fit,
    rate_obs   = p$og / ng,
    rate_sd    = sdo,
    outside_pointwise    = is.finite(z) & abs(z) > crit_point,
    outside_simultaneous = is.finite(z) & abs(z) > crit_sim,
    row.names  = NULL, stringsAsFactors = FALSE)

  structure(list(test = tst, groups = groups,
                 critical = list(pointwise = crit_point, simultaneous = crit_sim,
                                 level = level),
                 basis = basis, weights = weights, G = p$Gn, n = p$n,
                 adjusted = !p$naive, external = external, call = match.call()),
            class = "edge_belt")
}

## Internal: the grouped pieces, computed exactly as def.gof() computes them. Kept
## verbatim rather than refactored out of def.gof(), so nothing the battery depends on
## changes; tests/testthat/test-edge-belt.R checks that the two agree to machine precision.
.edge_belt_pieces <- function(object, predicted_probs, X, G, basis, external = FALSE) {
  if (inherits(object, "glm")) {
    if (object$family$family != "binomial")
      stop("'object' must be a binomial glm (or pass y, predicted_probs, X).")
    y   <- as.numeric(object$y)
    ph  <- pmin(pmax(as.numeric(object$fitted.values), 1e-6), 1 - 1e-6)
    eta <- as.numeric(object$linear.predictors)
    dmu <- object$family$mu.eta(eta)
    X   <- stats::model.matrix(object)
    naive <- FALSE
  } else {
    if (!is.numeric(object))
      stop("'object' must be a fitted binomial glm or a numeric (0/1) y vector.")
    y <- as.numeric(object)
    if (is.null(predicted_probs))
      stop("Provide 'predicted_probs' when 'object' is not a glm.")
    ph  <- pmin(pmax(as.numeric(predicted_probs), 1e-6), 1 - 1e-6)
    dmu <- ph * (1 - ph)
    naive <- is.null(X)
  }
  n <- length(y)
  if (!all(y %in% c(0, 1))) stop("EDGE needs a binary (0/1) response.")
  if (identical(G, "auto")) G <- .def_auto_G(n)
  if (G > n) stop("'G' cannot exceed the number of observations.")

  V <- ph * (1 - ph)
  w <- dmu^2 / V

  grp  <- pmin(ceiling(rank(ph, ties.method = "first") / (n / G)), G)
  idx  <- split(seq_len(n), grp)
  Gn   <- length(idx)
  og   <- vapply(idx, function(I) sum(y[I]),   numeric(1))
  eg   <- vapply(idx, function(I) sum(ph[I]),  numeric(1))
  Vg   <- vapply(idx, function(I) sum(V[I]),   numeric(1))
  pbar <- vapply(idx, function(I) mean(ph[I]), numeric(1))
  r    <- (og - eg) / sqrt(Vg)

  if (naive) {
    Omega <- diag(Gn)
  } else {
    U     <- t(vapply(idx, function(I) colSums(dmu[I] * X[I, , drop = FALSE]),
                      numeric(ncol(X)))) / sqrt(Vg)
    Omega <- diag(Gn) - U %*% solve(crossprod(X, w * X)) %*% t(U)
  }

  ## frozen predictions: the constant column joins the basis, reduced to full rank, as in the
  ## external mode of def.gof()
  Z <- if (external) cbind(1, .def_basis(pbar, basis)) else .def_basis(pbar, basis)
  Z <- Z[, colSums(abs(Z)) > 1e-8, drop = FALSE]
  if (ncol(Z) < 1)
    stop("The chosen basis is degenerate for this fit. Try basis = 'poly3' or a larger G.")
  Z <- Z / rep(sqrt(colSums(Z^2)), each = nrow(Z))
  if (external) {
    Q <- qr(Z)
    Z <- Z[, Q$pivot[seq_len(Q$rank)], drop = FALSE]
  }

  list(n = n, Gn = Gn, idx = idx, og = og, eg = eg, Vg = Vg, pbar = pbar,
       r = r, Omega = Omega, Z = Z, naive = naive)
}

## Internal: the simultaneous band. Under the fitted null r ~ N(0, Omega) to the same
## order as the test's own reference, so the band is the level quantile of the largest
## standardized residual, taken from draws of that exact law.
.edge_belt_simband <- function(Omega, sd_g, ok, level, nsim, seed) {
  if (!any(ok)) return(NA_real_)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  ev  <- eigen((Omega + t(Omega)) / 2, symmetric = TRUE)
  pos <- ev$values > 1e-10
  if (!any(pos)) return(NA_real_)
  L   <- ev$vectors[, pos, drop = FALSE] %*% diag(sqrt(ev$values[pos]), sum(pos))
  sim <- L %*% matrix(stats::rnorm(sum(pos) * nsim), sum(pos), nsim)
  m   <- apply(abs(sim[ok, , drop = FALSE] / sd_g[ok]), 2, max)
  as.numeric(stats::quantile(m, level, names = FALSE))
}

#' @param x An \code{"edge_belt"} object.
#' @param ... Ignored.
#' @rdname edge.belt
#' @export
print.edge_belt <- function(x, ...) {
  cat("EDGE-Belt\n")
  cat(sprintf("  n = %d, G = %d groups, basis = %s, form = %s%s\n",
              x$n, x$G, x$basis, x$weights,
              if (x$adjusted) "" else if (isTRUE(x$external))
                "  (frozen predictions: Omega = I is exact)" else
                "  (Omega = I: bands are conservative)"))
  cat(sprintf("  statistic = %.4g on df = %.3g,  p = %.4g\n",
              x$test$Test_Statistic, x$test$df, x$test$p_value))
  cat(sprintf("  bands at %g%%: pointwise +/-%.2f, simultaneous +/-%.2f\n",
              100 * x$critical$level, x$critical$pointwise, x$critical$simultaneous))
  ns <- sum(x$groups$outside_simultaneous)
  np <- sum(x$groups$outside_pointwise)
  cat(sprintf("  groups outside: %d simultaneous, %d pointwise (%.1f expected by chance)\n",
              ns, np, (1 - x$critical$level) * x$G))
  if (ns > 0) {
    g <- x$groups[x$groups$outside_simultaneous, ]
    cat("  flagged groups (simultaneous):\n")
    for (i in seq_len(nrow(g)))
      cat(sprintf("    group %2d  mean risk %.3f  observed %d vs expected %.1f  z = %+.2f\n",
                  g$group[i], g$risk_mean[i], as.integer(g$observed[i]), g$expected[i], g$z[i]))
  }
  invisible(x)
}

#' @param which \code{1} for the residual panel, \code{2} for the risk panel,
#'   \code{c(1, 2)} (default) for both.
#' @param main Optional title; the default names the basis, the grouping and the
#'   p-value.
#' @param palette Optional named character vector overriding the colours
#'   (\code{band}, \code{point_band}, \code{pt}, \code{flag}, \code{curve},
#'   \code{axis}, \code{grid}).
#' @rdname edge.belt
#' @importFrom graphics par plot polygon lines points abline text mtext legend box rect segments
#' @importFrom grDevices adjustcolor
#' @export
plot.edge_belt <- function(x, which = c(1, 2), main = NULL, palette = NULL, ...) {
  col <- c(band = "#cfe0f3", point_band = "#8fa3bb", pt = "#22384f", flag = "#c0392b",
           curve = "#e07b39", axis = "#4a5568", grid = "#e8ecf1")
  if (!is.null(palette)) col[names(palette)] <- palette
  g  <- x$groups
  cs <- x$critical$simultaneous
  cp <- x$critical$pointwise

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(1, length(which)), mar = c(4.4, 4.4, 3.4, 1.2),
                oma = c(2.6, 0, 2.2, 0), mgp = c(2.5, 0.7, 0), tcl = -0.3, las = 1,
                col.axis = col[["axis"]], col.lab = col[["axis"]], font.lab = 2,
                cex.axis = 0.9)

  if (1 %in% which) {
    yl <- max(cs * 1.25, max(abs(g$z), na.rm = TRUE) * 1.1, 3)
    plot(g$risk_mean, g$z, type = "n", ylim = c(-yl, yl),
         xlab = "mean predicted risk in the group", ylab = "standardized residual")
    graphics::abline(h = pretty(c(-yl, yl)), col = col[["grid"]], lwd = 1)
    graphics::rect(graphics::par("usr")[1], -cs, graphics::par("usr")[2], cs,
                   col = grDevices::adjustcolor(col[["band"]], 0.55), border = NA)
    graphics::abline(h = c(-cp, cp), lty = 2, col = col[["point_band"]], lwd = 1.3)
    graphics::abline(h = 0, col = col[["axis"]], lwd = 1.2)
    o <- order(g$risk_mean)
    graphics::lines(g$risk_mean[o], (g$direction / g$sd_null)[o], col = col[["curve"]], lwd = 2.6)
    fl <- g$outside_simultaneous
    graphics::points(g$risk_mean[!fl], g$z[!fl], pch = 21, bg = col[["pt"]],
                     col = "white", cex = 1.25, lwd = 1.1)
    if (any(fl)) {
      graphics::points(g$risk_mean[fl], g$z[fl], pch = 21, bg = col[["flag"]],
                       col = "white", cex = 1.7, lwd = 1.2)
      # name only the worst few, above a positive residual and below a negative one,
      # so a run of neighbouring flagged groups does not overprint itself
      lab <- which(fl)[order(-abs(g$z[fl]))]
      lab <- lab[seq_len(min(6L, length(lab)))]
      graphics::text(g$risk_mean[lab], g$z[lab], labels = g$group[lab],
                     pos = ifelse(g$z[lab] > 0, 3, 1), offset = 0.55,
                     cex = 0.72, col = col[["flag"]], font = 2)
    }
    graphics::box(col = col[["axis"]])
    graphics::mtext("group residuals and the shape EDGE tested", side = 3, line = 1.0,
                    adj = 0, font = 2, cex = 1.0, col = col[["pt"]])
  }

  if (2 %in% which) {
    o  <- order(g$risk_mean)
    lo <- (g$risk_mean - cs * g$rate_sd)[o]
    hi <- (g$risk_mean + cs * g$rate_sd)[o]
    rg <- range(c(lo, hi, g$rate_obs), na.rm = TRUE)
    plot(g$risk_mean, g$rate_obs, type = "n", xlim = rg, ylim = rg,
         xlab = "predicted risk", ylab = "observed rate")
    graphics::abline(h = pretty(rg), v = pretty(rg), col = col[["grid"]], lwd = 1)
    graphics::polygon(c(g$risk_mean[o], rev(g$risk_mean[o])), c(lo, rev(hi)),
                      col = grDevices::adjustcolor(col[["band"]], 0.55), border = NA)
    graphics::abline(0, 1, col = col[["axis"]], lwd = 1.2)
    fl <- g$outside_simultaneous
    graphics::segments(g$risk_mean, g$risk_mean, g$risk_mean, g$rate_obs,
                       col = grDevices::adjustcolor(col[["pt"]], 0.35), lwd = 1)
    graphics::points(g$risk_mean[!fl], g$rate_obs[!fl], pch = 21, bg = col[["pt"]],
                     col = "white", cex = 1.25, lwd = 1.1)
    if (any(fl))
      graphics::points(g$risk_mean[fl], g$rate_obs[fl], pch = 21, bg = col[["flag"]],
                       col = "white", cex = 1.7, lwd = 1.2)
    graphics::box(col = col[["axis"]])
    graphics::mtext("observed against predicted, same band", side = 3, line = 1.0,
                    adj = 0, font = 2, cex = 1.0, col = col[["pt"]])
  }

  ttl <- if (is.null(main))
    sprintf("EDGE-Belt   basis %s, G = %d, n = %d,  p = %.3g",
            x$basis, x$G, x$n, x$test$p_value) else main
  graphics::mtext(ttl, side = 3, line = 0.4, outer = TRUE, font = 2, cex = 1.15,
                  col = col[["pt"]])

  ## one key for the whole figure, in the outer margin, where nothing can collide with it
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  graphics::legend("bottom", horiz = TRUE, bty = "n", cex = 0.8, seg.len = 1.6,
                   legend = c(sprintf("simultaneous %g%%", 100 * x$critical$level),
                              sprintf("pointwise %g%%", 100 * x$critical$level),
                              "EDGE direction", "flagged group"),
                   fill = c(grDevices::adjustcolor(col[["band"]], 0.55), NA, NA, NA),
                   border = NA, lty = c(NA, 2, 1, NA), lwd = c(NA, 1.3, 2.6, NA),
                   pch = c(NA, NA, NA, 21), pt.bg = c(NA, NA, NA, col[["flag"]]),
                   pt.cex = c(NA, NA, NA, 1.4),
                   col = c(NA, col[["point_band"]], col[["curve"]], "white"),
                   text.col = col[["axis"]])
  invisible(x)
}
