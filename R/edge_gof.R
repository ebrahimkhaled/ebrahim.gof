#' EDGE: Directed Goodness-of-Fit Test for Binary Logistic Regression
#'
#' @description
#' \code{edge.gof()} is the primary interface to the EDGE test (Efficient Directed
#' Grouped Examination): a grouped, \emph{directed} goodness-of-fit test
#' for binary logistic regression under sparse data. EDGE projects the grouped
#' standardized residuals onto a small pre-specified basis of calibration shapes
#' (cubic \code{"poly3"} by default) and refers the resulting quadratic form to
#' its closed-form weighted chi-squared null distribution -- no refit, no
#' resampling, no tuning.
#'
#' \code{edge.gof()} computes exactly the same statistic as the legacy name
#' \code{\link{def.gof}} (retained for backward compatibility); the returned
#' \code{Test} label is \code{"EDGE"}.
#'
#' @details
#' \strong{Two partitions, side by side.} By default the test is reported at two
#' partitions, as the EDGE paper recommends. The \emph{default partition},
#' \code{G = "auto"} (\eqn{\max(10, \lceil n/25 \rceil)} groups), refines with the sample and has
#' more power, above all against misfit at the extremes of risk. \emph{Ten groups},
#' \code{G = 10}, keep many records in the extreme groups and so tolerate more
#' corrupted records, at some cost in power. Choose which one decides in the analysis
#' plan, from what is known about how the data were collected, not after seeing the
#' results. The first row is marked \code{Role = "verdict"} and the second
#' \code{Role = "check"}: with the default \code{G} the default partition decides and
#' ten groups are the robustness check; to let ten groups decide, give
#' \code{G = c(10, "auto")}. If only the default partition rejects, the signal sits in the
#' extreme groups: check those records. Give a single \code{G} to get one row.
#'
#' \strong{The picture.} \code{plot()} on the result draws the grouped residuals of one row
#' with the shape the test looked for and the bands that say which groups are out of line;
#' see \code{\link{plot.edge_gof}}.
#'
#' @inheritParams def.gof
#' @param G Number of groups: \code{"auto"}, a number, or a vector of them; the
#'   default \code{c("auto", 10)} reports both partitions.
#' @param y Optional alias for \code{object} when frozen predictions are tested:
#'   \code{edge.gof(y = y, predicted_probs = p, external = TRUE)}.
#'
#' @return A \code{data.frame} of class \code{c("edge_gof", "data.frame")} with one row per
#'   partition and columns \code{Test} (\code{"EDGE"}), \code{Partition}, \code{Role}
#'   (\code{"verdict"} for the first row, \code{"check"} for the second), \code{Basis},
#'   \code{Test_Statistic}, \code{df}, \code{Method} and \code{p_value}, as documented in
#'   \code{\link{def.gof}}. It is a data frame in every respect; the class only adds a
#'   \code{plot} method. The grouped pieces behind each row are kept in
#'   \code{attr(, "details")}, a list with one element per row holding \code{r} (the grouped
#'   residuals \eqn{(O_g - E_g)/\sqrt{V_g}}), \code{Pz_r} (their projection on the basis, whose
#'   sum of squares is the unit-form statistic), \code{O}, \code{E}, \code{V}, \code{size},
#'   \code{pbar} (the mean predicted risk of each group), \code{G}, \code{omega_diag} (the null
#'   variance of each residual: below one in-sample, one for frozen predictions),
#'   \code{external}, \code{adjusted}, \code{role}, \code{partition} and \code{statistic}.
#'   A row without a p-value, or with \code{basis = "ensemble"}, has \code{NULL} there.
#'   Rows selected or reordered with \code{x[i, ]} keep the attribute whole, and
#'   \code{plot()} finds the pieces of a row by its partition, role and statistic.
#'   \code{subset()} and \code{merge()} drop the attribute, and rows bound from another
#'   result have no pieces of their own; \code{plot()} then stops with a message.
#'
#' @references
#' Ebrahim EK, Hussein OAE-A, El-Kotory A (2026). "A Grouped Calibration Test for Logistic
#' Regression That Tolerates a Few Corrupted Records." arXiv:2608.20511 [stat.ME].
#' \doi{10.48550/arXiv.2608.20511}
#'
#' Ebrahim EK, Khattab IG, El-Kotory A (2026). "A Modified Hosmer-Lemeshow
#' Goodness-of-Fit Test for Asymmetric Links: Second-Order Power and Robustness."
#' arXiv:2607.15454 [stat.ME]. \doi{10.48550/arXiv.2607.15454}
#'
#' @author Ebrahim Khaled Ebrahim \email{ebrahimkhaled@@alexu.edu.eg}
#'
#' @examples
#' set.seed(1)
#' x <- runif(500, -3, 3)
#' y <- rbinom(500, 1, plogis(0.6 * x))
#' fit <- glm(y ~ x, family = binomial())
#' edge.gof(fit)                      # cubic basis, at the default partition and at ten groups
#' edge.gof(fit, G = 10)              # one partition only
#' edge.gof(fit, basis = "stukel")    # Stukel-shape basis
#' edge.gof(fit, basis = "sym")       # Stukel's symmetric direction, one column
#' edge.gof(fit, basis = "sym", weights = "score")   # its score form
#' edge.gof(fit, G = "auto")          # the default partition alone: max(10, ceiling(n / 25)) = 20 here
#' plot(edge.gof(fit))                # the verdict row, on the residual scale
#'
#' @seealso \code{\link{plot.edge_gof}}, \code{\link{def.gof}} (legacy name),
#'   \code{\link{edge.belt}}, \code{\link{ef.gof}}, \code{\link{def.ensemble.gof}},
#'   \code{\link{run.all.gof}}.
#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept EDGE
#' @concept directed test
#' @concept sparse data
#' @export
edge.gof <- function(object, predicted_probs = NULL, X = NULL, G = c("auto", 10),
                     basis = "poly3", method = "satterthwaite", weights = "unit",
                     external = FALSE, y = NULL) {
  if (missing(object)) {
    if (is.null(y)) stop("supply a fitted glm or the outcome vector")
    object <- y
  }
  n <- if (inherits(object, "glm")) length(object$y) else length(object)
  one <- function(g) {
    out <- .def_gof(object, predicted_probs = predicted_probs, X = X, G = g,
                    basis = basis, method = method, weights = weights, external = external,
                    keep = TRUE)
    det <- attr(out, "details")
    attr(out, "details") <- NULL
    out$Test <- "EDGE"
    gn <- if (identical(g, "auto")) .def_auto_G(n) else as.numeric(g)
    out$Partition <- if (identical(g, "auto")) sprintf("default (G = %d)", as.integer(gn)) else
      if (gn == 10) "ten groups (G = 10)" else sprintf("G = %d", as.integer(gn))
    list(row = out[, c("Test", "Partition", setdiff(names(out), c("Test", "Partition")))],
         details = det)
  }
  Gs  <- lapply(as.list(G), function(g) if (identical(g, "auto")) g else as.numeric(g))
  ## when n <= 250 the default partition is itself ten groups: report it once. The partitions
  ## are compared by their number of groups, so a row without a statistic (basis = "ensemble")
  ## is not taken for a repeat.
  gn  <- vapply(Gs, function(g) if (identical(g, "auto")) .def_auto_G(n) else g, numeric(1))
  Gs  <- Gs[!duplicated(gn)]
  res <- lapply(Gs, one)
  out <- do.call(rbind, lapply(res, `[[`, "row"))
  det <- lapply(res, `[[`, "details")
  out$Role <- c("verdict", "check", rep("check", max(0L, nrow(out) - 2L)))[seq_len(nrow(out))]
  out <- out[, c("Test", "Partition", "Role", setdiff(names(out), c("Test", "Partition", "Role")))]
  rownames(out) <- NULL
  det <- lapply(seq_along(det), function(i)
    if (is.null(det[[i]])) NULL else c(det[[i]], role = out$Role[i], partition = out$Partition[i],
                                       statistic = out$Test_Statistic[i]))
  attr(out, "details") <- det
  class(out) <- c("edge_gof", "data.frame")
  out
}

#' Plot an EDGE result: the grouped residuals, the shape tested and the bands
#'
#' Draws one row of an \code{\link{edge.gof}} result: the standardized residual of each risk
#' group against its mean predicted risk, the part of the residuals that the test squared
#' (\eqn{P_Z r}, the EDGE direction) as a curve, a pointwise band and a simultaneous band.
#' Groups beyond the band are drawn as filled triangles.
#'
#' @details
#' For group \eqn{g} with \eqn{O_g} events, \eqn{E_g = \sum \hat p_i} and
#' \eqn{V_g = \sum \hat p_i(1 - \hat p_i)}, the residual is
#' \eqn{r_g = (O_g - E_g)/\sqrt{V_g}}.
#'
#' \strong{Frozen predictions} (\code{external = TRUE}). Nothing was estimated from these
#' outcomes, so the residuals are independent with unit variance under a calibrated model and
#' the bands are the ones of the EDGE paper: pointwise at \eqn{\pm z_{1-(1-\mathrm{level})/2}}
#' and simultaneous at \eqn{\pm z_{1-(1-\mathrm{level}^{1/G})/2}}, the Sidak form of the
#' Bonferroni band, which has exact coverage \code{level} for \eqn{G} independent normal
#' residuals.
#'
#' \strong{In-sample} (a fitted glm, or \code{X} given). The fit removes part of every
#' residual, so \eqn{r} has covariance \eqn{\Omega = I - U(X'WX)^{-1}U'}, not \eqn{I}. Each
#' residual and the curve are divided by \eqn{\sqrt{\Omega_{gg}}} before the same constants
#' are applied. The pointwise band is then right to first order, and the simultaneous band is
#' conservative for correlated normal residuals (Sidak's inequality). Both rest on the normal
#' approximation and on an estimated \eqn{\Omega}, so the legend calls them approximate. When
#' only \code{y} and \code{predicted_probs} were given, without \code{X}, \eqn{\Omega = I} is
#' used and the bands are conservative.
#'
#' \strong{The rate scale.} \code{scale = "rate"} shows the observed event rate \eqn{O_g/n_g}
#' against the predicted rate \eqn{E_g/n_g} with the same band around the diagonal, half-width
#' times \eqn{\sqrt{V_g \Omega_{gg}}/n_g}, cut at 0 and 1. The bars join each group to the
#' diagonal: vermilion when more events were seen than predicted, blue when fewer.
#'
#' The band is a description of where the misfit lies. The verdict is the p-value of the row,
#' which does not depend on the band.
#'
#' @param x An \code{"edge_gof"} object from \code{\link{edge.gof}}.
#' @param row Which row to draw: \code{"verdict"} (default), \code{"check"}, or a row number.
#' @param scale \code{"residual"} (default), \code{"rate"}, or \code{c("residual", "rate")}
#'   for both panels side by side.
#' @param band \code{"simultaneous"} (default; the pointwise band is then drawn dashed),
#'   \code{"pointwise"} or \code{"none"}. Groups beyond the chosen band are marked.
#' @param colour \code{FALSE} draws in greys only.
#' @param level Coverage of the bands. Default \code{0.95}.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @examples
#' set.seed(1)
#' n <- 1000
#' x <- runif(n, -3, 3)
#' y <- rbinom(n, 1, 1 - exp(-exp(0.8 * x)))   # log-log truth, fitted as logit
#' fit <- glm(y ~ x, family = binomial())
#' e <- edge.gof(fit)
#' plot(e)
#' plot(e, row = "check", scale = c("residual", "rate"))
#'
#' ## frozen predictions: the bands are exact
#' xv <- runif(n, -3, 3)
#' yv <- rbinom(n, 1, plogis(0.3 + 0.6 * xv))
#' pv <- plogis(0.6 * xv)                       # the model misses the intercept shift
#' plot(edge.gof(yv, predicted_probs = pv, external = TRUE), scale = "rate", colour = FALSE)
#' @seealso \code{\link{edge.gof}}; \code{\link{edge.belt}} returns the same quantities as a
#'   table, with a simulated simultaneous band.
#' @importFrom graphics axis axTicks layout legend lines par plot points polygon rect segments abline box
#' @importFrom grDevices adjustcolor
#' @importFrom stats qnorm
#' @exportS3Method plot edge_gof
plot.edge_gof <- function(x, row = c("verdict", "check"), scale = c("residual", "rate"),
                          band = c("simultaneous", "pointwise", "none"), colour = TRUE,
                          level = 0.95, ...) {
  det <- attr(x, "details")
  if (is.null(det) || !nrow(x))
    stop("no grouped details on this object: plot the result of edge.gof(), or rows of it",
         call. = FALSE)
  if (is.numeric(row)) {
    k <- as.integer(row[1])
    if (is.na(k) || k < 1L || k > nrow(x)) stop("'row' must be between 1 and ", nrow(x), call. = FALSE)
  } else {
    asked <- !missing(row)
    row <- match.arg(row)
    k <- match(row, x$Role)
    if (is.na(k) && !asked && nrow(x) == 1L) {
      k <- 1L                              # one row taken with x[i, ]: draw it
    } else if (is.na(k)) stop(sprintf("this result has no '%s' row; it has: %s", row,
                                      paste(x$Role, collapse = ", ")), call. = FALSE)
  }
  ## A subset or a reordering of the rows keeps the attribute as it was, so the pieces are
  ## found by the row's partition, role and statistic, not by its position.
  d <- Filter(function(e) !is.null(e) && identical(e$partition, x$Partition[k]) &&
                identical(e$role, x$Role[k]) && isTRUE(all.equal(e$statistic, x$Test_Statistic[k])),
              det)
  if (!length(d))
    stop("this row has no grouped details: it has no p-value, has basis = 'ensemble', ",
         "or was bound from another result", call. = FALSE)
  d <- d[[1]]
  scale <- if (missing(scale)) "residual" else match.arg(scale, several.ok = TRUE)
  band  <- match.arg(band)
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1)
    stop("'level' must be a single number in (0, 1).", call. = FALSE)

  pal <- if (colour)
    list(band = "#E0E0E0", pw = "#737373", ink = "#404040", pt = "#404040", curve = "#D55E00",
         flag = "#0072B2", above = "#D55E00", below = "#0072B2", diag = "#8C8C8C")
  else
    list(band = "#E3E3E3", pw = "#6E6E6E", ink = "#262626", pt = "#8C8C8C", curve = "#000000",
         flag = "#000000", above = "#9E9E9E", below = "#9E9E9E", diag = "#8C8C8C")

  sdn <- sqrt(d$omega_diag)
  ok  <- sdn > 1e-8
  z   <- crv <- rep(NA_real_, d$G)
  z[ok]   <- d$r[ok] / sdn[ok]
  crv[ok] <- d$Pz_r[ok] / sdn[ok]
  cp  <- stats::qnorm(1 - (1 - level) / 2)
  cs  <- stats::qnorm(1 - (1 - level^(1 / max(1L, sum(ok)))) / 2)
  crit <- switch(band, simultaneous = cs, pointwise = cp, none = Inf)
  flag <- is.finite(z) & abs(z) > crit
  kind <- if (isTRUE(d$external)) "" else if (isTRUE(d$adjusted)) " (approximate)" else " (conservative)"
  xlab <- if (d$G == 10) c("Mean predicted risk in decile", "Mean predicted risk") else
    c(sprintf("Mean predicted risk in group (%d groups)", as.integer(d$G)),
      sprintf("Mean predicted risk (%d groups)", as.integer(d$G)), "Mean predicted risk")
  o    <- order(d$pbar)
  xr   <- range(d$pbar) + c(-1, 1) * 0.03 * diff(range(d$pbar))
  pct  <- function(side) {
    at <- graphics::axTicks(side)
    graphics::axis(side, at = at, labels = paste0(format(100 * at, trim = TRUE), "%"))
  }
  dot <- function(xx, yy) {
    graphics::points(xx[!flag], yy[!flag], pch = 21, bg = pal$pt, col = "white", cex = 1.15, lwd = 0.8)
    if (any(flag))
      graphics::points(xx[flag], yy[flag], pch = 24, bg = pal$flag, col = "white", cex = 1.3, lwd = 0.8)
  }

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::layout(matrix(seq_along(scale), 1))
  graphics::par(cex = 1)                     # layout() shrinks the text in a multi-panel figure
  nkey <- if (band == "none") 2L else if (band == "simultaneous") 5L else 4L
  nkey <- nkey - !("residual" %in% scale)
  graphics::par(mar = c(3.4, 3.7, 0.6, 0.8), oma = c(0.7 + 0.95 * ceiling(nkey / 2), 0, 0, 0),
                mgp = c(2.2, 0.55, 0), tcl = -0.3, las = 1, cex.axis = 0.85, cex.lab = 0.95)

  if ("residual" %in% scale) {
    yl <- max(3, 1.12 * max(abs(c(z, crv)), na.rm = TRUE), if (band == "none") 0 else 1.12 * crit)
    graphics::plot(NA, xlim = xr, ylim = c(-yl, yl), xaxt = "n", xlab = "",
                   ylab = "Standardized residual")
    graphics::title(xlab = .fit_label(xlab, 1))
    pct(1)
    u <- graphics::par("usr")
    if (band != "none") graphics::rect(u[1], -crit, u[2], crit, col = pal$band, border = NA)
    if (band == "simultaneous") graphics::abline(h = c(-cp, cp), lty = "22", col = pal$pw, lwd = 0.9)
    graphics::abline(h = 0, col = pal$ink, lwd = 0.9)
    graphics::lines(d$pbar[o], crv[o], col = pal$curve, lwd = 2)
    dot(d$pbar, z)
    graphics::box(col = pal$ink)
  }

  if ("rate" %in% scale) {
    ro <- d$O / d$size
    hw <- if (band == "none") 0 else crit * sqrt(d$V * d$omega_diag) / d$size
    lo <- pmax(0, d$pbar - hw)
    hi <- pmin(1, d$pbar + hw)
    yr <- range(c(ro, lo, hi, d$pbar), na.rm = TRUE)
    graphics::plot(NA, xlim = xr, ylim = yr + c(-1, 1) * 0.03 * diff(yr), xaxt = "n", yaxt = "n",
                   xlab = "", ylab = "Observed event rate")
    graphics::title(xlab = .fit_label(xlab, 1))
    pct(1); pct(2)
    if (band != "none")
      graphics::polygon(c(d$pbar[o], rev(d$pbar[o])), c(lo[o], rev(hi[o])), col = pal$band, border = NA)
    if (band == "simultaneous") {
      hp <- cp * sqrt(d$V * d$omega_diag) / d$size
      graphics::lines(d$pbar[o], pmax(0, d$pbar - hp)[o], lty = "22", col = pal$pw, lwd = 0.9)
      graphics::lines(d$pbar[o], pmin(1, d$pbar + hp)[o], lty = "22", col = pal$pw, lwd = 0.9)
    }
    graphics::abline(0, 1, col = pal$diag, lwd = 0.9)
    graphics::segments(d$pbar, d$pbar, d$pbar, ro, lwd = 3, lend = 1,
                       col = grDevices::adjustcolor(ifelse(ro > d$pbar, pal$above, pal$below), 0.5))
    dot(d$pbar, ro)
    graphics::box(col = pal$ink)
  }

  ## one key for the figure, in the outer margin, so it never covers a group
  graphics::par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  graphics::plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  if (band == "none") {
    lab <- c("EDGE direction", "Group"); fill <- c(NA, NA); lty <- c(1, NA)
    lc  <- c(pal$curve, "white"); pch <- c(NA, 21); bg <- c(NA, pal$pt)
  } else {
    pm  <- intToUtf8(177)                     # plus-minus sign; the source stays ASCII
    lab <- c(sprintf("%s %g%% band, %s%.2f%s", if (band == "simultaneous") "Simultaneous" else
                       "Pointwise", 100 * level, pm, crit, kind),
             if (band == "simultaneous") sprintf("Pointwise %g%% band, %s%.2f", 100 * level, pm, cp),
             "EDGE direction", "Group within the band", "Group beyond the band")
    nb   <- if (band == "simultaneous") 2L else 1L
    fill <- c(pal$band, rep(NA, length(lab) - 1L))
    lty  <- c(NA, if (nb == 2L) 2, 1, NA, NA)
    lc   <- c(NA, if (nb == 2L) pal$pw, pal$curve, "white", "white")
    pch  <- c(NA, if (nb == 2L) NA, NA, 21, 24)
    bg   <- c(NA, if (nb == 2L) NA, NA, pal$pt, pal$flag)
  }
  if (!any(scale == "residual")) {             # the curve is drawn on the residual scale only
    keep <- lab != "EDGE direction"
    lab <- lab[keep]; fill <- fill[keep]; lty <- lty[keep]; lc <- lc[keep]; pch <- pch[keep]; bg <- bg[keep]
  }
  graphics::legend("bottom", legend = lab, fill = fill, border = NA, lty = lty, col = lc, pch = pch,
                   pt.bg = bg, pt.cex = 1.1, lwd = ifelse(is.na(lty), NA, ifelse(lab == "EDGE direction", 2, 1)),
                   seg.len = 1.6, ncol = 2L, bty = "n", cex = 0.8, text.col = pal$ink,
                   x.intersp = 0.6, text.width = NULL, xpd = NA)
  invisible(x)
}

# Internal: the first of the axis titles that fits the width (side 1) or the height (side 2) of
# the current plot region, the last one whatever its length. Two panels side by side leave a
# long title too little room, and the plot would cut it at the figure edge.
.fit_label <- function(labs, side) {
  room <- graphics::par("pin")[if (side == 1) 1L else 2L]
  w <- graphics::strwidth(labs, units = "inches", cex = graphics::par("cex") * graphics::par("cex.lab"))
  labs[c(which(w <= room), length(labs))[1]]
}
