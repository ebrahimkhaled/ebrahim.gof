## Plot method for localization verdicts (localize.external(), localize.gof()). Two panels:
##   compass  one wedge per part; its length is the strength of evidence, -log10 p of the part's own test (capped at
##            p = 0.001); the dashed ring is alpha; the parts named are drawn in full colour
##   lattice  every set of the parts tested, with the p-value of its test; tinted when it rejects; a part is named when
##            every set that contains it rejects, and the links from a named part up to those sets are drawn dark
## Parts that cannot be tested in the setting (intercept and slope in-sample) are shown as empty wedges.

.lplot_pal <- function(colour) {
  if (colour) list(part = c(INTERCEPT = "#2E6FAC", SLOPE = "#2E9A6E", LINK = "#8E4A9E", COV = "#D2622A"),
                   tint = c(INTERCEPT = "#DCE8F4", SLOPE = "#DCF0E6", LINK = "#ECDFF0", COV = "#F8E3D8"),
                   rej = "#EEF1F5", ink = "#2B2B2B", soft = "#8A8F98", rule = "#C9CDD3")
  else list(part = stats::setNames(rep("#262626", 4), .LOC_ORDER), tint = stats::setNames(rep("#E3E3E3", 4), .LOC_ORDER),
            rej = "#EDEDED", ink = "#1F1F1F", soft = "#7A7A7A", rule = "#C4C4C4")
}
.lplot_lp <- function(p) pmin(-log10(pmax(p, 1e-3)), 3) / 3
.lplot_pfmt <- function(p) ifelse(p < .001, "p < .001", paste0("p = ", sub("^0", "", sprintf("%.3f", p))))
.lplot_split <- function(s) strsplit(s, "+", fixed = TRUE)[[1]]
.lplot_rrect <- function(x0, y0, x1, y1, r = .12, ...) {     # rectangle with rounded corners
  r <- min(r, (x1 - x0) / 2, (y1 - y0) / 2); t <- seq(0, pi / 2, length.out = 8)
  graphics::polygon(c(x1 - r + r * cos(t), x0 + r - r * sin(t), x0 + r - r * cos(t), x1 - r + r * sin(t)),
                    c(y1 - r + r * sin(t), y1 - r + r * cos(t), y0 + r - r * sin(t), y0 + r - r * cos(t)), ...)
}

.lplot_compass <- function(x, colour = TRUE, head = NULL) {
  pal <- .lplot_pal(colour); G <- .LOC_ORDER; alpha <- if (is.null(x$alpha)) 0.05 else x$alpha
  mid <- c(90, 0, 270, 180) * pi / 180                          # intercept top, slope right, link bottom, cov left
  graphics::plot(NA, xlim = c(-1.6, 1.6), ylim = c(-1.5, 1.85), asp = 1, axes = FALSE, xlab = "", ylab = "")
  if (!is.null(head)) graphics::text(-1.6, 1.85, head, adj = c(0, 1), font = 2, cex = .82, col = pal$ink)
  tt <- seq(0, 2 * pi, length.out = 300)
  for (q in c(.01, .001)) graphics::lines(.lplot_lp(q) * cos(tt), .lplot_lp(q) * sin(tt), col = pal$rule, lwd = .6)
  graphics::lines(.lplot_lp(alpha) * cos(tt), .lplot_lp(alpha) * sin(tt), col = pal$ink, lwd = .9, lty = "22")
  tested <- G %in% names(x$intersection); nm <- G %in% x$named
  for (k in seq_along(G)) {
    a <- seq(mid[k] - pi / 4 + .13, mid[k] + pi / 4 - .13, length.out = 40)
    if (!tested[k]) { graphics::polygon(c(0, 1.0 * cos(a)), c(0, 1.0 * sin(a)), col = NA, border = pal$rule, lty = 3); next }
    rr <- max(.lplot_lp(x$intersection[[G[k]]]), .025)
    graphics::polygon(c(0, rr * cos(a)), c(0, rr * sin(a)), col = if (nm[k]) pal$part[[G[k]]] else pal$tint[[G[k]]],
                      border = if (nm[k]) pal$part[[G[k]]] else grDevices::adjustcolor(pal$part[[G[k]]], .55), lwd = .8)
  }
  for (q in c(alpha, .01, .001))
    graphics::text(.lplot_lp(q) * cos(pi / 4) + .02, .lplot_lp(q) * sin(pi / 4) + .02, sub("^0", "", format(q)), adj = c(0, 0),
                   cex = .5, col = pal$soft)
  lx <- 1.33 * cos(mid); ly <- 1.33 * sin(mid) + c(.05, 0, -.05, 0)
  graphics::text(lx, ly + .07, tolower(G), cex = .8, font = ifelse(nm, 2, 1),
                 col = ifelse(tested, ifelse(nm, pal$part[G], pal$ink), pal$soft))
  graphics::text(lx, ly - .1, ifelse(tested, .lplot_pfmt(as.numeric(x$intersection[G])), "not tested"), cex = .58, col = pal$soft)
  invisible(x)
}

.lplot_lattice <- function(x, colour = TRUE, head = NULL) {
  pal <- .lplot_pal(colour); alpha <- if (is.null(x$alpha)) 0.05 else x$alpha
  sets <- names(x$intersection); lev <- vapply(sets, function(s) length(.lplot_split(s)), 1L); K <- max(lev)
  hh <- function(k) .12 * k + .17; hw <- .64
  Y <- .55; for (k in seq_len(K - 1)) Y[k + 1] <- Y[k] + hh(k) + hh(k + 1) + .62
  xy <- do.call(rbind, lapply(seq_len(K), function(k) { s <- sets[lev == k]; nk <- length(s); W <- min(1.55, 8.4 / max(nk, 1))
    data.frame(s = s, k = k, x = 5.2 + (seq_len(nk) - (nk + 1) / 2) * W, y = Y[k], stringsAsFactors = FALSE) }))
  top <- max(xy$y + hh(xy$k)); H <- 6.3; sh <- max(0, (H - top - .5) / 2); xy$y <- xy$y + sh; top <- top + sh
  graphics::plot(NA, xlim = c(0, 10), ylim = c(0, max(H, top + .5)), axes = FALSE, xlab = "", ylab = "")
  if (!is.null(head)) graphics::text(0, max(H, top + .5) - .03, head, adj = c(0, 1), font = 2, cex = .82, col = pal$ink)
  for (k in seq_len(K)) graphics::text(0, Y[k] + sh, if (k == 1) "each part" else if (k == K) (if (K == 2) "both" else paste("all", K)) else paste(k, "parts"),
                                       adj = 0, cex = .55, col = pal$soft, font = 3)
  sup <- function(a, b) all(.lplot_split(a) %in% .lplot_split(b))
  hot <- vapply(seq_len(nrow(xy)), function(j) any(vapply(x$named, function(g) sup(g, xy$s[j]), TRUE)), TRUE)
  for (i in seq_len(nrow(xy))) for (j in seq_len(nrow(xy)))
    if (xy$k[j] == xy$k[i] + 1 && sup(xy$s[i], xy$s[j])) {
      on <- hot[i] && hot[j] && any(vapply(x$named, function(g) sup(g, xy$s[i]), TRUE))
      graphics::segments(xy$x[i], xy$y[i] + hh(xy$k[i]), xy$x[j], xy$y[j] - hh(xy$k[j]),
                         col = if (on) pal$ink else pal$rule, lwd = if (on) 1.1 else .6)
    }
  for (i in seq_len(nrow(xy))) {
    s <- xy$s[i]; k <- xy$k[i]; p <- x$intersection[[s]]; rej <- p <= alpha; named <- k == 1 && s %in% x$named
    .lplot_rrect(xy$x[i] - hw, xy$y[i] - hh(k), xy$x[i] + hw, xy$y[i] + hh(k),
               col = if (named) pal$part[[s]] else if (rej) pal$rej else "white",
               border = if (named) pal$part[[s]] else if (rej) pal$soft else pal$rule, lwd = if (named) 1.5 else .8)
    nm <- .lplot_split(s)
    graphics::text(xy$x[i], xy$y[i] + hh(k) - .16 - (seq_along(nm) - 1) * .24, tolower(nm), cex = .56,
                   col = if (named) "white" else pal$part[nm], font = if (named) 2 else 1)
    graphics::text(xy$x[i], xy$y[i] - hh(k) + .13, .lplot_pfmt(p), cex = .47,
                   col = if (named) "white" else if (rej) pal$ink else pal$soft, font = if (rej) 2 else 1)
  }
  invisible(x)
}

#' Plot a localization verdict
#'
#' Draws the verdict of \code{\link{localize.external}} or \code{\link{localize.gof}} in two
#' panels. The \emph{misfit compass} has one wedge per part (intercept, slope, link, cov); the
#' length of a wedge is the strength of evidence, \eqn{-\log_{10} p} of the part's own test
#' (capped at \eqn{p = 0.001}), the dashed ring marks \code{alpha}, and the parts named are filled.
#' The \emph{lattice of closed tests} shows every set of the parts tested with the p-value of its
#' test, shaded when it rejects; a part is named when every set that contains it rejects, and the
#' links from a named part up to those sets are drawn dark. The highest part named and the next
#' step are written under the panels. Parts that cannot be tested in the setting (intercept and
#' slope in-sample) are shown as empty wedges marked "not tested".
#'
#' @param x A \code{gof_localize} object.
#' @param which \code{"both"} (default), \code{"compass"} or \code{"lattice"}.
#' @param colour \code{FALSE} draws in greys only, for journals that ask for black on white.
#' @param main Title; by default the setting and the sample size.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @examples
#' set.seed(1)
#' n <- 500
#' X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)
#' y <- rbinom(n, 1, plogis(qlogis(p) + 0.8 * X$x1 * X$x2))
#' r <- localize.external(y, p, X, M = 199, seed = 1)
#' plot(r)
#' plot(r, which = "compass", colour = FALSE)
#' @importFrom graphics plot
#' @exportS3Method plot gof_localize
plot.gof_localize <- function(x, which = c("both", "compass", "lattice"), colour = TRUE, main = NULL, ...) {
  which <- match.arg(which)
  op <- graphics::par(mar = c(.3, .3, .3, .3), oma = c(1.6, 0, 1.6, 0)); on.exit(graphics::par(op))
  if (which == "both") graphics::layout(matrix(1:2, 1), widths = c(.42, .58)) else graphics::layout(1)
  if (which != "lattice") .lplot_compass(x, colour, "Misfit compass")
  if (which != "compass") .lplot_lattice(x, colour, "Closed tests")
  if (is.null(main)) main <- paste0(toupper(substr(x$setting, 1, 1)), substring(x$setting, 2),
                                    if (is.null(x$n)) "" else paste0(", n = ", format(x$n, big.mark = ",")))
  graphics::mtext(main, side = 3, outer = TRUE, line = .3, adj = .01, cex = .85, font = 2)
  graphics::mtext(if (is.na(x$highest)) "Nothing named: keep the model; the verdict is limited by the sample size."
                  else paste0("Highest part named: ", tolower(x$highest), ". Next step: ", x$action, "."),
                  side = 1, outer = TRUE, line = .4, adj = .01, cex = .75)
  invisible(x)
}
