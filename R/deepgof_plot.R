#' Plot the residual map of a DeepGOF-1 test
#'
#' Draws the \code{K} by \code{K} residual map that the network scored in
#' \code{\link{deepgof1}} or \code{\link{deepgof1.external}}, with the value of every cell
#' written in it. A cell holds the observed minus the expected number of events among the
#' records whose two covariates fall in that pair of rank groups, divided by its null standard
#' deviation, so under a correct model the cells are roughly standard normal. Vermilion cells
#' had more events than the model predicted, blue cells fewer.
#'
#' The axes are the rank groups of the two covariates the map was drawn over (the
#' \code{axes} element of the result), 1 being the lowest. With one covariate the map is
#' \code{K * K} quantile cells along its ranks, read row by row from the top left. The
#' p-value is that of the test; the map says where the misfit sits, and a single large cell is
#' not a test.
#'
#' @param x A \code{deepgof1} object.
#' @param colour \code{FALSE} shades the cells in greys by the size of the value only; the sign
#'   is then read from the number in the cell.
#' @param digits Decimals shown in the cells. Default \code{1}.
#' @param ... Not used.
#' @return \code{x}, invisibly.
#' @examples
#' \donttest{
#' set.seed(1)
#' n  <- 400
#' x1 <- rnorm(n); x2 <- rnorm(n)
#' y  <- rbinom(n, 1, plogis(0.5 * x1 + 0.5 * x2 + 0.8 * x1 * x2))
#' fit <- glm(y ~ x1 + x2, family = binomial())
#' plot(deepgof1(fit, B = 49))
#' }
#' @seealso \code{\link{deepgof1}}, \code{\link{deepgof1.external}}.
#' @importFrom graphics axis box image layout lcm par text
#' @importFrom grDevices colorRampPalette
#' @exportS3Method plot deepgof1
plot.deepgof1 <- function(x, colour = TRUE, digits = 1, ...) {
  m <- x$map
  if (is.null(m) || !is.matrix(m)) stop("this deepgof1 object holds no residual map", call. = FALSE)
  K   <- nrow(m)
  lim <- max(2, max(abs(m), na.rm = TRUE))
  nc  <- 101L
  pal <- if (colour) grDevices::colorRampPalette(c("#0072B2", "#FFFFFF", "#D55E00"))(nc)
         else grDevices::colorRampPalette(c("#4D4D4D", "#FFFFFF", "#4D4D4D"))(nc)
  one <- length(x$axes) == 1L

  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  graphics::layout(matrix(1:2, 1), widths = c(1, graphics::lcm(2.1)))
  graphics::par(cex = 1)                     # layout() shrinks the text in a multi-panel figure
  graphics::par(mar = c(3.6, 3.6, 0.6, 0.4), mgp = c(2.2, 0.5, 0), tcl = -0.25, las = 1,
                cex.axis = 0.85, cex.lab = 0.95)

  ## image() puts z[i, j] at (i, j); map[i, j] is rank group i of the first axis and j of the
  ## second. With one covariate, cell (i - 1) K + j is drawn in row i from the top, column j.
  z <- if (one) t(m[K:1, , drop = FALSE]) else m
  graphics::image(seq_len(K), seq_len(K), z, zlim = c(-lim, lim), col = pal, axes = FALSE,
                  xlab = "", ylab = "")
  if (one) {
    graphics::title(xlab = .fit_label(c(sprintf("Quantile cells of %s, low to high, row by row", x$axes),
                                        sprintf("Cells of %s, row by row", x$axes),
                                        "Cells, row by row"), 1))
  } else {
    graphics::title(xlab = .fit_label(c(sprintf("%s (rank group)", x$axes[1]), x$axes[1]), 1),
                    ylab = .fit_label(c(sprintf("%s (rank group)", x$axes[2]), x$axes[2]), 2))
  }
  if (one) {
    graphics::axis(1, at = seq_len(K), labels = seq_len(K))
    graphics::axis(2, at = seq_len(K), labels = (K:1 - 1L) * K + 1L)
  } else {
    graphics::axis(1, at = seq_len(K))
    graphics::axis(2, at = seq_len(K))
  }
  gx <- rep(seq_len(K), K); gy <- rep(seq_len(K), each = K)
  v  <- z[cbind(gx, gy)]
  vr <- round(v, digits)
  vr[vr == 0] <- 0                           # no "-0.0"
  graphics::text(gx, gy, formatC(vr, format = "f", digits = digits), cex = 0.8,
                 col = ifelse(abs(v) > 0.6 * lim, "white", "#262626"))
  graphics::box(col = "#404040")

  ## the key: the same scale as the cells
  graphics::par(mar = c(3.6, 0.3, 0.6, 2.3))
  s <- seq(-lim, lim, length.out = nc)
  graphics::image(1, s, matrix(s, 1), col = pal, axes = FALSE, xlab = "", ylab = "")
  graphics::axis(4, at = pretty(c(-lim, lim)))
  graphics::box(col = "#404040")
  invisible(x)
}
