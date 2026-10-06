#' Summarise a goodness-of-fit battery by family
#'
#' Condenses a \code{\link{run.all.gof}} or \code{\link{run.all.external}} result to one line
#' per test family: how many tests ran, how many returned a p-value, how many reject at
#' \code{alpha}, and the smallest p-value with the test that gave it. It answers "which kind
#' of departure do the data point to" before the full table is read.
#'
#' The counts are descriptive. The tests within a family are strongly correlated and no
#' multiplicity adjustment is made, so a count of rejections is not itself a test; for an
#' error-controlled statement of which part of the model is wrong, use
#' \code{\link{localize.gof}} or \code{\link{localize.external}}.
#'
#' The rows that \code{\link{run.all.gof}} prints for comparison only (the sparse
#' \code{Pearson} and \code{Deviance} statistics and the \code{F-test}, which do not hold
#' their level on ungrouped binary data) are counted in \code{Tests} and \code{With_p} but
#' not in \code{Rejected}, and never give \code{Smallest_p}.
#'
#' @param object A \code{gof_battery} object.
#' @param alpha Level at which a p-value counts as a rejection. Default \code{0.05}.
#' @param ... Not used.
#' @return An object of class \code{"summary.gof_battery"}: a data frame with columns
#'   \code{Family}, \code{Tests} (rows in the battery), \code{With_p} (rows with a p-value),
#'   \code{Rejected} (p-values at most \code{alpha}, comparison rows left out),
#'   \code{Smallest_p} and \code{Smallest_test},
#'   in the family order of the printed battery, with the attribute \code{alpha}.
#' @examples
#' fit <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial)
#' summary(run.all.gof(fit, include_slow = FALSE, install = "no"))
#' @seealso \code{\link{run.all.gof}}, \code{\link{print.gof_battery}}.
#' @exportS3Method summary gof_battery
summary.gof_battery <- function(object, alpha = 0.05, ...) {
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha <= 0 || alpha >= 1)
    stop("'alpha' must be a single number in (0, 1).", call. = FALSE)
  d <- as.data.frame(object)
  fam <- unique(d$Family)
  ord <- match(fam, .battery_family_order); ord[is.na(ord)] <- 99L
  fam <- fam[order(ord, seq_along(fam))]
  out <- do.call(rbind, lapply(fam, function(f) {
    s <- d[d$Family == f, , drop = FALSE]
    p <- s$p_value
    has <- !is.na(p)
    use <- has & .gof_counted(s)
    j <- if (any(use)) which(use)[which.min(p[use])] else NA_integer_
    data.frame(Family = f, Tests = nrow(s), With_p = sum(has), Rejected = sum(use & p <= alpha),
               Smallest_p = if (is.na(j)) NA_real_ else p[j],
               Smallest_test = if (is.na(j)) NA_character_ else s$Test[j],
               stringsAsFactors = FALSE)
  }))
  if (is.null(out))                            # a battery with no rows
    out <- data.frame(Family = character(0), Tests = integer(0), With_p = integer(0),
                      Rejected = integer(0), Smallest_p = numeric(0), Smallest_test = character(0),
                      stringsAsFactors = FALSE)
  rownames(out) <- NULL
  attr(out, "alpha") <- alpha
  class(out) <- c("summary.gof_battery", "data.frame")
  out
}

#' @rdname summary.gof_battery
#' @param x A \code{summary.gof_battery} object.
#' @exportS3Method print summary.gof_battery
print.summary.gof_battery <- function(x, ...) {
  d <- as.data.frame(x)
  alpha <- attr(x, "alpha")
  fp <- vapply(d$Smallest_p, function(p)
    if (is.na(p)) "-" else if (p < 1e-4) sprintf("%.1e", p) else formatC(p, format = "f", digits = 4),
    character(1))
  w <- max(40L, getOption("width", 80L))
  by <- ifelse(is.na(d$Smallest_test), "", d$Smallest_test)
  ## the test that gives the smallest p goes last; a long name is cut so a row fits the console
  ## one leading space, then each column and its separating space; one column of slack
  room <- w - 2L - (max(nchar(d$Family), 6L) + 1L) - 6L - 7L - 7L - (max(nchar(fp), 5L) + 1L)
  room <- max(room, 12L)
  by <- ifelse(nchar(by) > room, paste0(substr(by, 1L, room - 3L), "..."), by)
  tab <- data.frame(Family = d$Family, Tests = d$Tests, `With p` = d$With_p,
                    Reject = d$Rejected, `Min p` = fp, `Given by` = by,
                    check.names = FALSE, stringsAsFactors = FALSE)
  cat("\n", paste(strwrap(sprintf(
    "Goodness-of-fit battery by family: %d %s, %d with a p-value, %d %s at %g",
    sum(d$Tests), ngettext(sum(d$Tests), "test", "tests"), sum(d$With_p), sum(d$Rejected),
    ngettext(sum(d$Rejected), "rejects", "reject"), alpha), width = w - 1L), collapse = "\n"),
    "\n\n", sep = "")
  print(tab, row.names = FALSE, right = FALSE)
  cat("\n", paste(strwrap(paste(
    "Counts are descriptive: the tests in a family are correlated and not adjusted for",
    "multiplicity. Rows printed for comparison only are not counted as rejections."),
    width = w - 1L), collapse = "\n"), "\n", sep = "")
  invisible(x)
}

# the family order of the printed battery, shared by print.gof_battery() and the summary
.battery_family_order <- c("Global", "Standardized", "Partition", "Covariate-space", "Directed",
                           "Smoothing", "GAM", "Bootstrap", "Calibration", "Ensemble")
