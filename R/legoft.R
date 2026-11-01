## legoft.R -- the pretrained detector and its two-domain localization.
##
## Companion software for the LEGofT manuscript. Two user-facing functions:
##   legoft(fit)           -- one p-value: does the model fit?
##   legoft.localize(fit)  -- which of two evidence domains carries the misfit, with
##                            familywise error control by closed testing.
##
## Both calibrate by a parametric bootstrap at the fitted parameters. Nothing is
## retrained online: the combination weights were fixed offline once and ship frozen,
## so two analysts running this on the same data obtain the same p-value.

## The eleven battery members the frozen rule reads, and their two domains.
## INDEX  -- departures expressible as f(x'beta): the grouping tests read fitted risk,
##           which is a monotone transform of the index, and the directed tests read
##           bends in that same index. They are one axis, not two.
## COV    -- structure in directions orthogonal to beta.
.LEGOFT_MEMBERS <- c("HL", "HL-equalwidth", "Pigeon-Heyse", "EF",
                     "DEF.poly2", "DEF.poly3", "DEF.stukel",
                     "Tsiatis", "Xie", "A.poly", "A.spline")

.LEGOFT_DOMAINS <- list(
  INDEX = c("HL", "HL-equalwidth", "Pigeon-Heyse", "EF",
            "DEF.poly2", "DEF.poly3", "DEF.stukel"),
  COV   = c("Tsiatis", "Xie", "A.poly", "A.spline"))

## Frozen weights (gamma = 3, aggregation = mean), selected offline on simulated
## misfit families and never refitted. See the manuscript's Section on training.
.LEGOFT_WEIGHTS <- c(
  "HL"            = 0.0098, "HL-equalwidth" = 0.0258, "Pigeon-Heyse" = 0.0099,
  "EF"            = 0.0104, "DEF.poly2"     = 0.0409, "DEF.poly3"    = 0.0401,
  "DEF.stukel"    = 0.0346, "Tsiatis"       = 0.1776, "Xie"          = 0.1994,
  "A.poly"        = 0.2255, "A.spline"      = 0.2260)

.clip_p  <- function(p) pmin(pmax(p, 1e-12), 1 - 1e-12)
.cauchy  <- function(p) tan((0.5 - .clip_p(p)) * pi)

## weighted Cauchy score over a named p-value vector
.legoft_score <- function(pv, w) sum(w[names(w)] * .cauchy(pv[names(w)]), na.rm = TRUE)

## Monte Carlo p-value: rank of the observed statistic among reference draws.
## The (1 + #{ref >= obs}) / (1 + #ref) form is the exact-rank convention of Barnard
## (1963); it is valid for any fixed statistic and never returns 0.
.mc_p <- function(obs, ref) (1 + sum(ref >= obs, na.rm = TRUE)) / (1 + sum(is.finite(ref)))

## Battery p-values for one fit, restricted to the members the rule reads.
## Nine members come from the registry; the two covariate-space directed tests
## (A.poly, A.spline) are NOT registry entries -- they are cdef.gof() with the
## polynomial and spline bases -- so they are fetched separately and named to match.
## Stukel, Osius-Rojek and le-Cessie are not read by the frozen rule, but they belong to
## the panels of the companion panel paper, so they can be named in `weights` or `panel`.
## le-Cessie is a slow registry test (it builds an n x n kernel), so it is switched on
## only when a panel asks for it.
.LEGOFT_REGISTRY_MEMBERS <- c("HL", "HL-equalwidth", "Pigeon-Heyse", "EF",
                              "DEF.poly2", "DEF.poly3", "DEF.stukel", "Tsiatis", "Xie",
                              "Stukel", "Osius-Rojek", "le-Cessie")

## The panels studied in the panel paper (Ebrahim, "From a panel of goodness-of-fit
## tests to one decision"). Each is combined with equal weights.
.LEGOFT_PANELS <- list(
  equal         = .LEGOFT_MEMBERS,
  full          = c("HL", "HL-equalwidth", "Pigeon-Heyse", "Tsiatis", "Xie", "EF",
                    "DEF.poly2", "DEF.poly3", "DEF.stukel", "Stukel", "Osius-Rojek",
                    "le-Cessie", "A.poly", "A.spline"),
  classical     = c("HL", "Osius-Rojek", "le-Cessie", "Stukel", "Tsiatis", "Xie"),
  one.per.group = c("HL", "Xie", "DEF.poly3", "Osius-Rojek", "le-Cessie"),
  protected     = c("HL", "EF", "Pigeon-Heyse", "le-Cessie"))

## Combination statistic, larger = more misfit. "cauchy" is the weighted Cauchy sum used
## by the frozen rule; the other three use the member p-values with equal weight.
.panel_score <- function(pv, w, rule) {
  p <- .clip_p(pv[names(w)])
  switch(rule,
         cauchy = sum(w * .cauchy(p), na.rm = TRUE),
         hmp    = mean(1 / p, na.rm = TRUE),
         fisher = -2 * sum(log(p), na.rm = TRUE),
         minp   = -min(p, na.rm = TRUE))
}

.legoft_battery <- function(object, members = .LEGOFT_MEMBERS) {
  reg <- intersect(members, .LEGOFT_REGISTRY_MEMBERS)
  pv  <- stats::setNames(rep(NA_real_, length(members)), members)
  if (length(reg)) {
    res <- run.all.gof(object, tests = reg, include_slow = "le-Cessie" %in% reg,
                       install = "no")
    got <- stats::setNames(as.numeric(res$p_value), as.character(res$Test))
    pv[names(got)[names(got) %in% members]] <- got[names(got) %in% members]
  }
  if ("A.poly"   %in% members)
    pv[["A.poly"]]   <- tryCatch(cdef.gof(object, basis = "poly")$p_value,
                                 error = function(e) NA_real_)
  if ("A.spline" %in% members)
    pv[["A.spline"]] <- tryCatch(cdef.gof(object, basis = "spline")$p_value,
                                 error = function(e) NA_real_)
  pv[members]
}

## parametric bootstrap reference: refit under the fitted model, recompute the battery
.legoft_reference <- function(object, B, members, seed) {
  if (!is.null(seed)) {             # seed the draws here, then give the caller's stream back
    .seed_restore <- .rng_snapshot()
    on.exit(.rng_restore(.seed_restore), add = TRUE)
    set.seed(seed)
  }
  ph  <- as.numeric(object$fitted.values)
  dat <- object$data
  if (is.null(dat)) dat <- stats::model.frame(object)
  frm <- stats::formula(object)
  yv  <- all.vars(frm)[1]
  out <- matrix(NA_real_, B, length(members), dimnames = list(NULL, members))
  for (b in seq_len(B)) {
    db <- dat
    db[[yv]] <- stats::rbinom(length(ph), 1L, ph)
    fb <- tryCatch(suppressWarnings(stats::glm(frm, data = db, family = stats::binomial())),
                   error = function(e) NULL)
    if (!is.null(fb)) out[b, ] <- tryCatch(.legoft_battery(fb, members),
                                           error = function(e) rep(NA_real_, length(members)))
  }
  out[stats::complete.cases(out), , drop = FALSE]
}

#' LEGofT: frozen-weight combination goodness-of-fit test for binary logistic regression
#'
#' Combines eleven classical and directed goodness-of-fit statistics with weights that
#' were fixed offline and ship frozen, and calibrates the combination by a parametric
#' bootstrap at the fitted parameters. Nothing is retrained when you call it.
#'
#' @param object a fitted binomial \code{glm}.
#' @param B number of parametric-bootstrap replicates for the reference distribution.
#'   199 gives a smallest attainable p-value of 0.005; raise it for smaller p-values.
#' @param seed optional integer for reproducibility.
#' @param weights optional named vector of member weights; defaults to the frozen rule.
#'   Supplying your own makes the result no longer the shipped procedure -- say so if
#'   you report it. Members may be any of the eleven frozen-rule tests plus
#'   \code{"Stukel"}, \code{"Osius-Rojek"} and \code{"le-Cessie"}.
#' @param panel which tests to combine. \code{"legoft"} (default) is the frozen rule.
#'   The others combine a panel of the panel paper with equal weights:
#'   \code{"equal"} (the eleven frozen-rule tests), \code{"full"} (all fourteen tests),
#'   \code{"classical"} (Hosmer-Lemeshow, Osius-Rojek, le Cessie-van Houwelingen, Stukel,
#'   Tsiatis, Xie), \code{"one.per.group"} (one test from each of the five groups of
#'   correlated tests) and \code{"protected"} (Hosmer-Lemeshow, EF, Pigeon-Heyse,
#'   le Cessie-van Houwelingen; robust to a few recording errors in a covariate).
#'   Ignored when \code{weights} is given.
#' @param rule combination statistic: \code{"cauchy"} (default), \code{"hmp"} (harmonic
#'   mean p-value), \code{"fisher"} or \code{"minp"} (minimum p-value). The last three use
#'   equal weights, so they need a panel other than \code{"legoft"} or equal
#'   \code{weights}. Every rule is calibrated by the same parametric bootstrap.
#' @return an object of class \code{"legoft"}: the statistic, its bootstrap p-value,
#'   the member p-values, the weights, the panel, the rule and \code{B_used}.
#' @details
#' The p-value is exact in finite samples when the null parameters are known. With the
#' parameters estimated -- the case here -- the calibration is asymptotic; simulation at
#' \eqn{n = 500} put the empirical size at 0.037 against a nominal 0.05.
#'
#' The panel paper recommends \code{panel = "equal"} or \code{"full"} with the default
#' Cauchy rule for a single decision from many tests, and \code{panel = "protected"}
#' when a few covariate values may be recorded wrongly. In its simulations the
#' bootstrap-calibrated combination kept the 5\% level in three null designs. Panels
#' that contain le Cessie-van Houwelingen are slower, because that test builds an
#' \eqn{n \times n} matrix for the data and for every bootstrap replicate.
#' @seealso \code{\link{legoft.localize}} for which domain of evidence carries the misfit.
#' @examples
#' \donttest{
#' set.seed(1)
#' x1 <- runif(300, -3, 3); x2 <- rnorm(300)
#' y  <- rbinom(300, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.25 * x1^2))
#' fit <- glm(y ~ x1 + x2, family = binomial())
#' legoft(fit, B = 99, seed = 1)                       # frozen LEGofT weights
#' legoft(fit, B = 99, seed = 1, panel = "equal")      # equal weights, same tests
#' legoft(fit, B = 99, seed = 1, panel = "protected")  # robust to recording errors
#' }
#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept LEGofT
#' @concept frozen weights
#' @concept familywise error
#' @export
legoft <- function(object, B = 199, seed = NULL, weights = NULL,
                   panel = c("legoft", "equal", "full", "classical", "one.per.group", "protected"),
                   rule = c("cauchy", "hmp", "fisher", "minp")) {
  if (!inherits(object, "glm") || !identical(stats::family(object)$family, "binomial"))
    stop("legoft() needs a fitted binomial glm")
  panel <- match.arg(panel)
  rule  <- match.arg(rule)
  if (!is.null(weights)) {
    unknown <- setdiff(names(weights), c(.LEGOFT_MEMBERS, .LEGOFT_REGISTRY_MEMBERS))
    if (is.null(names(weights)) || length(unknown))
      stop("weights must be named by battery members; unknown: ",
           paste(unknown, collapse = ", "))
    w <- weights; panel <- "custom"
  } else if (panel == "legoft") {
    w <- .LEGOFT_WEIGHTS
  } else {
    m <- .LEGOFT_PANELS[[panel]]
    w <- stats::setNames(rep(1 / length(m), length(m)), m)
  }
  if (rule != "cauchy" && length(unique(round(w, 12))) > 1)
    stop("rule = \"", rule, "\" combines with equal weights; choose a panel other than ",
         "\"legoft\" or give equal weights")
  members <- names(w)
  pv  <- .legoft_battery(object, members)
  if (all(is.na(pv))) stop("no battery member returned a p-value")
  REF <- .legoft_reference(object, B, members, seed)
  if (nrow(REF) < 20) stop("too few usable bootstrap replicates (", nrow(REF), ")")
  obs <- .panel_score(pv, w, rule)
  ref <- apply(REF, 1, function(r) .panel_score(r, w, rule))
  method <- if (panel == "legoft" && rule == "cauchy")
    "LEGofT (frozen weights, parametric-bootstrap calibration)" else
    sprintf("Panel combination (panel = %s, rule = %s, %d tests, parametric-bootstrap calibration)",
            panel, rule, length(w))
  structure(list(statistic = obs, p.value = .mc_p(obs, ref), members = pv,
                 B_used = nrow(REF), weights = w, panel = panel, rule = rule,
                 method = method),
            class = "legoft")
}

#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept LEGofT
#' @concept frozen weights
#' @concept familywise error
#' @export
print.legoft <- function(x, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat(sprintf("  statistic = %.4f    p-value = %.4f    (B = %d)\n\n",
              x$statistic, x$p.value, x$B_used))
  invisible(x)
}

#' Localize misspecification with familywise error control
#'
#' Reports which of two domains of evidence carries the misfit, using closed testing
#' over the domains. A domain is implicated only when every intersection containing it
#' rejects, so the probability of implicating any collection of domains whose pooled
#' evidence is jointly exchangeable with the reference draws' is at most \code{alpha}.
#'
#' @param object a fitted binomial \code{glm}.
#' @param B,seed as in \code{\link{legoft}}.
#' @param alpha familywise level.
#' @return an object of class \code{"legoft_localize"}: the verdict, the intersection
#'   p-values, and the member p-values.
#' @details
#' The two domains are defined by what the members can see, not by taxonomy.
#' \strong{INDEX} holds the seven statistics that read the linear index -- the grouping
#' tests, which stratify on fitted risk, and the directed tests, which examine bends in
#' that same index. \strong{COV} holds the four that read directions orthogonal to it.
#' The calibration and link readings are pooled deliberately: they are not separately
#' identifiable, because fitted risk is a monotone transform of the index.
#'
#' A verdict says where the \emph{evidence} lies. It does not say which part of the
#' model to repair: a departure of one kind can move members assigned to the other
#' domain, and an unimplicated domain is not thereby certified correct.
#' @seealso \code{\link{legoft}}
#' @examples
#' \donttest{
#' set.seed(1)
#' x1 <- runif(400, -3, 3); x2 <- rnorm(400)
#' y  <- rbinom(400, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.4 * x1 * x2))
#' fit <- glm(y ~ x1 + x2, family = binomial())
#' legoft.localize(fit, B = 99, seed = 1)
#' }
#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept LEGofT
#' @concept frozen weights
#' @concept familywise error
#' @export
legoft.localize <- function(object, B = 199, seed = NULL, alpha = 0.05) {
  if (!inherits(object, "glm") || !identical(stats::family(object)$family, "binomial"))
    stop("legoft.localize() needs a fitted binomial glm")
  members <- .LEGOFT_MEMBERS
  pv  <- .legoft_battery(object, members)
  REF <- .legoft_reference(object, B, members, seed)
  if (nrow(REF) < 20) stop("too few usable bootstrap replicates (", nrow(REF), ")")

  ## Three intersections. The TOP one uses the LEGofT statistic itself, which the
  ## closure proof permits (it asks only that each intersection statistic be a fixed
  ## measurable function of the pooled member p-values). Taking it makes the two
  ## procedures coherent by construction: no domain can be implicated unless legoft()
  ## rejects. The price is a familywise rate slightly below nominal.
  pool_c <- function(m) {
    o <- mean(.cauchy(pv[m]), na.rm = TRUE)
    r <- apply(REF[, m, drop = FALSE], 1, function(z) mean(.cauchy(z), na.rm = TRUE))
    .mc_p(o, r)
  }
  pI <- pool_c(.LEGOFT_DOMAINS$INDEX)
  pC <- pool_c(.LEGOFT_DOMAINS$COV)
  oT <- .legoft_score(pv, .LEGOFT_WEIGHTS)
  rT <- apply(REF, 1, function(r) .legoft_score(r, .LEGOFT_WEIGHTS))
  pT <- .mc_p(oT, rT)

  flag <- c(INDEX = unname(pI <= alpha && pT <= alpha),
            COV   = unname(pC <= alpha && pT <= alpha))
  structure(list(verdict = if (any(flag)) names(flag)[flag] else character(0),
                 intersection = c(INDEX = pI, COV = pC, BOTH = pT),
                 members = pv, alpha = alpha, B_used = nrow(REF),
                 domains = .LEGOFT_DOMAINS),
            class = "legoft_localize")
}

#' @concept goodness-of-fit
#' @concept calibration
#' @concept logistic regression
#' @concept model diagnostics
#' @concept LEGofT
#' @concept frozen weights
#' @concept familywise error
#' @export
print.legoft_localize <- function(x, ...) {
  cat("\nLEGofT-Localize: closed testing over two domains of evidence\n\n")
  cat(sprintf("  intersection p-values:  INDEX = %.4f   COV = %.4f   BOTH = %.4f\n",
              x$intersection[["INDEX"]], x$intersection[["COV"]], x$intersection[["BOTH"]]))
  v <- if (length(x$verdict)) paste(x$verdict, collapse = " + ") else "no domain implicated"
  cat(sprintf("  verdict at FWER %.2f:    %s\n\n", x$alpha, v))
  if (length(x$verdict))
    cat("  A verdict locates the evidence. It does not name the repair.\n\n")
  invisible(x)
}
