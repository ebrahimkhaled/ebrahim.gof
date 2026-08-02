#' Synthetic binary outcome data with a smooth calibration misfit
#'
#' A small, fully synthetic dataset for demonstrating the goodness-of-fit and
#' calibration battery. It was generated reproducibly (see
#' \code{data-raw/make_gof_demo.R}) from a logistic data-generating process whose
#' true linear predictor includes a \emph{quadratic} term in (standardized) age.
#' A model that regresses \code{outcome} on \code{age} linearly (together with
#' \code{bmi}, \code{sex} and \code{treatment}) is therefore mildly
#' misspecified, through a smooth, low-dimensional calibration distortion. This
#' is the regime in which the directed Ebrahim--Farrington / EDGE test
#' (\code{\link{edge.gof}}, \code{\link{def.gof}}) is designed to have more power
#' than classical omnibus tests such as Hosmer--Lemeshow.
#'
#' @format A data frame with 800 rows and 5 variables:
#' \describe{
#'   \item{outcome}{binary response, 0/1 (event rate about 0.27).}
#'   \item{age}{continuous covariate, years (range about 20--70). The true model
#'     depends on age quadratically.}
#'   \item{bmi}{continuous covariate, body mass index in kg/m^2.}
#'   \item{sex}{binary covariate, 0 = female, 1 = male.}
#'   \item{treatment}{binary covariate, 0 = control, 1 = treated.}
#' }
#'
#' @details The true data-generating linear predictor is
#'   \deqn{\eta = -0.6 + 0.8 z_a - 0.7 z_a^2 + 0.5 z_b + 0.4\,\mathrm{sex}
#'         - 0.3\,\mathrm{treatment},}
#'   where \eqn{z_a = (\mathrm{age} - 45)/14} and
#'   \eqn{z_b = (\mathrm{bmi} - 27)/4}, and
#'   \eqn{\Pr(\mathrm{outcome} = 1) = \mathrm{plogis}(\eta)}.
#'
#' @examples
#' data("gof_demo", package = "ebrahim.gof")
#' fit <- glm(outcome ~ age + bmi + sex + treatment,
#'            data = gof_demo, family = binomial)
#' edge.gof(fit)
#'
#' @source Simulated; see \code{data-raw/make_gof_demo.R} in the package sources.
#' @keywords datasets
"gof_demo"

#' Grouped-covariate companion to \code{gof_demo} (replicated covariate patterns)
#'
#' A companion dataset to \code{\link{gof_demo}} built from the \emph{same}
#' data-generating process and seed discipline (see
#' \code{data-raw/make_gof_demo.R}), except that the covariates are coarsened
#' \emph{before} the linear predictor is computed: \code{age} is rounded to
#' 10-year bins (20, 30, ..., 70) and \code{bmi} to whole integers. The recorded
#' covariates are therefore exactly the covariates the outcome was generated
#' from, and many observations share a covariate pattern (328 distinct patterns
#' among 800 observations, versus one pattern per observation in
#' \code{gof_demo}).
#'
#' Its purpose is to demonstrate the sparse-versus-grouped distinction that
#' \code{\link{run.all.gof}} surfaces: the battery reports the per-observation
#' ("sparse") and per-covariate-pattern ("grouped") forms of the
#' pattern-sensitive tests (Pearson, deviance, McCullagh) side by side, and on
#' replicated-pattern data such as this the two forms can disagree on the same
#' fitted model. On the all-continuous \code{gof_demo} (every observation its
#' own pattern) the two forms coincide -- the degenerate case.
#'
#' @format A data frame with 800 rows and 5 variables:
#' \describe{
#'   \item{outcome}{binary response, 0/1 (event rate about 0.27).}
#'   \item{age}{age in years, rounded to 10-year bins (20, 30, ..., 70). The
#'     true model depends on (binned) age quadratically.}
#'   \item{bmi}{body mass index in kg/m^2, rounded to whole integers.}
#'   \item{sex}{binary covariate, 0 = female, 1 = male.}
#'   \item{treatment}{binary covariate, 0 = control, 1 = treated.}
#' }
#'
#' @details The true data-generating linear predictor has the same form as for
#'   \code{\link{gof_demo}},
#'   \deqn{\eta = -0.6 + 0.8 z_a - 0.7 z_a^2 + 0.5 z_b + 0.4\,\mathrm{sex}
#'         - 0.3\,\mathrm{treatment},}
#'   with \eqn{z_a = (\mathrm{age} - 45)/14} and
#'   \eqn{z_b = (\mathrm{bmi} - 27)/4} computed from the \emph{binned}
#'   \code{age} and \code{bmi}, and
#'   \eqn{\Pr(\mathrm{outcome} = 1) = \mathrm{plogis}(\eta)}.
#'
#' @examples
#' data("gof_demo_grouped", package = "ebrahim.gof")
#' fit <- glm(outcome ~ age + bmi + sex + treatment,
#'            data = gof_demo_grouped, family = binomial)
#' # sparse and grouped forms reported side by side:
#' run.all.gof(fit, include_slow = FALSE, install = "no")
#'
#' @source Simulated; see \code{data-raw/make_gof_demo.R} in the package sources.
#' @seealso \code{\link{gof_demo}}, \code{\link{run.all.gof}}
#' @keywords datasets
"gof_demo_grouped"
