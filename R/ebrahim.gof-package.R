#' ebrahim.gof: Goodness-of-Fit and Calibration Tests for Logistic Regression
#'
#' A unified toolbox of goodness-of-fit and calibration tests for binary
#' logistic regression, callable in a single line via \code{\link{run.all.gof}}.
#' The package is aimed particularly at \emph{sparse} data, where the classical
#' Hosmer--Lemeshow test loses power.
#'
#' @section The author's own tests:
#' \itemize{
#'   \item \code{\link{ef.gof}} --- the omnibus Ebrahim--Farrington (EF) test for
#'     binary data with automatic grouping.
#'   \item \code{\link{def.gof}} / \code{\link{edge.gof}} --- the Directed EF
#'     (\dQuote{EDGE}) test, which spends its few degrees of freedom on the
#'     smooth calibration-shape directions where structured misfit concentrates.
#'   \item \code{\link{def.ensemble.gof}} --- a Cauchy-combination ensemble of
#'     the directed bases.
#' }
#'
#' @section Aggregated tests (for comparison):
#' \code{\link{run.all.gof}} also runs, in one call, a wide range of classical
#' and modern tests --- Hosmer--Lemeshow, McCullagh, Osius--Rojek,
#' le Cessie--van Houwelingen, Stute--Zhu, the binary-adaptive \pkg{BAGofT} test,
#' and the \pkg{givitiR} calibration test. Each aggregated test is obtained from
#' its own package (where installed) and is attributed to its authors; these are
#' provided for head-to-head comparison, not claimed as original to this package.
#'
#' @section Data:
#' \code{\link{gof_demo}} is a bundled example dataset with a documented,
#' reproducible misfit for illustrating the battery.
#'
#' @seealso The vignette
#'   \code{vignette("ebrahim-gof-toolbox", package = "ebrahim.gof")}.
#'
#' @docType package
#' @name ebrahim.gof-package
#' @aliases ebrahim.gof
#' @keywords internal
"_PACKAGE"
