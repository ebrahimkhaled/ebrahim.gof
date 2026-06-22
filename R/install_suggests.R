# Optional ("Suggests") packages used by the run.all.gof() battery, and a
# user-consented installer for them. CRAN forbids a package from installing
# anything silently, so nothing here ever installs without either an interactive
# yes/no confirmation (the default) or the caller explicitly setting ask = FALSE,
# and run.all.gof() only ever reaches the installer in an interactive session.

# Which optional package each test needs to run (only tests that are skipped
# entirely when the package is absent; Stukel still runs without statmod).
.GOF_TEST_PKGS <- list(
  "GiViTI"          = c("givitiR", "callr"),
  "GiViTI-external" = c("givitiR", "callr"),
  "HL-GAM"          = "mgcv",
  "PR-GAM"          = "mgcv",
  "Xie-GAM"         = "mgcv",
  "BAGofT"          = "BAGofT",
  "Lai-Liu-HL"      = "ResourceSelection"
)

# The full optional set, including the ones that only improve (not gate) a test.
.GOF_ALL_SUGGESTS <- c("givitiR", "callr", "mgcv", "BAGofT",
                       "ResourceSelection", "statmod", "CompQuadForm")

#' Install the optional packages used by \code{run.all.gof()}
#'
#' The slow tests in \code{\link{run.all.gof}} rely on optional packages that
#' live in \code{Suggests} (\pkg{givitiR} and \pkg{callr} for the GiViTI
#' calibration test, \pkg{mgcv} for the GAM tests, \pkg{BAGofT} for the adaptive
#' test, and \pkg{ResourceSelection} for the Lai-Liu test). Per CRAN policy the
#' package never installs them on its own; this helper installs the missing ones
#' for you, asking first.
#'
#' @param pkgs Optional character vector of package names to consider. Defaults
#'   to the full optional set used by the battery.
#' @param ask Logical; when \code{TRUE} (the default in interactive sessions)
#'   you are shown the list of missing packages and asked to confirm before
#'   anything is installed. Set \code{ask = FALSE} to install without a prompt
#'   (e.g. in a setup script you control).
#' @return Invisibly, the character vector of packages that were missing.
#' @examples
#' \dontrun{
#' # install whatever optional packages are missing, after confirming:
#' gof_install_suggests()
#'
#' # just the GiViTI dependencies, no prompt:
#' gof_install_suggests(c("givitiR", "callr"), ask = FALSE)
#' }
#' @seealso \code{\link{run.all.gof}}
#' @export
gof_install_suggests <- function(pkgs = NULL, ask = interactive()) {
  want    <- if (is.null(pkgs)) .GOF_ALL_SUGGESTS else pkgs
  missing <- want[!vapply(want, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) {
    message("All optional packages for run.all.gof() are already installed.")
    return(invisible(character(0)))
  }
  if (isTRUE(ask)) {
    if (!interactive()) {
      message("Optional package(s) not installed: ", paste(missing, collapse = ", "),
              ". Run gof_install_suggests() interactively, or with ask = FALSE, to install.")
      return(invisible(missing))
    }
    ans <- readline(sprintf(
      "Install %d optional package(s) for run.all.gof() now: %s? [y/N] ",
      length(missing), paste(missing, collapse = ", ")))
    if (!tolower(trimws(ans)) %in% c("y", "yes")) {
      message("Skipped. You can install them later with gof_install_suggests().")
      return(invisible(missing))
    }
  }
  utils::install.packages(missing)
  still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still))
    warning("Could not install: ", paste(still, collapse = ", "))
  invisible(missing)
}

# Called by run.all.gof(): when interactive and tests in the run need optional
# packages that aren't installed, offer to install them (consent-gated). Never
# does anything in a non-interactive session, so R CMD check / scripts / CRAN are
# untouched and the library is never modified without the user's say-so.
.gof_maybe_install <- function(sel, include_slow, install = "ask") {
  install <- match.arg(install, c("ask", "no", "yes"))
  if (identical(install, "no") || !interactive()) return(invisible())
  eff <- sel
  if (!isTRUE(include_slow))
    eff <- eff[!vapply(eff, function(nm) isTRUE(.GOF_REGISTRY[[nm]]$slow), logical(1))]
  needed  <- unique(unlist(.GOF_TEST_PKGS[intersect(eff, names(.GOF_TEST_PKGS))],
                           use.names = FALSE))
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible())
  gof_install_suggests(missing, ask = identical(install, "ask"))
  invisible()
}
