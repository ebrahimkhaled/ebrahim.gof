# Startup hint: when ebrahim.gof is attached, check whether the optional
# ("Suggests") packages that unlock the slow tests in run.all.gof() are
# installed, and, if any are missing, print a single suppressible hint naming
# them and the one command that installs them. Nothing is ever installed here;
# this is only a message. Uses packageStartupMessage so it is silenced by
# suppressPackageStartupMessages() and never prints during R CMD check output
# parsing. .GOF_ALL_SUGGESTS and .GOF_TEST_PKGS are defined in install_suggests.R.

.onAttach <- function(libname, pkgname) {
  opt <- get0(".GOF_ALL_SUGGESTS", envir = asNamespace(pkgname),
              ifnotfound = character(0))
  if (!length(opt)) return(invisible())
  missing <- opt[!vapply(opt, requireNamespace, logical(1), quietly = TRUE)]
  if (!length(missing)) return(invisible())

  packageStartupMessage(
    sprintf("ebrahim.gof: %d optional package%s not installed: %s.",
            length(missing), if (length(missing) == 1L) "" else "s",
            paste(missing, collapse = ", ")),
    "\n  The core battery still runs; these only add the extra slow tests",
    "\n  (GAM, BAGofT, GiViTI, Lai-Liu-HL). Install the missing ones with:",
    "\n      gof_install_suggests()",
    "\n  Or let run.all.gof(..., include_slow = TRUE) offer to install them.")
  invisible()
}
