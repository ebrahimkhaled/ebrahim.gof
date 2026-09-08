# Install the optional packages used by `run.all.gof()`

The slow tests in
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
rely on optional packages that live in `Suggests` (givitiR and callr for
the GiViTI calibration test, mgcv for the GAM tests, BAGofT for the
adaptive test, and ResourceSelection for the Lai-Liu test). Per CRAN
policy the package never installs them on its own; this helper installs
the missing ones for you, asking first.

## Usage

``` r
gof_install_suggests(pkgs = NULL, ask = interactive(), update = FALSE)
```

## Arguments

- pkgs:

  Optional character vector of package names to consider. Defaults to
  the full optional set used by the battery.

- ask:

  Logical; when `TRUE` (the default in interactive sessions) you are
  shown the list of packages to be installed/updated and asked to
  confirm first. Set `ask = FALSE` to proceed without a prompt (e.g. in
  a setup script you control).

- update:

  Logical; when `FALSE` (the default) only the *missing* packages are
  installed and anything already present is left untouched. When `TRUE`
  the function also checks (via
  [`old.packages`](https://rdrr.io/r/utils/update.packages.html)) which
  of the present packages are out of date and offers to update those
  too. The update check contacts your CRAN mirror, so it is a little
  slower.

## Value

Invisibly, the character vector of packages that were installed or
updated (empty if nothing was needed).

## See also

[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# install whatever optional packages are missing, after confirming:
gof_install_suggests()

# also update any that are out of date:
gof_install_suggests(update = TRUE)

# just the GiViTI dependencies, no prompt:
gof_install_suggests(c("givitiR", "callr"), ask = FALSE)
} # }
```
