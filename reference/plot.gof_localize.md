# Plot a localization verdict

Draws the verdict of
[`localize.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.external.md)
or
[`localize.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/localize.gof.md)
in two panels. The *misfit compass* has one wedge per part (intercept,
slope, link, cov); the length of a wedge is the strength of evidence,
\\-\log\_{10} p\\ of the part's own test (capped at \\p = 0.001\\), the
dashed ring marks `alpha`, and the parts named are filled. The *lattice
of closed tests* shows every set of the parts tested with the p-value of
its test, shaded when it rejects; a part is named when every set that
contains it rejects, and the links from a named part up to those sets
are drawn dark. The lattice writes the intercept as "intcpt", so the
label fits its node at 7 pt. The highest part named and the next step
are written under the panels. Parts that cannot be tested in the setting
(intercept and slope in-sample) are shown as empty wedges marked "not
tested".

## Usage

``` r
# S3 method for class 'gof_localize'
plot(
  x,
  which = c("both", "compass", "lattice"),
  colour = TRUE,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A `gof_localize` object.

- which:

  `"both"` (default), `"compass"` or `"lattice"`.

- colour:

  `FALSE` draws in greys only, for journals that ask for black on white.

- main:

  Title; by default the setting and the sample size.

- ...:

  Not used.

## Value

`x`, invisibly.

## Examples

``` r
set.seed(1)
n <- 500
X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
p <- plogis(-0.5 + 0.8 * X$x1 + 0.6 * X$x2)
y <- rbinom(n, 1, plogis(qlogis(p) + 0.8 * X$x1 * X$x2))
r <- localize.external(y, p, X, M = 199, seed = 1)
plot(r)

plot(r, which = "compass", colour = FALSE)
```
