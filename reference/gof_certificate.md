# Can DeepGOF-1 detect this departure? A certificate from the frozen weights

Decides, without running the test, whether
[`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md)
is consistent against a departure you name: its power tends to one
precisely when the certificate below is positive. The answer costs one
forward pass on the shipped weights – no bootstrap, no simulation, no
data beyond the fitted model.

## Usage

``` r
gof_certificate(
  fit,
  pi_star = NULL,
  eta_shift = NULL,
  direction = NULL,
  K = 6L
)
```

## Arguments

- fit:

  a fitted `glm` with `family = binomial()`: the null model whose
  adequacy is in question. Not needed if `direction` is given.

- pi_star:

  true success probabilities under the alternative, one per observation.

- eta_shift:

  the omitted term on the linear-predictor scale, one per observation;
  used only when `pi_star` is not given.

- direction:

  a `K` by `K` matrix, or a vector of length `K * K` in row-major order,
  giving the map direction directly.

- K:

  grid resolution. Leave at 6: the shipped weights are valid at no other
  resolution.

## Value

An object of class `"gof_certificate"`: `certificate` (the recession
slope at unit map norm; positive certifies consistency), `verdict`,
`mirror` (the slope of the same departure with its sign flipped, since a
test may see a departure one way and not the other), `direction`
(normalized), `axes` and `method`.

## Details

Give the departure in whichever form you have it. `pi_star` is the
vector of true success probabilities under the alternative, one per
observation of `fit`. `eta_shift` is the term the fitted linear
predictor is missing, so that the truth is `plogis(eta_hat + eta_shift)`
– an omitted quadratic is `eta_shift = c * x^2`. `direction` skips the
data entirely and takes a `K` by `K` map direction, which is what a
theoretical question ("a checkerboard at the cell scale") usually looks
like.

A non-positive certificate does not mean the test never rejects; it
means the guarantee of power tending to one is unavailable for that
direction, which for this network happens on a thin cone of
sub-resolution oscillation – the failure mode every binned statistic
shares. Because the grid reads two covariates, a departure that
integrates to zero within every cell has direction zero and is certified
blind, which is the representation limit of the method stated as a
computation rather than a caveat.

With a fitted model the direction is built on the map of the default
reading of
[`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md)
(`reading = "axes"`): the two covariates of its axis rule, or the 36
cells along the only covariate. For the all-pairs reading, name the
pair's map yourself through `direction`.

## References

Ebrahim EK, Hussein OAE-A, El-Kotory A (2026). "Where Does a Logistic
Risk Model Fail? An Audited Neural Goodness-of-Fit Test for Model
Development and External Validation." arXiv:2609.29575 \[stat.ME\].
[doi:10.48550/arXiv.2609.29575](https://doi.org/10.48550/arXiv.2609.29575)

## See also

[`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md)

## Examples

``` r
set.seed(1)
n  <- 200
x1 <- runif(n, -3, 3); x2 <- rnorm(n)
y  <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2))
fit <- glm(y ~ x1 + x2, family = binomial())

# an omitted quadratic in x1: is this test consistent against it?
gof_certificate(fit, eta_shift = 0.9 * (x1^2 - mean(x1^2)))
#> 
#>  DeepGOF-1 consistency certificate (recession slope of the frozen network)
#> 
#> certificate = 1.2534 (mirror 1.7450)
#> grid: 6 x 6 over ranks of x1 and x2
#> detectable: power tends to one
#> 

# a map direction named directly: a threshold departure along the first axis
gof_certificate(direction = matrix(rep(c(-1, -1, -1, 1, 1, 1), each = 6), 6, 6, byrow = TRUE))
#> 
#>  DeepGOF-1 consistency certificate (recession slope of the frozen network)
#> 
#> certificate = 2.1999 (mirror 1.0747)
#> detectable: power tends to one
#> 

# the blind cone is thin but not empty: multi-start minimization over the sphere reaches
# -0.281 for these weights, at a direction that oscillates cell by cell. Directions like
# that are the ones this function is for.
```
