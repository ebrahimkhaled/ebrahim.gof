# Plot the residual map of a DeepGOF-1 test

Draws the `K` by `K` residual map that the network scored in
[`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md)
or
[`deepgof1.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.external.md),
with the value of every cell written in it. A cell holds the observed
minus the expected number of events among the records whose two
covariates fall in that pair of rank groups, divided by its null
standard deviation, so under a correct model the cells are roughly
standard normal. Vermilion cells had more events than the model
predicted, blue cells fewer.

## Usage

``` r
# S3 method for class 'deepgof1'
plot(x, colour = TRUE, digits = 1, ...)
```

## Arguments

- x:

  A `deepgof1` object.

- colour:

  `FALSE` shades the cells in greys by the size of the value only; the
  sign is then read from the number in the cell.

- digits:

  Decimals shown in the cells. Default `1`.

- ...:

  Not used.

## Value

`x`, invisibly.

## Details

The axes are the rank groups of the two covariates the map was drawn
over (the `axes` element of the result), 1 being the lowest. With one
covariate the map is `K * K` quantile cells along its ranks, read row by
row from the top left. The p-value is that of the test; the map says
where the misfit sits, and a single large cell is not a test.

## See also

[`deepgof1`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.md),
[`deepgof1.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/deepgof1.external.md).

## Examples

``` r
# \donttest{
set.seed(1)
n  <- 400
x1 <- rnorm(n); x2 <- rnorm(n)
y  <- rbinom(n, 1, plogis(0.5 * x1 + 0.5 * x2 + 0.8 * x1 * x2))
fit <- glm(y ~ x1 + x2, family = binomial())
plot(deepgof1(fit, B = 49))

# }
```
