# Plot the GiViTI calibration belt from a goodness-of-fit battery

Draws the GiViTI calibration belt stored on a
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
result that was produced with `calibration_plot = TRUE`. The belt shows
the fitted calibration curve with a confidence region against the
45-degree line.

## Usage

``` r
# S3 method for class 'gof_battery'
plot(x, ...)
```

## Arguments

- x:

  A `gof_battery` object from
  [`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md).

- ...:

  Passed to the givitiR plot method.

## Value

`x`, invisibly.
