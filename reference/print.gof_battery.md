# Print a goodness-of-fit battery

Formats the
[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md)
result as a compact, readable table: rows grouped by test family,
p-values shown to four decimals (or scientific for very small values,
`"-"` when not available), and a significance flag. The object is still
a plain `data.frame` underneath, so all the raw columns remain available
for programmatic use.

## Usage

``` r
# S3 method for class 'gof_battery'
print(x, ...)
```

## Arguments

- x:

  A `gof_battery` object returned by
  [`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md).

- ...:

  Ignored.

## Value

`x`, invisibly.
