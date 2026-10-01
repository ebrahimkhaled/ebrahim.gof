# EDGE for Streaming Data: Calibration Monitoring Without Recomputation

Keeps the directed test of a frozen model up to date as new patients
arrive, without revisiting the records already seen. The external-mode
statistic of
[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md)
depends on the data only through four sums per risk group – observed
events, expected events, the binomial variance and the count – so each
new record updates one group in constant time, and the test is
recomputed from the \\G\\ group summaries alone.

## Usage

``` r
edge.stream(p_ref = NULL, breaks = NULL, G = 10, basis = "poly3")

# S3 method for class 'edge_stream'
update(object, y, p, ...)

# S3 method for class 'edge_stream'
summary(object, ...)

# S3 method for class 'edge_stream'
print(x, ...)
```

## Arguments

- p_ref:

  Optional numeric vector of reference predicted probabilities; the cut
  points are its `G`-quantiles.

- breaks:

  Optional increasing numeric vector of interior cut points in (0, 1)
  (`G - 1` of them); used instead of `p_ref`.

- G:

  Number of risk groups (default `10`).

- basis:

  Calibration basis: `"poly3"` (default), `"poly2"`, `"stukel"` or
  `"sym"`, as in
  [`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md).

- object:

  An `edge_stream` object.

- y:

  Binary (0/1) outcomes of the new records.

- p:

  Predicted probabilities of the new records, made without their
  outcomes.

- ...:

  Unused.

- x:

  An `edge_stream` object.

## Value

An object of class `edge_stream`. Add data with `update(object, y, p)`;
read the test with `summary(object)`, a one-row `data.frame` in the
format of `def.gof(..., external = TRUE)` with the number of records and
the smallest group added.

## Details

**The partition.** The groups are fixed in advance by cut points on the
predicted risk: either given as `breaks`, or the \\G\\-quantiles of a
reference set of predictions `p_ref` (for example the development data,
or the first batch). Because the cut points depend on predictions only
and never on outcomes, the statistic keeps its \\\chi^2\_{d+1}\\
reference under a calibrated model however the risk distribution of
later patients drifts: the groups need not stay of equal size. Small
groups weaken the protection against corrupted records that
equal-frequency groups give, so
[`summary()`](https://rdrr.io/r/base/summary.html) reports the smallest
group.

**Exactness.** After any sequence of updates the statistic equals the
one-shot external statistic computed on all records seen, with the same
partition: the sums are additive, so the order and the batching of the
updates do not matter.

**Repeated looks.** One test at any time is valid. Testing after every
batch and acting on the first rejection is a sequential procedure, and
the plain chi-squared reference does not control its overall false-alarm
rate; spend the level across looks (for example, Bonferroni over a
planned number of looks) or test at fixed times.

## See also

[`def.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/def.gof.md),
[`run.all.external`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.external.md).

## Author

Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>

## Examples

``` r
set.seed(1)
p_dev <- plogis(rnorm(5000, -1.5, 1))          # predictions on the development data
s <- edge.stream(p_ref = p_dev, G = 10)
for (day in 1:30) {                            # a deployed model, 100 patients a day
  p <- plogis(rnorm(100, -1.5, 1))
  y <- rbinom(100, 1, p)
  s <- update(s, y, p)
}
summary(s)
#>   Test Basis Test_Statistic df Method   p_value    n smallest_group
#> 1 EDGE poly3       5.732058  4 stream 0.2200719 3000            262
```
