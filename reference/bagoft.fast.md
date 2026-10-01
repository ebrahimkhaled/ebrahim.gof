# Fast BAGofT: the Binary Adaptive Goodness-of-Fit Test

`bagoft.fast()` computes the BAGofT test of Zhang, Ding and Yang (2023)
for a fitted binary `glm`, and returns the same p-values as the BAGofT
package (version 1.0.0) for the same random seed – identical, not
approximately equal. It is the call

    BAGofT::BAGofT(testGlmBi(formula, link), parRF(), data,
                   nsplits = 100, ne = floor(5 * sqrt(n)), nsim = 100)

with the package's own defaults, computed without its per-split overhead
(formula parsing, model frames, `xtabs`, `cut` labels and the `predict`
wrappers, repeated \\(nsim + 1) \times nsplits\\ times).

## Usage

``` r
bagoft.fast(
  object,
  data = NULL,
  nsplits = 100,
  nsim = 100,
  ne = NULL,
  ntree = 60,
  Kmax = NULL,
  nmin = NULL,
  mtry = NULL,
  maxnodes = NULL
)
```

## Arguments

- object:

  A fitted binary [`glm`](https://rdrr.io/r/stats/glm.html)
  (`family = binomial`, any link) with unit prior weights and no offset.

- data:

  Optional data frame to run the test on, as in `BAGofT(data = ...)`. It
  must contain the response (numeric 0/1) and the model's variables;
  every other column is offered to the forest. Default: the model frame
  of `object`.

- nsplits:

  Number of random splits (default 100, the package default).

- nsim:

  Number of simulated responses used to calibrate the p-value (default
  100, the package default). `nsim = 0` returns the statistics only.

- ne:

  Size of the held-out part of each split. Default `floor(5 * sqrt(n))`,
  the package default.

- ntree, Kmax, nmin, mtry, maxnodes:

  The random-forest partitioner's settings, as in
  [`BAGofT::parRF()`](https://rdrr.io/pkg/BAGofT/man/parRF.html): number
  of trees (default 60), maximum number of groups (default
  `floor(ne / nmin)`), minimum group size (default `ceiling(sqrt(ne))`),
  variables tried at each node (default 1, the value the package ends up
  with), and maximum number of terminal nodes (default
  `min(n - ne, 5 * ncol)`, `ncol` the number of covariate columns).

## Value

An object of class `"bagoft_fast"`: a list with `p.value`, `p.value2`,
`p.value3` (the calibrated p-values for the mean, median and minimum
split p-value), `pmean`, `pmedian`, `pmin` (the observed statistics),
`simRes` (the simulated statistics), and `settings`. The elements shared
with [`BAGofT::BAGofT()`](https://rdrr.io/pkg/BAGofT/man/BAGofT.html)
have the same names and values. `singleSplit.results` is not returned.

## Details

**The test.** Each split fits the model on \\n - n_e\\ observations,
grows a random forest of the Pearson residuals on the covariates, uses
the forest to choose an adaptive partition of the covariate space, and
computes a Hosmer-Lemeshow-type chi-squared statistic on the \\n_e\\
held-out observations. The split p-values are averaged over `nsplits`
splits, and the average is calibrated against `nsim` responses simulated
from the full-data fit. The calibrated p-value is `p.value`; `p.value2`
and `p.value3` calibrate the median and the minimum instead.

**Read `p.value`, not `pmean` or `pmin`.** `pmean`, `pmedian` and `pmin`
are the observed *statistics* (summaries of the split p-values), not
p-values, and are not uniform under the null.

**What the forest sees.** As in the package's default
`parRF(parVar = ".")`, the forest partitions on every column of `data`
except the response – not only the model's terms. By default `data` is
the model frame, so the forest sees the model's variables; pass a wider
`data` to let it look for structure in covariates the model leaves out,
exactly as the package does. With more than five such columns the
package's pre-selection runs first in each split (`dcPre()`: the five
columns with the largest distance correlation with the residuals,
computed by dcov), and so it does here.

**Differences from the package.** None in the result where the package
runs. Where it does not, this function still does: with a single
covariate BAGofT 1.0.0 stops (its `parRF()` drops the data frame to a
vector), while here the forest simply uses the one column; and a formula
with transformed terms (e.g. `log(x)`) is evaluated through the fitted
model's own design matrix. Offsets and non-unit prior weights are not
supported, as in the package.

**Speed.** The work that remains is the forests themselves: \\(nsim + 1)
\times nsplits = 10{,}100\\ forests with the defaults. At \\n = 200\\
with two covariates the default test took 99 s against 239 s for the
package, and 188 s at \\n = 500\\ (package about 360 s); the forests
themselves are about half of it. Reduce `nsim` or `nsplits` for a
quicker, noisier answer.

## References

Zhang J, Ding J, Yang Y (2023). "Is a classification procedure good
enough? A goodness-of-fit assessment tool for classification learning."
*Journal of the American Statistical Association*, 118(542), 1115–1125.
[doi:10.1080/01621459.2021.1979010](https://doi.org/10.1080/01621459.2021.1979010)

Zhang J, Ding J, Yang Y (2021). BAGofT: A Binary Regression Adaptive
Goodness-of-Fit Test. R package version 1.0.0.
<https://CRAN.R-project.org/package=BAGofT>

## See also

[`run.all.gof`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/run.all.gof.md),
where this is the engine of the row `"BAGofT"`.

## Author

The procedure and the code it is adapted from are by Jiawei Zhang, Jie
Ding and Yuhong Yang (package BAGofT, GPL-3). The fast re-expression is
by Ebrahim Khaled Ebrahim <ebrahimkhaled@alexu.edu.eg>.

## Examples

``` r
# \donttest{
if (requireNamespace("randomForest", quietly = TRUE)) {
  set.seed(1)
  n  <- 100
  x1 <- rnorm(n); x2 <- rnorm(n)
  y  <- rbinom(n, 1, plogis(0.5 * x1 + 0.5 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  set.seed(2)
  r <- bagoft.fast(fit, nsplits = 20, nsim = 20)    # reduced for speed
  r
  ## the same numbers from the BAGofT package, same seed:
  ## set.seed(2)
  ## BAGofT::BAGofT(BAGofT::testGlmBi(y ~ x1 + x2, link = "logit"),
  ##                BAGofT::parRF(), data = data.frame(y, x1, x2),
  ##                nsplits = 20, nsim = 20)$p.value
}
#> 
#>  BAGofT (Zhang, Ding and Yang 2023), fast implementation
#> 
#> model: y ~ x1 + x2
#> nsplits = 20, nsim = 20, ne = 50, ntree = 60
#> p-value = 0.45   (mean split p-value; median: 0.3, minimum: 0.75)
#> statistics: pmean = 0.3929, pmedian = 0.3126, pmin = 0.005661 (not p-values)
#> 
# }
```
