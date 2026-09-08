# Localize misspecification with familywise error control

Reports which of two domains of evidence carries the misfit, using
closed testing over the domains. A domain is implicated only when every
intersection containing it rejects, so the probability of implicating
any collection of domains whose pooled evidence is jointly exchangeable
with the reference draws' is at most `alpha`.

## Usage

``` r
legoft.localize(object, B = 199, seed = NULL, alpha = 0.05)
```

## Arguments

- object:

  a fitted binomial `glm`.

- B, seed:

  as in
  [`legoft`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/legoft.md).

- alpha:

  familywise level.

## Value

an object of class `"legoft_localize"`: the verdict, the intersection
p-values, and the member p-values.

## Details

The two domains are defined by what the members can see, not by
taxonomy. **INDEX** holds the seven statistics that read the linear
index – the grouping tests, which stratify on fitted risk, and the
directed tests, which examine bends in that same index. **COV** holds
the four that read directions orthogonal to it. The calibration and link
readings are pooled deliberately: they are not separately identifiable,
because fitted risk is a monotone transform of the index.

A verdict says where the *evidence* lies. It does not say which part of
the model to repair: a departure of one kind can move members assigned
to the other domain, and an unimplicated domain is not thereby certified
correct.

## See also

[`legoft`](https://ebrahimkhaled.github.io/ebrahim.gof/reference/legoft.md)

## Examples

``` r
# \donttest{
set.seed(1)
x1 <- runif(400, -3, 3); x2 <- rnorm(400)
y  <- rbinom(400, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 0.4 * x1 * x2))
fit <- glm(y ~ x1 + x2, family = binomial())
legoft.localize(fit, B = 99, seed = 1)
#> 
#> LEGofT-Localize: closed testing over two domains of evidence
#> 
#>   intersection p-values:  INDEX = 0.0100   COV = 0.0200   BOTH = 0.0200
#>   verdict at FWER 0.05:    INDEX + COV
#> 
#>   A verdict locates the evidence. It does not name the repair.
#> 
# }
```
