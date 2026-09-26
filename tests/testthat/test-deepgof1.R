## deepgof1(): the shipped weights must reproduce the training framework, and the test
## must behave (valid p-value, correct grid, sensible calls rejected).

test_that("the base-R forward pass reproduces the training framework", {
  ## deepgof1_verify holds standardized maps and the scores PyTorch produced for them.
  got <- apply(deepgof1_verify$z, 1, function(z) ebrahim.gof:::.dg_score(z))
  expect_equal(got, deepgof1_verify$score, tolerance = 1e-6)
})

test_that("deepgof1() returns a valid p-value on a correctly specified model", {
  skip_on_cran()
  set.seed(11)
  n <- 120
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  r <- deepgof1(fit, B = 49)
  expect_s3_class(r, "deepgof1")
  expect_true(r$p.value >= 1 / 50 && r$p.value <= 1)
  expect_equal(r$p.value * 50, round(r$p.value * 50))   # lies on the 1/(B+1) grid
  expect_identical(sort(r$axes), c("x1", "x2"))
})

test_that("deepgof1() detects a strong omitted quadratic", {
  skip_on_cran()
  set.seed(12)
  n <- 200
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2 + 1.2 * (x1^2 - mean(x1^2))))
  fit <- glm(y ~ x1 + x2, family = binomial())
  expect_lt(deepgof1(fit, B = 99)$p.value, 0.05)
})

## Tied covariates (binary, integer) used to be put into map cells by row order, so a data
## file sorted by the outcome leaked y into the map. Ties are now broken at random, once per
## call, and the same random order serves the observed and the bootstrap maps.
tied_sorted_data <- function(seed, n = 200) {
  set.seed(seed)
  x1 <- sample(1:4, n, replace = TRUE); x2 <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, plogis(-1 + 0.4 * x1 + 0.8 * x2))
  d <- data.frame(y, x1, x2)
  d[order(d$y), ]                                   # all events last, as in many clinical files
}

test_that("the map does not depend on row order once the tie-break travels with the rows", {
  d <- tied_sorted_data(21)
  n <- nrow(d)
  set.seed(22); tb <- sample.int(n); perm <- sample.int(n)
  f1 <- glm(y ~ x1 + x2, data = d, family = binomial())
  f2 <- glm(y ~ x1 + x2, data = d[perm, ], family = binomial())
  expect_equal(ebrahim.gof:::.dg_map(f1, 6L, tb)$map,
               ebrahim.gof:::.dg_map(f2, 6L, tb[perm])$map)
})

test_that("deepgof1() does not reject a correct model because the rows are sorted by y", {
  skip_on_cran()
  d <- tied_sorted_data(23)
  fit <- glm(y ~ x1 + x2, data = d, family = binomial())
  p <- vapply(1:8, function(s) { set.seed(100 + s); deepgof1(fit, B = 19)$p.value }, 0)
  ## with row-order ties this model was rejected at the smallest attainable p on every seed
  expect_gt(mean(p), 0.2)
  expect_lt(sum(p <= 0.05), 5)
  ## the shuffled file gives p-values from the same law: compare the observed statistics,
  ## whose spread over seeds comes only from the random tie-break
  ds <- d[sample.int(nrow(d)), ]
  fs <- glm(y ~ x1 + x2, data = ds, family = binomial())
  s_sorted <- vapply(1:8, function(s) { set.seed(200 + s); deepgof1(fit, B = 1)$statistic }, 0)
  s_shuf   <- vapply(1:8, function(s) { set.seed(200 + s); deepgof1(fs, B = 1)$statistic }, 0)
  expect_gt(suppressWarnings(stats::wilcox.test(s_sorted, s_shuf))$p.value, 0.01)
})

test_that("without ties deepgof1() draws nothing extra, so earlier results reproduce", {
  set.seed(24)
  n <- 80
  x1 <- rnorm(n); x2 <- runif(n)
  y <- rbinom(n, 1, plogis(0.5 * x1 - x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  ph <- fitted(fit)
  set.seed(25); invisible(deepgof1(fit, B = 3)); a <- runif(1)
  set.seed(25); for (b in 1:3) rbinom(n, 1L, ph); z <- runif(1)
  expect_identical(a, z)
})

test_that("deepgof1() rejects calls it cannot serve", {
  set.seed(13)
  n <- 60
  x1 <- runif(n); y <- rbinom(n, 1, 0.5)
  expect_error(deepgof1(glm(y ~ x1, family = binomial())), "at least two covariates")
  expect_error(deepgof1(glm(y ~ x1, family = gaussian())), "family = binomial")
  x2 <- rnorm(n)
  f2 <- glm(y ~ x1 + x2, family = binomial())
  expect_error(deepgof1(f2, K = 8), "valid only at K")
})
