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
  expect_error(deepgof1(glm(y ~ 1, family = binomial())), "at least one covariate")
  expect_error(deepgof1(glm(y ~ x1, family = gaussian())), "family = binomial")
  x2 <- rnorm(n)
  f2 <- glm(y ~ x1 + x2, family = binomial())
  expect_error(deepgof1(f2, K = 8), "valid only at K")
  expect_error(deepgof1(f2, reading = "allpairs", covariates = c("x1", "x9")), "must name covariates")
  expect_error(deepgof1(f2, covariates = "x1"), "all-pairs and combined")
})

## ---- the all-pairs reading (2.9.0) -------------------------------------------------------

test_that("with two covariates the all-pairs reading is the default rule, bit for bit", {
  set.seed(31)
  n <- 100
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.4 * x1 - 0.6 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  set.seed(32); a <- deepgof1(fit, B = 29, reading = "allpairs")
  set.seed(32); b <- deepgof1(fit, B = 29)
  expect_identical(a$statistic, b$statistic)
  expect_identical(a$boot, b$boot)
  expect_identical(a$p.value, b$p.value)
})

test_that("the statistic is the largest pair score and the map belongs to that pair", {
  set.seed(33)
  n <- 120
  x1 <- runif(n, -3, 3); x2 <- rnorm(n); x3 <- rnorm(n); x4 <- runif(n)
  y <- rbinom(n, 1, plogis(0.3 * x1 + 0.3 * x2 + 0.8 * (x1^2 - 3)))
  fit <- glm(y ~ x1 + x2 + x3 + x4, family = binomial())
  r <- deepgof1(fit, B = 9, reading = "allpairs")
  expect_equal(nrow(r$pairs), choose(4, 2))
  expect_identical(r$statistic, max(r$pairs$score))
  top <- r$pairs[which.max(r$pairs$score), ]
  expect_identical(unname(r$axes), c(top$axis1, top$axis2))
  expect_equal(dim(r$map), c(6L, 6L))
  ## rescoring the returned map gives the statistic back
  M <- ebrahim.gof:::deepgof1_weights
  expect_equal(ebrahim.gof:::.dg_score((as.vector(t(r$map)) - M$mu) / M$sd), r$statistic)
  ## 'covariates' restricts the pairs
  rc <- deepgof1(fit, B = 9, reading = "allpairs", covariates = c("x2", "x3", "x4"))
  expect_equal(nrow(rc$pairs), 3L)
  expect_false("x1" %in% c(rc$pairs$axis1, rc$pairs$axis2))
})

test_that("with one covariate the map is 36 quantile cells along it", {
  set.seed(34)
  n <- 90
  x <- runif(n, -3, 3)
  y <- rbinom(n, 1, plogis(-1 + x))
  fit <- glm(y ~ x, family = binomial())
  cells <- ebrahim.gof:::.dg_cells(matrix(x, dimnames = list(NULL, "x")))
  ## the cells follow the ranks of x: the smallest values in cell 1, the largest in cell 36
  expect_identical(cells$idx[[1]][which.min(x)], 1L)
  expect_identical(cells$idx[[1]][which.max(x)], 36L)
  expect_true(all(tabulate(cells$idx[[1]], 36) %in% 2:3))
  set.seed(1); r <- deepgof1(fit, B = 19)
  set.seed(1); ra <- deepgof1(fit, B = 19, reading = "allpairs")
  expect_identical(r$axes, "x")
  expect_identical(r$p.value, ra$p.value)     # one covariate: both readings are the same map
  expect_equal(r$p.value * 20, round(r$p.value * 20))
})

test_that("the all-pairs p-value holds its level on a correct model", {
  skip_on_cran()
  ## a coarse check of calibration with three covariates, where the maximum over pairs is taken
  p <- vapply(1:40, function(s) {
    set.seed(400 + s)
    n <- 80
    x1 <- runif(n, -3, 3); x2 <- rnorm(n); x3 <- rt(n, 4)
    y <- rbinom(n, 1, plogis(-0.5 + 0.4 * x1 + 0.4 * x2 + 0.4 * x3))
    deepgof1(glm(y ~ x1 + x2 + x3, family = binomial()), B = 19, reading = "allpairs")$p.value
  }, 0)
  ## uniform on the 1/20 grid under the null: mean near .525, few small values
  expect_gt(mean(p), 0.35)
  expect_lt(sum(p <= 0.05), 8)
})

test_that("the maps are drawn over covariates, not over model-matrix columns", {
  set.seed(35)
  n <- 150
  d <- data.frame(x1 = runif(n, -3, 3), x2 = rnorm(n), g = factor(sample(c("a", "b", "c"), n, TRUE)))
  d$y <- rbinom(n, 1, plogis(0.5 * d$x1 - 0.4 * d$x2))
  ## a spline in x1 and an interaction: still two covariates, so one pair, ranked by x1 itself
  f <- glm(y ~ splines::ns(x1, 3) + x2 + x1:x2, data = d, family = binomial())
  r <- deepgof1(f, B = 9, reading = "allpairs")
  expect_identical(nrow(r$pairs), 1L)
  expect_identical(unname(r$axes), c("x1", "x2"))
  ## the same map as a model with x1 untransformed would be drawn on
  X <- ebrahim.gof:::.dg_covariates(f)
  expect_identical(colnames(X), c("x1", "x2"))
  expect_identical(unname(X[, "x1"]), d$x1)
  ## a three-level factor is one covariate, not two dummy columns
  fg <- glm(y ~ x1 + x2 + g, data = d, family = binomial())
  expect_identical(colnames(ebrahim.gof:::.dg_covariates(fg)), c("x1", "x2", "g"))
  ## rows dropped for missing values are dropped from the covariates too
  d2 <- d; d2$x2[c(3, 10)] <- NA
  f2 <- glm(y ~ I(x1^2) + x2, data = d2, family = binomial())
  X2 <- ebrahim.gof:::.dg_covariates(f2)
  expect_equal(nrow(X2), n - 2)
  expect_identical(unname(X2[, "x1"]), d$x1[-c(3, 10)])
})

test_that("a model fitted without a data argument still finds its covariates", {
  set.seed(36)
  n <- 80
  u <- runif(n); w <- rnorm(n)
  yy <- rbinom(n, 1, plogis(u - w))
  f <- glm(yy ~ log(u) + w, family = binomial())
  expect_identical(colnames(ebrahim.gof:::.dg_covariates(f)), c("u", "w"))
})

## Up to 2.8.0 the bootstrap refitted the formula on the model frame, which fails for every term
## that transforms a covariate: each replicate scored +Inf and the p-value was 1 whatever the data.
test_that("transformed terms refit in every bootstrap replicate", {
  set.seed(37)
  n <- 300
  d <- data.frame(x1 = runif(n, -3, 3), x2 = rnorm(n))
  d$y <- rbinom(n, 1, plogis(0.5 * d$x1 + 1.2 * (d$x1^2 - 3) - 0.5 * d$x2))
  f_log <- glm(y ~ log(x1 + 4) + x2, data = d, family = binomial())
  f_ns  <- glm(y ~ splines::ns(x2, 3) + x1, data = d, family = binomial())
  for (f in list(f_log, f_ns)) for (rd in c("axes", "allpairs", "combined", "columns")) {
    r <- deepgof1(f, B = 19, reading = rd)
    expect_true(all(is.finite(r$boot)))
    ## both models omit the quadratic in x1. The 2.8.0 column rule ("columns") can choose two spline
    ## columns of x2 and miss it; the axis rule over covariates reads x1 and x2 themselves
    if (rd != "columns") expect_lt(r$p.value, 0.1)
  }
})

test_that("the axis rule reads covariates, not spline or dummy columns", {
  set.seed(39)
  n <- 250
  d <- data.frame(x1 = runif(n, -3, 3), x2 = rnorm(n), g = factor(sample(letters[1:3], n, TRUE)))
  d$y <- rbinom(n, 1, plogis(0.4 * d$x1 - 0.8 * d$x2))
  f <- glm(y ~ splines::ns(x2, 3) + x1 + g, data = d, family = binomial())
  r <- deepgof1(f, B = 9)
  expect_true(all(r$axes %in% c("x1", "x2", "g")))
  ## on a model with untransformed numeric covariates the axis rule is the 2.8.0 rule, bit for bit
  f2 <- glm(y ~ x1 + x2, data = d, family = binomial())
  set.seed(40); a <- deepgof1(f2, B = 19)
  set.seed(40); b <- deepgof1(f2, B = 19, reading = "columns")
  expect_identical(a$boot, b$boot)
  expect_identical(a$p.value, b$p.value)
})

test_that("the combined reading returns a valid p-value and its two components", {
  set.seed(41)
  n <- 150
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.5 * x1 - 0.5 * x2))
  r <- deepgof1(glm(y ~ x1 + x2 + x3, family = binomial()), B = 19, reading = "combined")
  expect_true(r$p.value >= 1 / 20 && r$p.value <= 1)
  expect_equal(r$p.value * 20, round(r$p.value * 20))
  expect_named(r$components, c("axes", "allpairs"))
  ## the minimum of two p-values is never below the smaller of them after calibration
  expect_gte(r$p.value, min(r$components))
})

test_that("grouped or weighted binomial fits are refused", {
  set.seed(38)
  x1 <- rnorm(40); x2 <- rnorm(40); k <- rbinom(40, 5, 0.4)
  expect_error(deepgof1(glm(cbind(k, 5 - k) ~ x1 + x2, family = binomial())), "0/1 outcome")
  y <- rbinom(40, 1, 0.5)
  expect_error(deepgof1(glm(y ~ x1 + x2, family = binomial(), weights = rep(2, 40))), "0/1 outcome")
})
