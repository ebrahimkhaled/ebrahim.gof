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
