make_fit <- function(seed = 1, n = 600, link = "logit") {
  set.seed(seed)
  x <- runif(n, -3, 3)
  eta <- 0.6 * x
  p <- if (link == "cloglog") 1 - exp(-exp(eta)) else 1 / (1 + exp(-eta))
  glm(rbinom(n, 1, p) ~ x, family = binomial())
}

test_that("def.gof returns the documented shape", {
  res <- def.gof(make_fit())
  expect_s3_class(res, "data.frame")
  expect_equal(names(res), c("Test", "Basis", "Test_Statistic", "df", "Method", "p_value"))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_identical(res$Basis, "poly3")
})

test_that("def.gof accepts a model or (y, predicted_probs, X) identically", {
  fit <- make_fit(2)
  a <- def.gof(fit, basis = "poly2")
  b <- def.gof(as.numeric(fit$y), stats::fitted(fit),
               X = stats::model.matrix(fit), basis = "poly2")
  expect_equal(a$p_value, b$p_value, tolerance = 1e-8)
})

test_that("def.gof warns and stays computable without X", {
  fit <- make_fit(2)
  expect_warning(def.gof(as.numeric(fit$y), stats::fitted(fit)), "conservative")
})

test_that("def.gof detects a wrong link (cloglog truth, logit fit)", {
  res <- def.gof(make_fit(7, 1500, "cloglog"), basis = "poly2")
  expect_lt(res$p_value, 0.05)
})

test_that("all three bases run and give valid p-values", {
  fit <- make_fit()
  for (b in c("poly2", "poly3", "stukel")) {
    expect_true(def.gof(fit, basis = b)$p_value >= 0)
  }
})

test_that("def.gof errors on bad input", {
  set.seed(11)
  expect_error(def.gof(glm(rpois(30, 1) ~ rnorm(30), family = poisson())), "binomial")
  expect_error(def.gof(make_fit(), G = 2), "G' must be a single integer >= 3")
})
