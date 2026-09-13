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

test_that("all four bases run and give valid p-values", {
  fit <- make_fit()
  for (b in c("poly2", "poly3", "stukel", "sym")) {
    expect_true(def.gof(fit, basis = b)$p_value >= 0)
  }
})

test_that("basis 'sym' and the score form on gof_demo", {
  data("gof_demo", package = "ebrahim.gof")
  wrong <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial())
  right <- glm(outcome ~ poly(age, 2) + bmi + sex + treatment, data = gof_demo, family = binomial())
  u <- def.gof(wrong, basis = "sym")
  expect_equal(u$Test_Statistic, 1.5188084, tolerance = 1e-6)
  expect_equal(u$df, 1)
  expect_equal(u$p_value, 0.0095796, tolerance = 1e-4)
  s <- def.gof(wrong, basis = "sym", weights = "score")
  expect_equal(names(s), c("Test", "Basis", "Test_Statistic", "df", "Method", "p_value"))
  expect_identical(s$Method, "score")
  expect_equal(s$Test_Statistic, 8.4571259, tolerance = 1e-6)
  expect_equal(s$df, 1)
  expect_equal(s$p_value, 0.0036362, tolerance = 1e-4)
  p3 <- def.gof(wrong, weights = "score")
  expect_equal(p3$Test_Statistic, 10.4736625, tolerance = 1e-6)
  expect_equal(p3$df, 3)
  expect_equal(p3$p_value, 0.0149406, tolerance = 1e-4)
  st <- def.gof(wrong, basis = "stukel", weights = "score")        # the same three columns as the unit form
  expect_equal(st$Test_Statistic, 8.6849505, tolerance = 1e-6)
  expect_equal(st$df, 3)
  expect_gt(def.gof(right, basis = "sym")$p_value, 0.5)                       # 0.5803836
  expect_gt(def.gof(right, basis = "sym", weights = "score")$p_value, 0.4)    # 0.4437683
})

test_that("the score form is the Rao score test for the grouped step covariate", {
  data("gof_demo", package = "ebrahim.gof")
  fit <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial())
  ph <- pmin(pmax(fitted(fit), 1e-6), 1 - 1e-6); n <- length(ph)
  g <- pmin(ceiling(rank(ph, ties.method = "first") / (n / 10)), 10)
  e <- qlogis(as.numeric(tapply(ph, g, mean))); s <- (e * abs(e))[g]
  y <- fit$y; X <- model.matrix(fit)
  rao <- anova(glm(y ~ X - 1, family = binomial()), glm(y ~ X + s - 1, family = binomial()),
               test = "Rao")$Rao[2]
  expect_equal(def.gof(fit, basis = "sym", weights = "score")$Test_Statistic, rao, tolerance = 1e-4)
})

test_that("the score form is identical from a glm and from (y, predicted_probs, X)", {
  fit <- make_fit(2)
  for (b in c("sym", "poly3")) {
    a <- def.gof(fit, basis = b, weights = "score")
    d <- def.gof(as.numeric(fit$y), stats::fitted(fit), X = stats::model.matrix(fit),
                 basis = b, weights = "score")
    expect_equal(a$p_value, d$p_value, tolerance = 1e-10)
  }
})

test_that("G = 'auto' is max(10, round(n / 25))", {
  fit <- make_fit()                                                 # n = 600, so G = 24
  expect_identical(def.gof(fit, G = "auto"), def.gof(fit, G = 24))
  expect_identical(edge.gof(fit, G = "auto", basis = "sym"), edge.gof(fit, G = 24, basis = "sym"))
  expect_identical(def.gof(make_fit(n = 200), G = "auto"), def.gof(make_fit(n = 200), G = 10))
  expect_error(def.gof(fit, G = "many"), "or 'auto'")
  res <- run.all.gof(fit, tests = "DEF.sym", control = list(DEF.sym = list(weights = "score", G = "auto")))
  expect_equal(res$p_value, def.gof(fit, G = 24, basis = "sym", weights = "score")$p_value)
  expect_match(res$Note, "G = 24 \\(auto\\)")
  expect_match(res$Note, "score form")
})

test_that("def.gof errors on bad input", {
  set.seed(11)
  expect_error(def.gof(glm(rpois(30, 1) ~ rnorm(30), family = poisson())), "binomial")
  expect_error(def.gof(make_fit(), G = 2), "G' must be a single integer >= 3")
})
