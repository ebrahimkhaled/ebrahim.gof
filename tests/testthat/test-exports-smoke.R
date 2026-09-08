# B3-03: ten exports had no test at all, calm.gof among them, and it is the function this
# release exists to add. These are smoke tests, not a substitute for the simulation studies
# behind each method: they assert the contract each function documents -- the class and fields
# it returns, that p-values are probabilities, that the documented guards fire -- so that a
# refactor cannot silently change the shape of what a user gets back.

test_that("calm.gof returns the documented object and a valid p-value", {
  skip_if_not_installed("CompQuadForm")
  set.seed(11)
  n <- 240; p <- 6
  X <- matrix(rnorm(n * p), n, p)
  y <- rbinom(n, 1, plogis(0.2 + X %*% c(0.8, -0.5, 0.3, 0, 0, 0)))

  res <- calm.gof(X, y, lambda = 30)
  expect_s3_class(res, "calm.gof")
  expect_true(all(c("SC.HL", "SC.EDGE.adaptive", "lambda", "lambda_glmnet",
                    "G", "kappa", "prevalence", "rho_hat") %in% names(res)))
  for (nm in c("SC.HL", "SC.EDGE", "SC.EDGE.adaptive")) {
    expect_gte(res[[nm]]$p.value, 0)
    expect_lte(res[[nm]]$p.value, 1)
  }
  expect_equal(res$kappa, p / n)
  expect_equal(res$prevalence, mean(y))
  # the two penalty scales differ by exactly n, which is the commonest error with this method
  expect_equal(res$lambda, res$lambda_glmnet * n)
  expect_output(print(res), "CALM")
})

test_that("calm.gof refuses p >= n and non-binary responses", {
  skip_if_not_installed("CompQuadForm")
  set.seed(12)
  X <- matrix(rnorm(20 * 25), 20, 25)
  y <- rbinom(20, 1, 0.5)
  expect_error(calm.gof(X, y, lambda = 10), "p < n")

  X2 <- matrix(rnorm(100 * 3), 100, 3)
  expect_error(calm.gof(X2, rnorm(100), lambda = 10), "only 0 and 1")
})

test_that("calm.gof warns when the outcome is strongly unbalanced", {
  skip_if_not_installed("CompQuadForm")
  set.seed(13)
  n <- 300; p <- 20
  X <- matrix(rnorm(n * p), n, p)
  y <- rbinom(n, 1, plogis(-2.6 + X %*% c(0.8, -0.5, 0.3, rep(0, p - 3))))
  expect_lt(mean(y), 0.35)                      # the design really is unbalanced
  expect_warning(calm.gof(X, y, lambda = 40), "prevalence")
})

test_that("shrink.gof returns the documented object and prints as a method", {
  set.seed(14)
  n <- 200
  X <- matrix(rnorm(n * 4), n, 4)
  y <- rbinom(n, 1, plogis(0.3 + X %*% c(0.8, -0.5, 0.3, 0)))

  res <- shrink.gof(X, y, lambda = 20, basis = "decile", B = 39, seed = 1)
  expect_s3_class(res, "shrink.gof")
  expect_gte(res$SC.HL$p.value, 0)
  expect_lte(res$SC.HL$p.value, 1)
  # B1-06: the print method is registered, so this dispatches rather than
  # falling through to print.default and dumping the raw list
  expect_output(print(res), "Shrinkage-corrected")
})

test_that("the directed tests fire on gof_demo's documented misfit and stand down when corrected", {
  data("gof_demo", package = "ebrahim.gof")
  wrong <- glm(outcome ~ age + bmi + sex + treatment,
               data = gof_demo, family = binomial())
  right <- glm(outcome ~ poly(age, 2) + bmi + sex + treatment,
               data = gof_demo, family = binomial())

  expect_lt(def.gof(wrong)$p_value, 0.05)
  expect_gt(def.gof(right)$p_value, 0.05)
  expect_lt(def.ensemble.gof(wrong)$p_value, 0.05)
})

test_that("edge.gof and edges.gof return probabilities on a fitted model", {
  data("gof_demo", package = "ebrahim.gof")
  wrong <- glm(outcome ~ age + bmi + sex + treatment,
               data = gof_demo, family = binomial())
  for (f in list(edge.gof, edges.gof)) {
    p <- f(wrong)$p_value
    expect_true(is.numeric(p) && length(p) == 1)
    expect_gte(p, 0); expect_lte(p, 1)
  }
})

test_that("cdef.gof returns a probability on both bases", {
  data("gof_demo", package = "ebrahim.gof")
  wrong <- glm(outcome ~ age + bmi + sex + treatment,
               data = gof_demo, family = binomial())
  for (b in c("poly", "spline")) {
    p <- suppressWarnings(cdef.gof(wrong, basis = b)$p_value)
    expect_gte(p, 0); expect_lte(p, 1)
  }
})

test_that("gof.features honours its tests argument (B1-08) and returns the documented length", {
  data("gof_demo", package = "ebrahim.gof")
  fit <- glm(outcome ~ age + bmi + sex + treatment,
             data = gof_demo, family = binomial())

  full <- suppressWarnings(gof.features(fit))
  expect_length(full, 11)                       # nine panel tests + CDEF.poly + CDEF.spline
  expect_true(all(is.finite(full)))

  # a subset really is a subset: the argument is used, not ignored
  few <- suppressWarnings(gof.features(fit, tests = c("HL", "EF")))
  expect_length(few, 4)                         # two panel tests + the two CDEF entries
  expect_true(all(c("HL", "EF", "CDEF.poly", "CDEF.spline") %in% names(few)))
})

test_that("deploy.gof accepts a user-supplied feature function", {
  data("gof_demo", package = "ebrahim.gof")
  fit <- glm(outcome ~ age + bmi + sex + treatment,
             data = gof_demo, family = binomial())
  # feature_fn is the documented hook; the default resolves gof.features in the namespace
  # meta is the pre-trained scorer: here the simplest one that satisfies the contract,
  # a function of the feature vector returning a scalar misfit score
  out <- suppressWarnings(deploy.gof(fit, meta = function(feat) max(abs(feat)), B = 19,
                                     feature_fn = function(f) gof.features(f, tests = c("HL", "EF"))))
  expect_s3_class(out, "data.frame")
  expect_true(all(c("Test", "Score", "B", "p_value") %in% names(out)))
  expect_gte(out$p_value, 0); expect_lte(out$p_value, 1)
})
