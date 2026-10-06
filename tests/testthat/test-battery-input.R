## Inputs outside the scope of the battery, the rows it leaves out, and the width of what it prints.

sim_fit <- function(n = 600, seed = 11) {
  set.seed(seed)
  d <- data.frame(x1 = rnorm(n), x2 = runif(n, -2, 2))
  d$y <- rbinom(n, 1, plogis(-0.3 + 0.8 * d$x1 - 0.5 * d$x2))
  glm(y ~ x1 + x2, family = binomial(), data = d)
}

test_that("weighted and grouped (cbind) fits are refused, another link warns", {
  f <- sim_fit()
  d <- f$data
  w <- glm(y ~ x1 + x2, family = binomial(), data = d, weights = rep(1:3, length.out = nrow(d)))
  expect_error(run.all.gof(w, include_slow = FALSE), "one row per binary")
  set.seed(2)
  g <- data.frame(x = rep(seq(-2, 2, length.out = 40), each = 1), n = 20)
  g$s <- rbinom(40, 20, plogis(0.5 * g$x))
  cb <- glm(cbind(s, n - s) ~ x, family = binomial(), data = g)
  expect_error(run.all.gof(cb, include_slow = FALSE), "one row per binary")
  pr <- glm(y ~ x1 + x2, family = binomial(link = "probit"), data = d)
  expect_warning(run.all.gof(pr, tests = "HL"), "probit")
})

test_that("separation and very few events are flagged", {
  set.seed(5)
  x <- c(rnorm(40, -3), rnorm(40, 3)); y <- rep(0:1, each = 40)
  sep <- suppressWarnings(glm(y ~ x, family = binomial()))
  expect_warning(run.all.gof(sep, tests = "HL"), "separated")
  small <- sim_fit(n = 30, seed = 3)
  expect_warning(run.all.gof(small, tests = "HL"), "events")
})

test_that("the (y, p, X) form with a data frame reproduces the glm path", {
  f <- sim_fit(n = 1000)
  X <- f$data[c("x1", "x2")]
  a <- run.all.gof(f, tests = c("DEF.poly2", "DEF.stukel", "DEF.sym", "Tsiatis"), install = "no")
  b <- run.all.gof(f$y, f$fitted.values, X = X,
                   tests = c("DEF.poly2", "DEF.stukel", "DEF.sym", "Tsiatis"), install = "no")
  expect_equal(b$p_value, a$p_value, tolerance = 1e-8)
  m <- run.all.gof(f$y, f$fitted.values, X = as.matrix(X), tests = "DEF.sym", install = "no")
  expect_equal(m$p_value, a$p_value[a$Test == "DEF.sym"], tolerance = 1e-8)
})

test_that("the full battery leaves out repeated rows and GiViTI-external", {
  f <- sim_fit(n = 600)
  res <- run.all.gof(f, include_slow = FALSE, install = "no")
  expect_false("DEF.poly3" %in% res$Test)             # = EDGE.G10 at G = 10
  expect_true("EDGE.G10" %in% res$Test)
  expect_true("DEF.poly3" %in% run.all.gof(f, include_slow = FALSE, install = "no", G = 12)$Test)
  expect_true("DEF.poly3" %in% run.all.gof(f, tests = "DEF.poly3")$Test)
  small <- run.all.gof(sim_fit(n = 200), include_slow = FALSE, install = "no")
  expect_false("EDGE.G10" %in% small$Test)              # the default partition has ten groups
  expect_false("GiViTI-external" %in% names(Filter(function(e) !isTRUE(e$on_request),
                                                    ebrahim.gof:::.GOF_REGISTRY)))
  expect_false(any(grepl(" \\(grouped\\)$", res$Test)))  # every record its own pattern
})

test_that("comparison rows are not counted, and EF-normal runs on many groups", {
  f <- sim_fit(n = 600)
  res <- run.all.gof(f, include_slow = FALSE, install = "no")
  out <- capture.output(print(res))
  k <- sum(res$p_value < 0.05 & !res$Test %in% c("Pearson", "Deviance", "F-test"), na.rm = TRUE)
  expect_match(out[2], sprintf("(%d reject at 0.05)", k), fixed = TRUE)
  s <- summary(res)
  expect_false(any(s$Smallest_test %in% c("Pearson", "Deviance", "F-test")))
  expect_identical(res$Family[res$Test == "EF-normal"], "Partition")
  expect_match(res$Note[res$Test == "EF-normal"], "G = 24")
  expect_true(is.na(res$df[res$Test == "Copas-RSS"]))
  one <- capture.output(print(run.all.gof(f, tests = "HL")))
  expect_true(any(grepl("battery: 1 test ", one, fixed = TRUE)))
})

test_that("printed lines fit a 70-character console", {
  f <- sim_fit(n = 3000)
  old <- options(width = 70); on.exit(options(old))
  res <- run.all.gof(f, include_slow = FALSE, install = "no")
  out <- c(capture.output(print(res)), capture.output(print(summary(res))))
  expect_true(all(nchar(out) <= 70))
  expect_true(any(grepl("Family.*Given by", out)))        # the summary table is not wrapped
  long <- structure(data.frame(Test = c("Hosmer-Lemeshow (external)", "Ensemble.Univ(3DEF+EF)"),
                               Family = c("Covariate-space", "Ensemble"), Statistic = NA_real_,
                               df = NA_real_, p_value = c(1e-5, 0.0123), Note = ""),
                    class = c("gof_battery", "data.frame"))
  lo <- capture.output(print(summary(long)))
  expect_true(any(grepl("Family.*Given by", lo)))
  expect_true(all(nchar(lo) <= 70))
  expect_true(any(grepl(" 2997 ", out, fixed = TRUE)))   # whole-number df in full
})

test_that("Stute-Zhu never reports a p-value of 0", {
  set.seed(1)
  x <- runif(400, -3, 3); y <- rbinom(400, 1, plogis(x + 0.8 * x^2))
  bad <- glm(y ~ x, family = binomial())
  set.seed(1)
  r <- run.all.gof(bad, tests = "Stute-Zhu", control = list("Stute-Zhu" = list(B = 19)))
  expect_gt(r$p_value, 0)                           # (1 + 0) / (1 + valid replicates)
  expect_lte(r$p_value, 1 / 15)
})

test_that("a seed argument leaves the caller's random number stream alone", {
  f <- sim_fit(n = 300)
  set.seed(10); invisible(runif(1)); ref <- runif(3)
  set.seed(10); invisible(runif(1))
  invisible(localize.gof(f, B = 19, seed = 1, plot = FALSE))
  expect_identical(runif(3), ref)
  set.seed(10); invisible(runif(1))
  invisible(localize.external(f$y, f$fitted.values * 0.9, f$data[c("x1", "x2")], M = 19, seed = 1,
                              plot = FALSE))
  expect_identical(runif(3), ref)
})

test_that("localize.external() warns on a model's own fitted values", {
  f <- sim_fit(n = 300)
  expect_warning(localize.external(f$y, f$fitted.values, f$data[c("x1", "x2")], M = 19, seed = 1,
                                   plot = FALSE), "localize.gof")
})

test_that("calm.gof() checks lambda and missing values", {
  skip_if_not_installed("CompQuadForm")
  set.seed(4)
  X <- matrix(rnorm(200 * 10), 200, 10); y <- rbinom(200, 1, plogis(X[, 1]))
  expect_error(calm.gof(X, y, lambda = -1), "non-negative")
  expect_error(calm.gof(X, y, lambda = c(1, 2)), "single")
  X[3, 2] <- NA
  expect_error(calm.gof(X, y, lambda = 10), "missing")
  expect_error(shrink.gof(X, y, lambda = 10, B = 19), "missing")
})

test_that("plot() of a single row taken from edge.gof() draws it", {
  f <- sim_fit(n = 600)
  e <- edge.gof(f)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(plot(e[2, ]))
  expect_error(plot(e[2, ], row = "verdict"), "no 'verdict' row")
})

test_that("the GiViTI row leaves the session's random number stream alone", {
  skip_on_cran()
  skip_if_not_installed("givitiR")
  skip_if_not_installed("callr")
  f <- sim_fit(n = 300)
  set.seed(3); ref <- runif(2)
  set.seed(3); invisible(run.all.gof(f, tests = "GiViTI"))
  expect_identical(runif(2), ref)
})
