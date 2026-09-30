test_that("run.all.external returns one row per test in the battery format", {
  set.seed(11)
  n <- 600
  x <- rnorm(n)
  p <- plogis(-0.8 + 0.7 * x)
  y <- rbinom(n, 1, p)
  out <- run.all.external(y, p)
  expect_s3_class(out, "gof_battery")
  expect_named(out, c("Test", "Family", "Statistic", "df", "p_value", "Note"))
  expect_true(all(c("EDGE (G=10)", "Cox recalibration", "Calibration-in-the-large", "Spiegelhalter-z",
                    "Hosmer-Lemeshow (external)", "Stukel (offset)", "O/E ratio", "Calibration slope",
                    "c-statistic (AUC)") %in% out$Test))
  expect_true(any(grepl("^EDGE \\(G=auto, 24\\)$", out$Test)))
  expect_false("le Cessie (external)" %in% out$Test)
  inferential <- out[out$Family != "Descriptive" & out$Test != "GiViTI", ]
  expect_true(all(is.finite(inferential$p_value)))
  expect_true(all(inferential$p_value >= 0 & inferential$p_value <= 1))
})

test_that("the directed rows equal def.gof(external = TRUE)", {
  set.seed(12)
  n <- 700
  x <- rnorm(n)
  p <- plogis(-1 + 0.9 * x)
  y <- rbinom(n, 1, plogis(-1 + 0.6 * x))
  out <- run.all.external(y, p)
  e10 <- def.gof(y, predicted_probs = p, G = 10, external = TRUE)
  ea  <- def.gof(y, predicted_probs = p, G = "auto", external = TRUE)
  expect_equal(out$Statistic[out$Test == "EDGE (G=10)"], e10$Test_Statistic, tolerance = 1e-12)
  expect_equal(out$p_value[out$Test == "EDGE (G=10)"], e10$p_value, tolerance = 1e-12)
  expect_equal(out$p_value[grepl("G=auto", out$Test)], ea$p_value, tolerance = 1e-12)
})

test_that("the closed-form rows match hand computations", {
  set.seed(13)
  n <- 500
  p <- plogis(rnorm(n, -0.7, 0.9))
  y <- rbinom(n, 1, p)
  out <- run.all.external(y, p)
  lp <- qlogis(p); v <- p * (1 - p)
  dev0 <- -2 * sum(y * log(p) + (1 - y) * log(1 - p))
  cox <- dev0 - glm(y ~ lp, family = binomial())$deviance
  expect_equal(out$Statistic[out$Test == "Cox recalibration"], cox, tolerance = 1e-10)
  zs <- sum((y - p) * (1 - 2 * p)) / sqrt(sum((1 - 2 * p)^2 * v))
  expect_equal(out$p_value[out$Test == "Spiegelhalter-z"], 2 * pnorm(-abs(zs)), tolerance = 1e-12)
  g <- pmin(ceiling(rank(p, ties.method = "first") / (n / 10)), 10)
  O <- tapply(y, g, sum); E <- tapply(p, g, sum); ng <- tabulate(g, 10)
  expect_equal(out$Statistic[out$Test == "Hosmer-Lemeshow (external)"],
               sum((O - E)^2 / (E * (1 - E / ng))), tolerance = 1e-10)
  expect_equal(out$df[out$Test == "Hosmer-Lemeshow (external)"], 10)
  expect_equal(out$Statistic[out$Test == "O/E ratio"], sum(y) / sum(p), tolerance = 1e-12)
})

test_that("an overfitted model is flagged by the recalibration and directed tests", {
  set.seed(14)
  n <- 3000
  x <- rnorm(n)
  p <- plogis(-0.5 + 1.2 * x)                  # published slope 1.2
  y <- rbinom(n, 1, plogis(-0.5 + 0.6 * x))    # true slope 0.6
  out <- run.all.external(y, p)
  expect_lt(out$p_value[out$Test == "Cox recalibration"], 1e-6)
  expect_lt(out$p_value[out$Test == "EDGE (G=10)"], 1e-3)
  expect_lt(out$Statistic[out$Test == "Calibration slope"], 0.8)
})

test_that("le Cessie's external row runs only on request", {
  set.seed(15)
  n <- 300
  X <- data.frame(a = rnorm(n), b = rnorm(n))
  p <- plogis(-0.5 + 0.7 * X$a)
  y <- rbinom(n, 1, p)
  held <- run.all.external(y, p, X = X)
  expect_true(is.na(held$p_value[held$Test == "le Cessie (external)"]))
  expect_match(held$Note[held$Test == "le Cessie (external)"], "include_slow")
  ran <- run.all.external(y, p, X = X, include_slow = TRUE)
  pv <- ran$p_value[ran$Test == "le Cessie (external)"]
  expect_true(is.finite(pv) && pv > 0 && pv < 1)
})

test_that("bad input is refused", {
  expect_error(run.all.external(c(0, 1, 2), c(.2, .5, .7)), "binary")
  expect_error(run.all.external(c(0, 1), c(.2, .5, .7)), "lengths")
  expect_error(run.all.external(c(0, 1), c(0, .5)), "strictly between")
  expect_error(run.all.external(c(0, 1), c(.3, .5), X = matrix(1:6, 3)), "one row")
})
