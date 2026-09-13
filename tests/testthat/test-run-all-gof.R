make_fit <- function(seed = 1, n = 500) {
  set.seed(seed)
  x <- runif(n, -3, 3)
  glm(rbinom(n, 1, 1 / (1 + exp(-(0.6 * x)))) ~ x, family = binomial())
}

test_that("run.all.gof returns a tidy battery for a model", {
  res <- run.all.gof(make_fit(), include_slow = FALSE)
  expect_s3_class(res, "data.frame")
  expect_equal(names(res), c("Test", "Family", "Statistic", "df", "p_value", "Note"))
  expect_true(all(c("EF", "DEF.poly3", "DEF.sym", "HL", "Stukel") %in% res$Test))
  expect_true(any(grepl("Ensemble", res$Test)))
  pv <- res$p_value[is.finite(res$p_value)]
  expect_true(all(pv >= 0 & pv <= 1))
})

test_that("run.all.gof selects a subset", {
  res <- run.all.gof(make_fit(), tests = c("EF", "DEF.poly3", "HL"))
  expect_setequal(res$Test, c("EF", "DEF.poly3", "HL"))
})

test_that("prediction-only input runs the ph-only tests with a message", {
  fit <- make_fit(2)
  expect_message(res <- run.all.gof(as.numeric(fit$y), stats::fitted(fit), include_slow = FALSE),
                 "prediction-based")
  expect_true(all(is.na(res$p_value[res$Test %in% c("DEF.poly2", "Stukel")])))
  expect_true(is.finite(res$p_value[res$Test == "EF"]))
})

test_that("Osius-Rojek and Copas match the thesis sources (golden values)", {
  # Osius verified vs LogisticDx::gof.glm; Copas vs rms::lrm resid(., "gof").
  set.seed(42)
  n <- 400
  x1 <- rnorm(n); x2 <- runif(n, -2, 2)
  y  <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  res <- run.all.gof(fit, include_slow = FALSE)
  o <- res[res$Test == "Osius-Rojek", ]
  cps <- res[res$Test == "Copas-RSS", ]
  expect_equal(o$Statistic, 1.341183, tolerance = 1e-4)   # Osius z (LogisticDx form)
  expect_equal(o$p_value,  0.179861, tolerance = 1e-4)
  expect_equal(cps$p_value, 0.466596, tolerance = 1e-4)   # Copas p
})

test_that("Information-Matrix matches its source (golden)", {
  # Verified vs IMtest_fast (IM (infromation Matrix).R), diff = 0.
  set.seed(42)
  n <- 400
  x1 <- rnorm(n); x2 <- runif(n, -2, 2)
  y  <- rbinom(n, 1, plogis(0.3 + 0.8 * x1 - 0.5 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  im <- run.all.gof(fit, include_slow = FALSE)
  im <- im[im$Test == "Information-Matrix", ]
  expect_equal(im$Statistic, 1.198150, tolerance = 1e-4)
  expect_equal(im$df, 3)
  expect_equal(im$p_value, 0.753448, tolerance = 1e-4)
})

test_that("EF appears in both chisq and normal forms", {
  set.seed(8)
  n <- 400; x <- rnorm(n)
  fit <- glm(rbinom(n, 1, plogis(0.5 * x)) ~ x, family = binomial())
  res <- run.all.gof(fit, include_slow = FALSE)
  expect_true(all(c("EF", "EF-normal") %in% res$Test))
  expect_equal(res$p_value[res$Test == "EF-normal"],
               ef.gof(fit, method = "normal")$p_value, tolerance = 1e-8)
})

test_that("Pigeon-Heyse matches its source (deterministic golden)", {
  set.seed(77)
  n <- 600
  x1 <- rnorm(n); x2 <- runif(n, -2, 2)
  y <- rbinom(n, 1, plogis(0.1 + 0.8 * x1 - 0.4 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  ph <- run.all.gof(fit, include_slow = FALSE)
  ph <- ph[ph$Test == "Pigeon-Heyse", ]
  expect_equal(ph$Statistic, 8.658528, tolerance = 1e-4)
  expect_equal(ph$df, 9)
  expect_equal(ph$p_value, 0.469373, tolerance = 1e-4)
})

test_that("covariate-space tests appear and return valid results", {
  set.seed(5)
  n <- 400
  x1 <- rnorm(n); d <- factor(sample(c("A", "B"), n, replace = TRUE))
  y  <- rbinom(n, 1, plogis(0.3 + 0.6 * x1 + ifelse(d == "B", 0.5, 0)))
  fit <- glm(y ~ x1 + d, family = binomial())
  res <- run.all.gof(fit, include_slow = FALSE)
  expect_true(all(c("Tsiatis", "Xie", "Pulkstenis-Robinson") %in% res$Test))
  for (tt in c("Tsiatis", "Xie", "Pulkstenis-Robinson")) {
    p <- res$p_value[res$Test == tt]
    expect_true(is.na(p) || (p >= 0 && p <= 1))
  }
})

test_that("Pulkstenis-Robinson matches its source (deterministic golden)", {
  # PR uses no clustering, so this is platform-stable. Verified vs PR_test_only.R.
  set.seed(99)
  n <- 500
  x1 <- rnorm(n); x2 <- runif(n, -2, 2); d <- factor(sample(c("A", "B", "C"), n, replace = TRUE))
  eta <- 0.2 + 0.7 * x1 - 0.5 * x2 + ifelse(d == "B", 0.6, ifelse(d == "C", -0.4, 0))
  y <- rbinom(n, 1, plogis(eta))
  fit <- glm(y ~ x1 + x2 + d, family = binomial())
  pr <- run.all.gof(fit, include_slow = FALSE)
  pr <- pr[pr$Test == "Pulkstenis-Robinson", ]
  expect_equal(pr$Statistic, 1.910442, tolerance = 1e-4)
  expect_equal(pr$p_value,  0.591201, tolerance = 1e-4)
})

test_that("Pulkstenis-Robinson returns NA when no categorical covariate", {
  set.seed(3)
  n <- 300; x1 <- rnorm(n); x2 <- rnorm(n)
  fit <- glm(rbinom(n, 1, plogis(0.5 * x1)) ~ x1 + x2, family = binomial())
  res <- run.all.gof(fit, include_slow = FALSE)
  expect_true(is.na(res$p_value[res$Test == "Pulkstenis-Robinson"]))
})

test_that("Tier-2 tests are opt-in and return valid results", {
  skip_on_cran()                 # the full slow battery is heavy; keep CRAN fast
  skip_if_not_installed("mgcv")
  set.seed(202)
  n <- 300
  x1 <- runif(n, -3, 3); x2 <- rnorm(n); d <- factor(sample(c("A", "B"), n, replace = TRUE))
  y <- rbinom(n, 1, plogis(0.3 + 0.7 * x1 - 0.4 * x2 + ifelse(d == "B", 0.5, 0)))
  fit <- glm(y ~ x1 + x2 + d, family = binomial())
  slow_names <- c("HL-GAM", "PR-GAM", "Xie-GAM", "Stute-Zhu", "eHL", "BAGofT", "Lai-Liu-HL")
  expect_false(any(slow_names %in% run.all.gof(fit, include_slow = FALSE)$Test))   # not in default battery
  set.seed(1)
  res <- run.all.gof(fit, include_slow = TRUE,
                     control = list("Stute-Zhu" = list(B = 30), "BAGofT" = list(nsim = 10),
                                    "Lai-Liu-HL" = list(k = 30)))
  for (tt in slow_names) {
    p <- res$p_value[res$Test == tt]
    expect_true(length(p) == 1 && (is.na(p) || (p >= 0 && p <= 1)))
  }
  # Lai-Liu reports an accept/reject decision (no p-value)
  ll <- res[res$Test == "Lai-Liu-HL", ]
  expect_true(is.na(ll$p_value))
  expect_match(ll$Note, "decision:")
})

test_that("eHL internal reimplementation matches the source value", {
  # set.seed(5): .gof_ehl matches the marius-cp/eHL source eHL() to ~1e-11.
  set.seed(202)
  n <- 400; x <- runif(n, -3, 3)
  y <- rbinom(n, 1, plogis(0.3 + 0.7 * x))
  P <- plogis(0.3 + 0.7 * x)
  set.seed(5)
  expect_equal(ebrahim.gof:::.gof_ehl(y, P, boot = 10, s = 0.5), 0.068584, tolerance = 1e-3)
})

test_that("le Cessie is opt-in (slow) and uses the CORRECTED moment reference", {
  skip_on_cran()                 # include_slow = TRUE runs the full slow battery; keep CRAN fast
  set.seed(321)
  n <- 150
  x1 <- rnorm(n); x2 <- runif(n, -2, 2)
  y <- rbinom(n, 1, plogis(0.2 + 0.7 * x1 - 0.5 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  expect_false("le-Cessie" %in% run.all.gof(fit, include_slow = FALSE)$Test)         # not in the default battery
  res <- run.all.gof(fit, include_slow = TRUE)
  lc <- res[res$Test == "le-Cessie", ]
  ## These are the values AFTER the 2.4.1 transpose fix, and they deliberately do NOT
  ## match the inherited lecessie1995.r source: that implementation formed the moment
  ## reference as (I-H) R (I-H), and in a weighted fit H is not symmetric, so the correct
  ## form is (I-H)' R (I-H). Disagreeing with the source here is the point of the fix.
  ## The pre-fix values were Statistic 16.794627, p 0.149928.
  ## Note the reported Statistic is Q rescaled by the matched moments (Q * 2E[Q]/Var[Q]),
  ## so it moves with the correction even though the raw quadratic form Q does not.
  expect_equal(lc$Statistic, 15.895427, tolerance = 1e-4)
  expect_equal(lc$df,        11.648746, tolerance = 1e-4)
  expect_equal(lc$p_value,    0.177213, tolerance = 1e-4)
})

test_that("run.all.gof flags a cloglog link misfit", {
  set.seed(7)
  x <- runif(1500, -3, 3)
  fit <- glm(rbinom(1500, 1, 1 - exp(-exp(0.6 * x))) ~ x, family = binomial())
  res <- run.all.gof(fit, include_slow = FALSE)
  expect_lt(res$p_value[res$Test == "EF"], 0.05)
  expect_lt(res$p_value[res$Test == "Stukel"], 0.05)
})

## Stukel's two directions, built the way gof_stukel builds them, for the reference fits below
stukel_cols <- function(fit) {
  eta <- predict(fit, type = "link"); ph <- fitted(fit)
  cbind(za = 0.5 * eta^2 * (ph >= 0.5), zb = -0.5 * eta^2 * (ph < 0.5))
}

test_that("Stukel is the joint score test and agrees with anova(test = 'Rao')", {
  fit <- make_fit()
  res <- run.all.gof(fit, tests = "Stukel")
  expect_equal(res$Statistic, 1.270855, tolerance = 1e-5)
  expect_equal(res$df, 2)
  expect_equal(res$p_value, 0.5297091, tolerance = 1e-5)
  expect_identical(res$Note, "")
  Z <- stukel_cols(fit); y <- fit$y; X <- model.matrix(fit)
  rao <- anova(glm(y ~ X - 1, family = binomial()), glm(y ~ X + Z - 1, family = binomial()),
               test = "Rao")$Rao[2]
  expect_equal(res$Statistic, rao, tolerance = 1e-5)
})

test_that("the marginal Stukel sum is kept for reproducing 2.7.0 and earlier", {
  skip_if_not_installed("statmod")
  res <- run.all.gof(make_fit(), tests = "Stukel", control = list(Stukel = list(form = "marginal")))
  expect_equal(res$Statistic, 0.394475, tolerance = 1e-5)
  expect_equal(res$df, 2)
  expect_equal(res$p_value, 0.8209956, tolerance = 1e-5)
  expect_match(res$Note, "marginal")
})

test_that("Stukel falls back to one df when every fitted risk is on one side of 0.5", {
  set.seed(5); x <- runif(400, -3, 3); y <- rbinom(400, 1, plogis(-4 + 0.5 * x))
  lo <- run.all.gof(glm(y ~ x, family = binomial()), tests = "Stukel")
  expect_equal(lo$df, 1)
  expect_equal(lo$Statistic, 0.6245748, tolerance = 1e-5)
  expect_equal(lo$p_value, 0.4293523, tolerance = 1e-5)
  expect_match(lo$Note, "no fitted risk at or above 0.5")
  set.seed(6); x <- runif(400, -3, 3); y <- rbinom(400, 1, plogis(4 + 0.5 * x))
  hi <- run.all.gof(glm(y ~ x, family = binomial()), tests = "Stukel")
  expect_equal(hi$df, 1)
  expect_equal(hi$Statistic, 0.0356688, tolerance = 1e-4)
  expect_equal(hi$p_value, 0.8502010, tolerance = 1e-5)
  expect_match(hi$Note, "no fitted risk below 0.5")
})

test_that("Stukel form = 'lr' is the deviance drop of the augmented fit", {
  fit <- make_fit()
  res <- run.all.gof(fit, tests = "Stukel", control = list(Stukel = list(form = "lr")))
  Z <- stukel_cols(fit); y <- fit$y; X <- model.matrix(fit)
  f1 <- glm(y ~ X + Z - 1, family = binomial())
  expect_equal(res$Statistic, deviance(fit) - deviance(f1), tolerance = 1e-6)
  expect_equal(res$df, 2)
  expect_equal(res$p_value, pchisq(res$Statistic, 2, lower.tail = FALSE), tolerance = 1e-10)

  set.seed(5); x <- runif(400, -3, 3); y <- rbinom(400, 1, plogis(-4 + 0.5 * x))
  lo  <- glm(y ~ x, family = binomial())
  zb  <- stukel_cols(lo)[, "zb"]
  res <- run.all.gof(lo, tests = "Stukel", control = list(Stukel = list(form = "lr")))
  expect_equal(res$df, 1)
  expect_equal(res$Statistic, deviance(lo) - deviance(glm(y ~ x + zb, family = binomial())),
               tolerance = 1e-6)
  expect_match(res$Note, "0.5")
})

test_that("the joint and LR Stukel forms refuse a non-logit fit", {
  set.seed(3); x <- runif(300, -3, 3)
  fit <- glm(rbinom(300, 1, pnorm(0.5 * x)) ~ x, family = binomial("probit"))
  for (f in c("joint", "lr")) {
    res <- run.all.gof(fit, tests = "Stukel", control = list(Stukel = list(form = f)))
    expect_true(is.na(res$p_value))
    expect_match(res$Note, "logit")
  }
})
