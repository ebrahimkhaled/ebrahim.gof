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

test_that("a near-empty Stukel half-column is kept and scaled, not left to a singular solve", {
  set.seed(4)
  ph <- c(sort(runif(475, 0.02, 0.45)), rep(0.5001, 25))            # G = 20: the top group's mean risk is 0.5001
  X  <- cbind(1, qlogis(ph)); y <- rbinom(500, 1, ph)
  res <- def.gof(y, ph, X = X, G = 20, basis = "stukel")
  expect_true(is.finite(res$p_value))
  ## the unit form built by hand with all three columns kept; the e^2 (e >= 0) column is about 1.6e-7
  ## long, so each column is scaled to unit length before the solve (S does not depend on that scale)
  g  <- pmin(ceiling(rank(ph, ties.method = "first") / 25), 20)
  V  <- ph * (1 - ph); Vg <- as.numeric(tapply(V, g, sum))
  r  <- as.numeric(tapply(y - ph, g, sum)) / sqrt(Vg)
  e  <- qlogis(as.numeric(tapply(ph, g, mean)))
  Z  <- cbind(e, e^2 * (e >= 0), -e^2 * (e < 0))
  Z  <- sweep(Z, 2, sqrt(colSums(Z^2)), "/")
  U  <- rowsum(V * X, g) / sqrt(Vg)
  Om <- diag(20) - U %*% solve(crossprod(X, V * X), t(U))
  S  <- drop(crossprod(r, Z %*% solve(crossprod(Z), crossprod(Z, r))))
  lam <- Re(eigen(solve(crossprod(Z), crossprod(Z, Om %*% Z)), only.values = TRUE)$values)
  lam <- lam[lam > 1e-9]
  p  <- pchisq(S / (sum(lam^2) / sum(lam)), sum(lam)^2 / sum(lam^2), lower.tail = FALSE)
  expect_equal(res$Test_Statistic, S, tolerance = 1e-10)
  expect_equal(res$p_value, p, tolerance = 1e-10)
  expect_equal(def.gof(y, ph, X = X, G = 20, basis = "stukel", weights = "score")$df, 3)

  ## a glm fit whose only group at or above 0.5 has mean risk 0.5008: the score form keeps all three
  ## columns and equals anova(..., test = "Rao") for the grouped step covariates, on a warm-restarted fit
  set.seed(392)
  dat <- data.frame(xa = runif(500, -3, 3), db = rbinom(500, 1, 0.5))
  dat$out <- rbinom(500, 1, plogis(-2 + 0.6 * dat$xa + 0.5 * dat$db))
  ctl <- glm.control(epsilon = 1e-13, maxit = 100)
  fit <- glm(out ~ xa + db, family = binomial(), data = dat, control = ctl)
  fit <- glm(out ~ xa + db, family = binomial(), data = dat, control = ctl, start = coef(fit))
  pf <- pmin(pmax(fitted(fit), 1e-6), 1 - 1e-6)
  gf <- pmin(ceiling(rank(pf, ties.method = "first") / 25), 20)
  ef <- qlogis(as.numeric(tapply(pf, gf, mean)))
  expect_equal(sum(ef >= 0), 1)
  s   <- cbind(ef, ef^2 * (ef >= 0), -ef^2 * (ef < 0))[gf, ]
  aug <- suppressWarnings(glm(out ~ xa + db + s, family = binomial(), data = dat, control = ctl))
  rao <- anova(fit, aug, test = "Rao")
  sc  <- def.gof(fit, G = 20, basis = "stukel", weights = "score")
  expect_equal(rao$Df[2], 3)
  expect_equal(sc$df, 3)
  expect_equal(sc$Test_Statistic, rao$Rao[2], tolerance = 1e-8)
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

test_that("def.gof warns, but still returns a p-value, with fewer events than groups", {
  set.seed(9)
  x <- sort(runif(300, -3, 3))
  y <- integer(300); y[sample(251:300, 7)] <- 1L                   # 7 events
  fit <- glm(y ~ x, family = binomial())
  expect_warning(res <- def.gof(fit), "unreliable with fewer events than groups")
  expect_true(is.finite(res$p_value))
  expect_silent(def.gof(fit, G = 5))
  yy <- 1L - y
  expect_warning(def.gof(glm(yy ~ x, family = binomial()), basis = "sym"), "fewer non-events than groups")

  count <- function(expr) {
    k <- 0L
    suppressWarnings(withCallingHandlers(expr, def_few_events = function(w) k <<- k + 1L))
    k
  }
  expect_equal(count(def.ensemble.gof(fit)), 1L)                    # once, not once per basis
  expect_equal(count(bat <- run.all.gof(fit, include_slow = FALSE, install = "no")), 0L)
  expect_match(bat$Note[bat$Test == "DEF.poly3"], "fewer events than groups")
})

test_that("a sample with no event or no non-event gives NA, with a def_degenerate warning", {
  set.seed(20)
  x1 <- as.numeric(scale(rchisq(200, 4))); x2 <- runif(200, -3, 3); d2 <- rbinom(200, 1, 0.5)
  sets <- list(all0 = data.frame(x = x1, y = rep(0, 200)),
               all1 = data.frame(x = x1, y = rep(1, 200)),
               all0_two = data.frame(x = x2, d = d2, y = rep(0, 200)))
  for (nm in names(sets)) {
    fit <- suppressWarnings(glm(y ~ ., family = binomial(), data = sets[[nm]]))
    for (b in c("poly3", "poly2", "stukel", "sym")) for (w in c("unit", "score")) {
      expect_warning(res <- def.gof(fit, basis = b, weights = w), class = "def_degenerate")
      expect_true(is.na(res$p_value))
      expect_identical(res$Method, if (w == "score") "score" else "satterthwaite")
    }
    expect_warning(res <- def.gof(as.numeric(fit$y), fitted(fit), X = model.matrix(fit), basis = "stukel",
                                  weights = "score"), class = "def_degenerate")
    expect_true(is.na(res$p_value))
    expect_warning(res <- edge.gof(fit, basis = "sym", weights = "score"), class = "def_degenerate")
    expect_true(is.na(res$p_value))
    k <- 0L
    en <- withCallingHandlers(def.ensemble.gof(fit, add_ef = TRUE),
                              def_degenerate = function(w) { k <<- k + 1L; invokeRestart("muffleWarning") })
    expect_equal(k, 1L)                                            # once, not once per basis
    expect_true(is.na(en$p_value))
    for (f in c("joint", "lr", "marginal")) {
      st <- run.all.gof(fit, tests = "Stukel", install = "no", control = list(Stukel = list(form = f)))
      expect_true(is.na(st$p_value))
      expect_match(st$Note, "no (events|non-events) in the sample")
    }
    dr <- run.all.gof(fit, tests = "DEF.stukel", install = "no")
    expect_true(is.na(dr$p_value))
    expect_match(dr$Note, "no maximum-likelihood fit")
  }
})

test_that("def.gof errors on bad input", {
  set.seed(11)
  expect_error(def.gof(glm(rpois(30, 1) ~ rnorm(30), family = poisson())), "binomial")
  expect_error(def.gof(make_fit(), G = 2), "G' must be a single integer >= 3")
})
