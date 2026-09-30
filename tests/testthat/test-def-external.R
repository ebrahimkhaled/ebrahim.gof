## External mode (frozen predictions) and the G = "auto" rule.

## The reference implementation of the EDGE paper-3 external mode (simulations/edge_external.R), with one
## line changed: the groups. The reference cuts at quantile(type = 1) breaks, which puts 101 and 99 rows in
## some groups at n = 1000, G = 10 (n * k / G is not exact in floating point); the package keeps its own
## equal-frequency rule, ceiling(rank / (n / G)), in both modes. Everything after the grouping is verbatim.
edge_external_ref <- function(y, p, G = 10, basis = c("poly3", "stukel", "sym")) {
  basis <- match.arg(basis)
  stopifnot(length(y) == length(p), all(y %in% c(0, 1)), all(p > 0 & p < 1))
  if (identical(G, "auto")) G <- max(10L, as.integer(ceiling(length(y) / 25)))
  n <- length(y)
  g <- pmin(ceiling(rank(p, ties.method = "first") / (n / G)), G)          # the package's groups
  M <- rowsum(cbind(y, p, p * (1 - p), 1), g, reorder = TRUE)
  r <- (M[, 1] - M[, 2]) / sqrt(M[, 3])
  pb <- M[, 2] / M[, 4]
  eta <- stats::qlogis(pb)
  Z <- switch(basis,
    poly3  = cbind(1, pb, pb^2, pb^3),
    stukel = cbind(1, eta, eta^2 * (eta >= 0), -eta^2 * (eta < 0)),
    sym    = cbind(1, eta * abs(eta)))
  Q <- qr(Z)
  Z <- Z[, Q$pivot[seq_len(Q$rank)], drop = FALSE]
  S <- sum(qr.fitted(qr(Z), r)^2)
  list(statistic = S, df = ncol(Z), p_value = stats::pchisq(S, ncol(Z), lower.tail = FALSE))
}

test_that("external mode reproduces the reference statistic, df and p-value", {
  for (n in c(1000, 1234)) {
    set.seed(n)
    p <- plogis(rnorm(n, -0.5, 1.2))
    y <- rbinom(n, 1, plogis(1.2 * qlogis(p)))                      # a frozen model that is too flat
    for (G in list(10, "auto")) for (b in c("poly3", "stukel", "sym")) {
      e <- edge_external_ref(y, p, G, b)
      a <- expect_silent(edge.gof(y, predicted_probs = p, G = G, basis = b, external = TRUE))
      expect_equal(a$Test_Statistic, e$statistic, tolerance = 1e-10)
      expect_identical(a$df, e$df)
      expect_equal(a$p_value, e$p_value, tolerance = 1e-10)
      expect_identical(a$Method, "external")
      expect_identical(a$Test, "EDGE")
    }
  }
})

test_that("external mode has d + 1 degrees of freedom and a one-sided Stukel basis loses a column", {
  set.seed(3)
  p <- plogis(rnorm(800, -0.3, 1)); y <- rbinom(800, 1, p)
  expect_identical(def.gof(y, p, basis = "poly3", external = TRUE)$df, 4L)
  expect_identical(def.gof(y, p, basis = "poly2", external = TRUE)$df, 3L)
  expect_identical(def.gof(y, p, basis = "stukel", external = TRUE)$df, 4L)
  expect_identical(def.gof(y, p, basis = "sym", external = TRUE)$df, 2L)
  lo <- plogis(rnorm(800, -3, 0.4)); ylo <- rbinom(800, 1, lo)       # every group mean risk below 0.5
  a <- def.gof(ylo, lo, basis = "stukel", external = TRUE)
  expect_identical(a$df, 3L)
  expect_equal(a$Test_Statistic, edge_external_ref(ylo, lo, 10, "stukel")$statistic, tolerance = 1e-10)
})

test_that("external mode on a glm uses its y and fitted probabilities as frozen predictions", {
  set.seed(8)
  x <- runif(500, -3, 3); y <- rbinom(500, 1, plogis(0.6 * x))
  fit <- glm(y ~ x, family = binomial())
  expect_identical(edge.gof(fit, external = TRUE),
                   edge.gof(as.numeric(fit$y), predicted_probs = fitted(fit), external = TRUE))
})

test_that("external mode: no conservative warning, and clear errors for what it does not support", {
  set.seed(5)
  p <- plogis(rnorm(400)); y <- rbinom(400, 1, p)
  expect_silent(def.gof(y, p, external = TRUE))
  expect_warning(def.gof(y, p), "conservative")                    # the default mode still warns
  expect_error(def.gof(y, p, weights = "score", external = TRUE), "weights = 'unit' only")
  expect_error(def.gof(y, p, basis = "ensemble", external = TRUE), "not 'ensemble'")
  expect_error(def.gof(y, external = TRUE), "predicted_probs")
  expect_error(def.gof(y, p, external = NA), "TRUE or FALSE")
  expect_warning(def.gof(y, p, X = cbind(1, qlogis(p)), external = TRUE), "ignored")
})

test_that("the default mode is unchanged: values saved before external mode was added", {
  make_fit <- function(seed = 1, n = 600, link = "logit") {
    set.seed(seed)
    x <- runif(n, -3, 3); z <- rnorm(n)
    eta <- 0.6 * x + 0.3 * z
    p <- if (link == "cloglog") 1 - exp(-exp(eta)) else 1 / (1 + exp(-eta))
    glm(rbinom(n, 1, p) ~ x + z, family = binomial())
  }
  data("gof_demo", package = "ebrahim.gof")
  fits <- list(
    a = make_fit(1), b = make_fit(2, 437), c = make_fit(7, 1500, "cloglog"), d = make_fit(5, 250),
    wrong = glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial()),
    right = glm(outcome ~ poly(age, 2) + bmi + sex + treatment, data = gof_demo, family = binomial()))
  saved <- list(                                                   # statistic, df, p-value (2.9.0 before the edit)
    "def|a|poly3|10|unit|satterthwaite" = c(1.883451544571022, 2.0463037055988722, 0.39132446051372227),
    "def|a|stukel|20|unit|imhof" = c(0.89918826303508415, 1.9362271316866078, 0.57397592252967833),
    "def|b|poly2|7|unit|satterthwaite" = c(0.08585282999794229, 1.0596508982459689, 0.78629393922040292),
    "def|c|sym|10|score|satterthwaite" = c(6.8652515692383087, 1, 0.008788787436423964),
    "yp|wrong|poly3|20|unit|satterthwaite" = c(12.084462094679509, 2.0268885863004025, 0.002216057067383002),
    "naive|d|stukel|10|unit|satterthwaite" = c(3.7729907684211774, 3, 0.28704342684822476),
    "edge|right|sym|10|unit|satterthwaite" = c(0.12493535309096107, 1, 0.58038362951152556),
    "def|wrong|poly3|10|score|satterthwaite" = c(10.473662508945077, 3, 0.014940628000779695),
    "auto|a|poly3" = c(0.50477806116029589, 2.0103957066865861, 0.77538029664447239),
    "auto|d|sym" = c(0.17281070111127222, 1, 0.23171210815173343))
  for (k in names(saved)) {
    s <- strsplit(k, "|", fixed = TRUE)[[1]]
    f <- fits[[s[2]]]
    if (s[1] == "auto") {
      r <- edge.gof(f, G = "auto", basis = s[3])
    } else {
      if (s[6] == "imhof" && !requireNamespace("CompQuadForm", quietly = TRUE)) next
      G <- as.numeric(s[4])
      r <- suppressWarnings(switch(s[1],
        def   = def.gof(f, G = G, basis = s[3], weights = s[5], method = s[6]),
        edge  = edge.gof(f, G = G, basis = s[3], weights = s[5], method = s[6]),
        yp    = def.gof(as.numeric(f$y), fitted(f), X = model.matrix(f), G = G, basis = s[3],
                        weights = s[5], method = s[6]),
        naive = def.gof(as.numeric(f$y), fitted(f), G = G, basis = s[3], weights = s[5], method = s[6])))
    }
    expect_equal(c(r$Test_Statistic, r$df, r$p_value), saved[[k]], tolerance = 1e-12, info = k)
  }
})

test_that("G = 'auto' is max(10, ceiling(n / 25))", {
  expect_identical(.def_auto_G(200), 10)
  expect_identical(.def_auto_G(250), 10)
  expect_identical(.def_auto_G(251), 11)                             # round() gave 10
  expect_identical(.def_auto_G(610), 25)                             # round() gave 24
  expect_identical(.def_auto_G(600), 24)
  set.seed(12)
  x <- runif(610, -3, 3); y <- rbinom(610, 1, plogis(0.6 * x))
  fit <- glm(y ~ x, family = binomial())
  expect_identical(def.gof(fit, G = "auto"), def.gof(fit, G = 25))
  expect_identical(edge.gof(fit, G = "auto", external = TRUE), edge.gof(fit, G = 25, external = TRUE))
  expect_identical(def.ensemble.gof(fit, G = "auto"), def.ensemble.gof(fit, G = 25))
})
