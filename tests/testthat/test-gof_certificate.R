## gof_certificate(): the recession slope must reproduce the published certificates, agree with
## the direction deposited as the worst on the sphere, and stay homogeneous of degree zero after
## the unit-norm rescaling.

## the map directions of the paper's certificate table: a value per bin of the first axis,
## repeated across the second (row-major, the order .dg_score() reads)
vec_rows <- function(val, K = 6L) as.numeric(sapply(seq_len(K), function(i) rep(val[i], K)))

test_that("the published certificates are reproduced exactly", {
  step   <- c(-1, -1, -1, 1, 1, 1)
  stripe <- c(1.2, -0.4, -1.1, -1.1, -0.4, 1.2)
  g1 <- gof_certificate(direction = vec_rows(step))
  g2 <- gof_certificate(direction = vec_rows(stripe))
  expect_equal(round(g1$certificate, 3), 2.200)
  expect_equal(round(g1$mirror, 3), 1.075)
  expect_equal(round(g2$certificate, 3), 1.729)
  expect_equal(round(g2$mirror, 3), 1.644)
  expect_match(g1$verdict, "detectable")
})

test_that("the certificate is reported per unit map norm, so scaling a direction cannot move it", {
  v <- vec_rows(c(-1, -1, -1, 1, 1, 1))
  expect_equal(gof_certificate(direction = v)$certificate,
               gof_certificate(direction = 1000 * v)$certificate)
  ## and a matrix argument means the same thing as its row-major vector
  expect_equal(gof_certificate(direction = matrix(v, 6, 6, byrow = TRUE))$certificate,
               gof_certificate(direction = v)$certificate)
})

test_that("a departure that cancels inside every cell has no direction at all", {
  skip_on_cran()
  set.seed(3)
  n <- 150
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.2 + 0.7 * x1 - 0.4 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  g <- gof_certificate(fit, pi_star = as.numeric(fitted(fit)))
  expect_true(g$flat)
  expect_match(g$verdict, "cancels")
})

test_that("naming the departure on the two scales gives the same certificate", {
  skip_on_cran()
  set.seed(4)
  n <- 150
  x1 <- runif(n, -3, 3); x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.2 + 0.7 * x1 - 0.4 * x2))
  fit <- glm(y ~ x1 + x2, family = binomial())
  sh <- 0.8 * (x1^2 - mean(x1^2))
  a <- gof_certificate(fit, eta_shift = sh)
  b <- gof_certificate(fit, pi_star = plogis(qlogis(fitted(fit)) + sh))
  expect_identical(a$certificate, b$certificate)
  expect_gt(a$certificate, 0)              # an omitted quadratic is detectable
  expect_equal(a$axes, c("x1", "x2"))
})

test_that("the direction is taken on the map deepgof1() reads", {
  skip_on_cran()
  set.seed(6)
  n <- 300
  x1 <- runif(n, -3, 3); x2 <- rnorm(n); x3 <- rnorm(n)
  y <- rbinom(n, 1, plogis(0.2 + 0.7 * x1 - 0.4 * x2 + 0.3 * x3))
  ## with pi_star = y the direction is the observed residual map itself (no ties, so no random
  ## tie-break inside deepgof1())
  for (fit in list(glm(y ~ x1 + x2 + x3, family = binomial()),
                   glm(y ~ x1 + I(x1^2) + x2, family = binomial()),
                   glm(y ~ x1, family = binomial()))) {
    g <- gof_certificate(fit, pi_star = fit$y)
    d <- deepgof1(fit, B = 1L)
    expect_identical(g$axes, d$axes)
    expect_equal(g$direction, d$map / sqrt(sum(d$map^2)), tolerance = 1e-12)
  }
})

test_that("bad calls are refused by name", {
  set.seed(5)
  y <- rbinom(40, 1, 0.5); x1 <- rnorm(40); x2 <- rnorm(40)
  fit <- glm(y ~ x1 + x2, family = binomial())
  expect_error(gof_certificate(fit), "name the departure")
  expect_error(gof_certificate(direction = 1:10), "36 entries")
  expect_error(gof_certificate(fit, pi_star = rep(2, 40)), "\\[0, 1\\]")
  expect_error(gof_certificate(fit, eta_shift = 1:3), "one entry per observation")
  expect_error(gof_certificate(direction = rep(1, 36), K = 8L), "valid only at K")
})
