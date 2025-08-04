test_that("ef.gof works with binary data", {
  # Set up test data
  set.seed(123)
  n <- 200
  x <- rnorm(n)
  linpred <- 0.5 + 1.2 * x
  prob <- plogis(linpred)
  y <- rbinom(n, 1, prob)
  
  model <- glm(y ~ x, family = binomial())
  predicted_probs <- fitted(model)
  
  # Test with default parameters
  result <- ef.gof(y, predicted_probs)
  
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 3)
  expect_equal(names(result), c("Test", "Test_Statistic", "p_value"))
  expect_equal(result$Test, "Ebrahim-Farrington")
  expect_true(is.numeric(result$Test_Statistic))
  expect_true(is.numeric(result$p_value))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that("ef.gof works with different group numbers", {
  set.seed(456)
  n <- 100
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(x))
  predicted_probs <- plogis(x)
  
  # Test with different G values
  result_4 <- ef.gof(y, predicted_probs, G = 4)
  result_10 <- ef.gof(y, predicted_probs, G = 10)
  result_20 <- ef.gof(y, predicted_probs, G = 20)
  
  expect_equal(result_4$Test, "Ebrahim-Farrington")
  expect_equal(result_10$Test, "Ebrahim-Farrington")
  expect_equal(result_20$Test, "Ebrahim-Farrington")
  
  # All should return valid p-values
  expect_true(all(c(result_4$p_value, result_10$p_value, result_20$p_value) >= 0))
  expect_true(all(c(result_4$p_value, result_10$p_value, result_20$p_value) <= 1))
})

test_that("input validation works correctly", {
  set.seed(789)
  n <- 50
  y <- rbinom(n, 1, 0.5)
  predicted_probs <- runif(n, 0.1, 0.9)
  
  # Test length mismatch
  expect_error(
    ef.gof(y[1:10], predicted_probs),
    "Length mismatch"
  )
  
  # Test invalid predicted probabilities
  expect_error(
    ef.gof(y, c(predicted_probs[1:49], 1.5)),
    "must be between 0 and 1"
  )
  
  expect_error(
    ef.gof(y, c(predicted_probs[1:49], -0.1)),
    "must be between 0 and 1"
  )
  
  # Test invalid G
  expect_error(
    ef.gof(y, predicted_probs, G = 1),
    "must be a single integer >= 2"
  )
  
  expect_error(
    ef.gof(y, predicted_probs, G = n + 10),
    "cannot exceed sample size"
  )
  
  # Test non-binary y without m
  expect_error(
    ef.gof(c(0, 1, 2), c(0.1, 0.5, 0.9)),
    "must contain only 0s and 1s"
  )
})

test_that("grouped data (original Farrington) works", {
  # Simulate grouped data
  set.seed(202)
  n_groups <- 20
  m_trials <- sample(5:15, n_groups, replace = TRUE)
  x <- rnorm(n_groups)
  prob_true <- plogis(0.2 + 0.8 * x)
  y_successes <- rbinom(n_groups, m_trials, prob_true)
  
  # Create data frame and fit model
  data_grouped <- data.frame(
    successes = y_successes,
    trials = m_trials,
    x = x
  )
  
  model_grouped <- glm(
    cbind(successes, trials - successes) ~ x,
    data = data_grouped,
    family = binomial()
  )
  
  predicted_probs <- fitted(model_grouped)
  
  # Test original Farrington (should work with grouped data)
  result <- ef.gof(
    y_successes,
    predicted_probs,
    model = model_grouped,
    m = m_trials,
    G = NULL
  )
  
  expect_equal(result$Test, "Farrington-Original")
  expect_true(is.numeric(result$Test_Statistic))
  expect_true(result$p_value >= 0 && result$p_value <= 1)
})

test_that("error handling for original Farrington without required parameters", {
  set.seed(303)
  y <- rbinom(50, 1, 0.5)
  predicted_probs <- rep(0.5, 50)
  
  # Should error when G is NULL but m is not provided
  expect_error(
    ef.gof(y, predicted_probs, G = NULL),
    "you must provide the 'm' vector"
  )
  
  # Should error when m is provided but model is not
  expect_error(
    ef.gof(y, predicted_probs, m = rep(1, 50), G = NULL),
    "you must provide the fitted 'model' object"
  )
})

test_that("handles edge cases gracefully", {
  # Very small dataset
  set.seed(404)
  y_small <- c(0, 1, 0, 1)
  pred_small <- c(0.2, 0.8, 0.3, 0.7)
  
  result_small <- ef.gof(y_small, pred_small, G = 2)
  expect_true(is.numeric(result_small$p_value))
  
  # Extreme probabilities (but within bounds)
  y_extreme <- c(0, 0, 1, 1)
  pred_extreme <- c(0.001, 0.001, 0.999, 0.999)
  
  result_extreme <- ef.gof(y_extreme, pred_extreme, G = 2)
  expect_true(is.numeric(result_extreme$p_value))
})

test_that("test produces consistent results with same seed", {
  set.seed(505)
  n <- 100
  x <- rnorm(n)
  y1 <- rbinom(n, 1, plogis(x))
  predicted_probs1 <- plogis(x)
  
  set.seed(505)
  y2 <- rbinom(n, 1, plogis(x))
  predicted_probs2 <- plogis(x)
  
  result1 <- ef.gof(y1, predicted_probs1, G = 10)
  result2 <- ef.gof(y2, predicted_probs2, G = 10)
  
  expect_identical(result1, result2)
}) 