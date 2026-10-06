test_that("summary.gof_battery() counts each family from the battery itself", {
  data("gof_demo", package = "ebrahim.gof")
  fit <- glm(outcome ~ age + bmi + sex + treatment, data = gof_demo, family = binomial)
  res <- run.all.gof(fit, tests = c("HL", "Pigeon-Heyse", "EF", "EDGE", "EDGE.G10", "Stukel"),
                     install = "no")
  s <- summary(res)
  expect_s3_class(s, "summary.gof_battery")
  expect_s3_class(s, "data.frame")
  expect_identical(names(s), c("Family", "Tests", "With_p", "Rejected", "Smallest_p", "Smallest_test"))
  expect_identical(sum(s$Tests), nrow(res))
  expect_setequal(s$Family, unique(res$Family))
  for (f in s$Family) {
    r <- res[res$Family == f, ]
    k <- s[s$Family == f, ]
    expect_identical(k$Rejected, sum(!is.na(r$p_value) & r$p_value <= 0.05))
    expect_equal(k$Smallest_p, min(r$p_value, na.rm = TRUE))
    expect_identical(k$Smallest_test, r$Test[which.min(r$p_value)])
  }
  # a stricter level can only lower the counts
  expect_true(all(summary(res, alpha = 0.001)$Rejected <= s$Rejected))
  expect_identical(attr(s, "alpha"), 0.05)
  # families come in the order the battery prints them
  expect_identical(s$Family, intersect(.battery_family_order, s$Family))
  expect_output(print(s), "Goodness-of-fit battery by family")
  expect_error(summary(res, alpha = 2), "'alpha'")
})

test_that("a family whose rows have no p-value is listed without one", {
  b <- data.frame(Test = c("a", "b", "c"), Family = c("Partition", "Partition", "Descriptive"),
                  Statistic = c(1, 2, 0.9), df = c(8, 8, NA), p_value = c(0.2, 0.01, NA),
                  Note = NA_character_, stringsAsFactors = FALSE)
  class(b) <- c("gof_battery", "data.frame")
  s <- summary(b)
  expect_identical(s$Family, c("Partition", "Descriptive"))
  expect_identical(s$With_p, c(2L, 0L))
  expect_identical(s$Rejected, c(1L, 0L))
  expect_true(is.na(s$Smallest_p[2]))
  expect_output(print(s), "Descriptive")
})

test_that("an empty battery gives an empty summary", {
  b <- data.frame(Test = character(0), Family = character(0), Statistic = numeric(0),
                  df = numeric(0), p_value = numeric(0), Note = character(0))
  class(b) <- c("gof_battery", "data.frame")
  s <- summary(b)
  expect_identical(nrow(s), 0L)
  expect_output(print(s), "0 tests")
})
