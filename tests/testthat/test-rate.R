test_that("rate.1s.ci returns stable exact intervals", {
  result <- rate.1s.ci(5, T = 10)

  expect_s3_class(result, "ci")
  expect_output(print(result), "95% confidence interval")
  expect_equal(unname(result$estimate), 0.5)
  expect_equal(unname(result$conf.int), c(0.1623486, 1.1668332),
               tolerance = 1e-7)
})

test_that("rate.1s.ci methods return stable intervals", {
  methods <- c("score", "wh", "wald", "log")
  expected <- list(
    score = c(0.2135701, 1.1705758),
    wh = c(0.1611348, 1.1668164),
    wald = c(0.06173873, 0.93826127),
    log = c(0.2081139, 1.2012652)
  )

  for (method in methods) {
    result <- rate.1s.ci(5, T = 10, method = method, correct = FALSE)
    expect_equal(unname(result$conf.int), expected[[method]],
                 tolerance = 1e-7)
  }
})

test_that("rate.1s.ci Wilson-Hilferty continuity correction uses midpoint", {
  result <- rate.1s.ci(5, T = 10, method = "wh", correct = TRUE)

  expect_equal(unname(result$conf.int), c(0.1896388, 1.0959559),
               tolerance = 1e-7)
})

test_that("rate.2s.ci returns log-Wald rate ratio intervals", {
  result <- rate.2s.ci(x = c(151, 55), T = c(57518.1, 74573.5))

  expect_s3_class(result, "ci")
  expect_equal(unname(result$estimate["rate ratio"]), 3.559543,
               tolerance = 1e-6)
  expect_equal(as.numeric(result$conf.int), c(2.614178, 4.846780),
               tolerance = 1e-6)
})

test_that("rate.2s.ci returns transformed binomial intervals", {
  score <- rate.2s.ci(x = c(151, 55), T = c(57518.1, 74573.5),
                      method = "score")
  exact <- rate.2s.ci(x = c(9, 12), T = c(1817.6, 7496.3),
                      method = "exact")
  poisson <- poisson.test(x = c(9, 12), T = c(1817.6, 7496.3))

  expect_equal(as.numeric(score$conf.int), c(2.617350, 4.840907),
               tolerance = 1e-6)
  expect_equal(as.numeric(exact$conf.int), as.numeric(poisson$conf.int),
               tolerance = 1e-12)
})

test_that("rate.test handles one- and two-sample tests", {
  one_sample <- rate.test(x = 411, T = 25800, r = 0.0119,
                          correct = FALSE)
  two_sample <- rate.test(x = c(12, 5), T = c(100, 80),
                          correct = FALSE)

  expect_s3_class(one_sample, "htest")
  expect_s3_class(two_sample, "htest")
  expect_equal(unname(one_sample$estimate), 411 / 25800)
  expect_equal(two_sample$p.value, 0.2122691, tolerance = 1e-6)
})

test_that("rate.test reports a log-Wald rate ratio interval for two samples", {
  result <- rate.test(x = c(151, 55), T = c(57518.1, 74573.5))

  expect_equal(unname(result$estimate["rate ratio"]), 3.559543,
               tolerance = 1e-6)
  expect_equal(as.numeric(result$conf.int), c(2.614178, 4.846780),
               tolerance = 1e-6)
  expect_equal(unname(result$null.value), 1)
  expect_named(result$null.value, "rate ratio")
  expect_match(result$method, "rate ratio by log-Wald")
})

test_that("rate.test uses score interval for one-sample rate", {
  greater <- rate.test(x = 411, T = 25800, r = 0.0119,
                       alternative = "greater")
  less <- rate.test(x = 411, T = 25800, r = 0.0119,
                    alternative = "less")
  score_ci <- rate.1s.ci(411, T = 25800, method = "score")$conf.int

  expect_equal(as.numeric(greater$conf.int), as.numeric(score_ci))
  expect_equal(as.numeric(less$conf.int), as.numeric(score_ci))
  expect_equal(as.numeric(greater$conf.int), c(0.01444433596, 0.01756689433),
               tolerance = 1e-7)
  expect_equal(attr(greater$conf.int, "conf.level"), 0.95)
})

test_that("rate.test handles one-sided alternatives", {
  greater <- rate.test(x = 411, T = 25800, r = 0.0119,
                       alternative = "greater", correct = FALSE)
  less <- rate.test(x = 411, T = 25800, r = 0.0119,
                    alternative = "less", correct = FALSE)

  expect_equal(greater$p.value, 1.47588e-09, tolerance = 1e-6)
  expect_equal(greater$p.value + less$p.value, 1, tolerance = 1e-12)
  expect_equal(greater$alternative, "greater")
  expect_equal(less$alternative, "less")
})

test_that("rate functions validate inputs", {
  expect_error(rate.1s.ci(-1), "non-negative integer")
  expect_error(rate.1s.ci(1, T = 0), "positive")
  expect_error(rate.2s.ci(c(1, 2, 3)), "length-2")
  expect_error(rate.2s.ci(c(1, 2), T = c(1, 0)), "positive exposures")
  expect_error(rate.2s.ci(c(0, 0)), "at least one event")
  expect_error(rate.2s.ci(c(0, 2), method = "log"),
               "requires positive event counts")
  expect_equal(as.numeric(rate.2s.ci(c(0, 2), method = "score")$conf.int[1]),
               0)
  expect_error(rate.test(x = c(1, 2), T = 1), "same length")
  expect_error(rate.test(x = c(0, 0), T = c(1, 1)),
               "at least one event")
  expect_error(rate.test(x = c(0, 2), T = c(1, 1)),
               "requires positive event counts")
})
