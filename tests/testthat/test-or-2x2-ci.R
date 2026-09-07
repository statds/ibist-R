test_that("or.2x2.ci returns stable Wald intervals", {
  tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)

  wald <- or.2x2.ci(tab)
  adjusted <- or.2x2.ci(tab, method = "adjusted")

  expect_s3_class(wald, "ci")
  expect_equal(wald$estimate, c("odds ratio" = 8.555556), tolerance = 1e-6)
  expect_equal(as.numeric(wald$conf.int), c(0.9904903, 73.9003),
               tolerance = 1e-6)
  expect_equal(as.numeric(adjusted$conf.int), c(0.9827961, 37.7486),
               tolerance = 1e-6)
})

test_that("or.2x2.ci returns Baptista-Pike mid-p interval", {
  tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)

  bp <- or.2x2.ci(tab, method = "bp")

  expect_s3_class(bp, "ci")
  expect_equal(bp$estimate, c("odds ratio" = 8.555556), tolerance = 1e-6)
  expect_equal(as.numeric(bp$conf.int), c(1.3277, 98.8359),
               tolerance = 1e-4)
})

test_that("or.2x2.ci validates inputs", {
  expect_error(or.2x2.ci(matrix(1:6, nrow = 2)), "2 x 2")
  expect_error(or.2x2.ci(matrix(c(1, 2, 3, -1), nrow = 2)), "non-negative")
  expect_error(or.2x2.ci(matrix(c(1, 2, 3, 4.5), nrow = 2)), "integer")
  expect_error(
    or.2x2.ci(matrix(c(1, 2, 3, 4), nrow = 2), conf.level = 1),
    "conf.level"
  )
  expect_error(
    or.2x2.ci(matrix(c(1, 0, 3, 4), nrow = 2)),
    "positive"
  )
})
