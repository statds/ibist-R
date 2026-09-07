test_that("prop.paired.ci returns stable Wald-family intervals", {
  result <- prop.paired.ci(b = 8, c = 25, n = 180, method = "wald")
  corrected <- prop.paired.ci(b = 8, c = 25, n = 180, method = "waldcc")
  agresti_min <- prop.paired.ci(b = 8, c = 25, n = 180, method = "agresti-min")

  expect_s3_class(result, "ci")
  expect_equal(result$estimate, c("proportion difference" = -17 / 180))
  expect_equal(as.numeric(result$conf.int), c(-0.155455, -0.03343391),
               tolerance = 1e-5)
  expect_equal(as.numeric(corrected$conf.int), c(-0.1610106, -0.02787836),
               tolerance = 1e-5)
  expect_equal(as.numeric(agresti_min$conf.int), c(-0.1547155, -0.03209658),
               tolerance = 1e-5)
})

test_that("prop.paired.ci score interval matches book exercise values", {
  ci1 <- prop.paired.ci(b = 160, c = 3, n = 320)$conf.int
  ci2 <- prop.paired.ci(b = 154, c = 3, n = 320)$conf.int
  ci3 <- prop.paired.ci(b = 43, c = 23, n = 320)$conf.int
  ci4 <- prop.paired.ci(b = 23, c = 96, n = 320)$conf.int

  expect_equal(as.numeric(ci1), c(0.4334417, 0.5466968), tolerance = 1e-4)
  expect_equal(as.numeric(ci2), c(0.4149951, 0.5281011), tolerance = 1e-4)
  expect_equal(as.numeric(ci3), c(0.01306769, 0.11320331), tolerance = 1e-4)
  expect_equal(as.numeric(ci4), c(-0.2902044, -0.1657910), tolerance = 1e-4)
})

test_that("prop.paired.ci validates inputs", {
  expect_error(prop.paired.ci(1, 2, 2), "must not exceed")
  expect_error(prop.paired.ci(1.5, 2, 5), "integer")
  expect_error(prop.paired.ci(1, -2, 5), "non-negative")
  expect_error(prop.paired.ci(1, 2, 5, conf.level = 1), "conf.level")
})

test_that("prop.paired.ci supports Wang's exact interval", {
  result <- prop.paired.ci(
    b = 3, c = 0, n = 4, method = "wang", precision = 0.0001
  )

  expect_s3_class(result, "ci")
  expect_equal(as.numeric(result$conf.int), c(-0.2494, 0.9937),
               tolerance = 1e-4)
  expect_match(result$method, "Wang exact")
})
