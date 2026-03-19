#' Confidence Intervals for a One-Sample Poisson Rate
#'
#' Computes confidence intervals for a Poisson rate parameter
#' \eqn{\lambda = x / T} using several methods.
#'
#' @param x Non-negative integer. Observed number of events.
#' @param T Positive numeric. Exposure time (or total time at risk).
#' @param conf.level Confidence level. Default is 0.95.
#' @param method Method for confidence interval. One of:
#'   \itemize{
#'     \item \code{"exact"}: Exact (Garwood) interval
#'     \item \code{"score"}: Score interval (inversion of score test)
#'     \item \code{"wh"}: Wilson-Hilferty interval
#'     \item \code{"wald"}: Wald interval (normal approximation)
#'     \item \code{"log"}: Log-transformed interval
#'   }
#' @param correct Logical. Apply continuity correction for
#'   \code{"score"} and \code{"wald"} methods. Default is TRUE.
#' @param ... Reserved for future extensions.
#'
#' @details
#' This function provides several confidence intervals for the Poisson
#' rate \eqn{\lambda = x / T}.
#'
#' \strong{Exact (Garwood):}
#' Based on inversion of the Poisson test using chi-square quantiles:
#' \deqn{
#'   \left[
#'     \frac{1}{2T} \chi^2_{2x, \alpha/2},
#'     \frac{1}{2T} \chi^2_{2(x+1), 1-\alpha/2}
#'   \right].
#' }
#'
#' \strong{Score:}
#' Obtained by inverting the score test. Computed from solving quadratic
#' equations for boundaries.
#'
#' \strong{Wilson-Hilferty:}
#' The Wilson–Hilferty (WH) interval is based on a cube-root transformation
#' of a chi-square approximation, which reduces skewness and yields an
#' approximately normal pivot, leading to the cubic form of the limits.
#' The lower limit is
#' \deqn{
#' \frac{x}{T}
#' \left(1 - \frac{1}{9x} + \frac{z_{\alpha/2}}{3\sqrt{x}}\right)^3
#' }
#' and the upper limit is
#' \deqn{
#' \frac{x + 1}{T}
#' \left(1 - \frac{1}{9(x + 1)} +
#' \frac{z_{1-\alpha/2}}{3\sqrt{x + 1}}\right)^3.
#' }
#'
#' \strong{Wald:}
#' For \eqn{x > 0},
#' \deqn{\hat{\lambda} \pm z_{\alpha/2} \sqrt{\hat{\lambda}/T}.}
#' For \eqn{x = 0}, lower limit 0, upper limit \eqn{-\log(\alpha/2)/T}.
#'
#' \strong{Log:}
#' For \eqn{x > 0}, applies normal approximation to \eqn{\log(\lambda)}
#' and transforms back. For \eqn{x = 0}, lower limit 0, upper limit
#' \eqn{-\log(\alpha/2)/T}.
#'
#' When `correct = TRUE`, continuity correction is applied by replacing
#' `x` with `x - 0.5` in the lower limit and `x + 1 + 0.5` in the upper limit.
#' 
#' @return An object of class \code{"htest"}.
#'
#' @examples
#' rate.1s.ci(5, 10)
#' rate.1s.ci(5, 10, method = "score")
#' rate.1s.ci(0, 10, method = "exact")
#'
#' @importFrom stats qchisq

#' @export
rate.1s.ci <- function(
  x, T = 1.0,
  conf.level = 0.95,
  method = c("exact", "score", "wh", "wald", "log"),
  correct = TRUE,
  ...
) {

  if (length(x) != 1 || x < 0 || x != floor(x)) {
    stop("x must be a non-negative integer.")
  }

  if (length(T) != 1 || T <= 0) {
    stop("T must be positive.")
  }

  if (conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be in (0, 1).")
  }

  method <- match.arg(method)

  ci_methods <- list(
    exact = ci_exact,
    score = ci_score,
    wald  = ci_wald,
    log   = ci_log,
    wh    = ci_wh
  )

  ci <- ci_methods[[method]](x, T, conf.level, correct, ...)

  structure(
    list(
      conf.int = ci,
      estimate = c(rate = x / T),
      conf.level = conf.level,
      method = paste(method, "CI for Poisson rate"),
      data.name = paste0("x = ", x, ", T = ", T)
    ),
    class = "htest"
  )
}

ci_exact <- function(x, T, conf.level, ...) {
  alpha <- 1 - conf.level
  lower <- if (x == 0) 0 else stats::qchisq(alpha / 2, 2 * x) / (2 * T)
  upper <- qchisq(1 - alpha / 2, 2 * (x + 1)) / (2 * T)
  c(lower, upper)
}

ci_score <- function(x, T, conf.level, correct = TRUE, ...) {
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)

  if (correct) {
    ## optional: adjust x → x ± 0.5
    x_lower <- max(0, x - 0.5)
    x_upper <- x + 0.5
  } else {
    x_lower <- x_upper <- x
  }

  lower <- (2 * x_lower + z^2 -
            z * sqrt(z^2 + 4 * x_lower)) / (2 * T)

  upper <- (2 * x_upper + z^2 +
            z * sqrt(z^2 + 4 * x_upper)) / (2 * T)

  c(lower, upper)
}

ci_wh <- function(x, T, conf.level, correct = TRUE, ...) {
  alpha <- 1 - conf.level
  zL <- qnorm(alpha / 2)
  zU <- qnorm(1 - alpha / 2)

  # continuity correction offset
  d <- if (correct) 0.5 else 0

  lower_x <- x - d
  upper_x <- x + 1 + d

  lower <- if (lower_x <= 0) {
    0
  } else {
    lower_x * (1 - 1 / (9 * lower_x) + zL / (3 * sqrt(lower_x)))^3
  }

  upper <- upper_x * (1 - 1 / (9 * upper_x) + zU / (3 * sqrt(upper_x)))^3

  c(lower, upper) / T
}

ci_wald <- function(x, T, conf.level, correct = TRUE, ...) {
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)

  if (x == 0 && !correct) {
    return(c(0, -log(alpha / 2) / T))
  }

  # continuity correction offset
  d <- if (correct) 0.5 else 0

  lower_x <- max(0, x - d)
  upper_x <- x + d

  lower <- max(0, (lower_x - z * sqrt(lower_x)) / T)
  upper <- (upper_x + z * sqrt(upper_x)) / T

  c(lower, upper)
}

ci_log <- function(x, T, conf.level, ...) {
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)

  if (x == 0) return(c(0, -log(alpha / 2) / T))

  rate <- x / T
  se <- 1 / sqrt(x)

  lower <- rate * exp(-z * se)
  upper <- rate * exp(z * se)

  c(lower, upper)
}
