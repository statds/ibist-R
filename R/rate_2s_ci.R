#' Confidence Intervals for a Two-Sample Poisson Rate Ratio
#'
#' Computes confidence intervals for the rate ratio
#' \eqn{\rho = \lambda_1 / \lambda_2} from two independent Poisson counts with
#' known exposures.
#'
#' @param x Length-2 vector of non-negative integer event counts.
#' @param T Length-2 positive numeric vector of exposures.
#' @param conf.level Confidence level. Default is 0.95.
#' @param method Method for confidence interval. One of:
#'   \itemize{
#'     \item \code{"log"}: Log-Wald interval
#'     \item \code{"score"}: Transformed Wilson score interval
#'     \item \code{"exact"}: Transformed exact Clopper-Pearson interval
#'   }
#' @param ... Reserved for future extensions.
#'
#' @details
#' For two independent Poisson variables
#' \eqn{X_i \sim \mathrm{Pois}(\lambda_i T_i)}, the point estimate of the
#' rate ratio is
#' \deqn{
#'   \hat\rho = \frac{X_1/T_1}{X_2/T_2}.
#' }
#'
#' The \code{"log"} method uses the large-sample interval
#' \deqn{
#'   \log(\hat\rho) \pm z_{1-\alpha/2}\sqrt{1/X_1 + 1/X_2}
#' }
#' and exponentiates the endpoints. This method requires both counts to be
#' positive.
#'
#' The \code{"score"} and \code{"exact"} methods use the conditional binomial
#' representation. Given \eqn{n = X_1 + X_2},
#' \eqn{X_1 \sim \mathrm{Binom}(n, \pi)}, where
#' \deqn{
#'   \pi = \frac{T_1\rho}{T_1\rho + T_2}.
#' }
#' Confidence limits \eqn{L} and \eqn{U} for \eqn{\pi} are transformed to
#' confidence limits for \eqn{\rho} as
#' \deqn{
#'   \left(
#'   \frac{L T_2}{(1 - L)T_1},
#'   \frac{U T_2}{(1 - U)T_1}
#'   \right).
#' }
#' The \code{"score"} method uses Wilson score limits for \eqn{\pi}; the
#' \code{"exact"} method uses the exact Clopper-Pearson limits returned by
#' \code{\link[stats]{binom.test}}.
#'
#' @return An object of class \code{"ci"} containing the estimate and
#'   confidence limits.
#'
#' @examples
#' rate.2s.ci(c(151, 55), T = c(57518.1, 74573.5))
#' rate.2s.ci(c(151, 55), T = c(57518.1, 74573.5), method = "score")
#' rate.2s.ci(c(9, 12), T = c(1817.6, 7496.3), method = "exact")
#'
#' @export
rate.2s.ci <- function(
  x,
  T = c(1.0, 1.0),
  conf.level = 0.95,
  method = c("log", "score", "exact"),
  ...
) {
  if (length(x) != 2L || any(!is.finite(x)) ||
      any(x < 0) || any(abs(x - round(x)) > 0)) {
    stop("x must be a length-2 vector of non-negative integers.")
  }

  if (length(T) != 2L || any(!is.finite(T)) || any(T <= 0)) {
    stop("T must be a length-2 vector of positive exposures.")
  }

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      conf.level <= 0 || conf.level >= 1) {
    stop("conf.level must be a single number in (0, 1).")
  }

  if (sum(x) == 0) {
    stop("at least one event is required.")
  }

  method <- match.arg(method)
  rate1 <- x[1] / T[1]
  rate2 <- x[2] / T[2]
  rate_ratio <- rate1 / rate2

  ci_methods <- list(
    log = ci_rate_ratio_log,
    score = ci_rate_ratio_score,
    exact = ci_rate_ratio_exact
  )
  ci <- ci_methods[[method]](x, T, conf.level, ...)

  structure(
    list(
      conf.int = structure(ci, conf.level = conf.level),
      estimate = c("rate 1" = rate1, "rate 2" = rate2,
                   "rate ratio" = rate_ratio),
      conf.level = conf.level,
      method = paste(method, "CI for Poisson rate ratio"),
      data.name = paste0(
        "x = c(", x[1], ", ", x[2], "), ",
        "T = c(", T[1], ", ", T[2], ")"
      )
    ),
    class = "ci"
  )
}

ci_rate_ratio_log <- function(x, T, conf.level, ...) {
  if (x[1] == 0 || x[2] == 0) {
    stop("log-Wald confidence interval for the rate ratio requires positive event counts")
  }

  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)
  rate_ratio <- (x[1] / T[1]) / (x[2] / T[2])
  se_log_ratio <- sqrt(1 / x[1] + 1 / x[2])

  exp(log(rate_ratio) + c(-1, 1) * z * se_log_ratio)
}

ci_rate_ratio_score <- function(x, T, conf.level, ...) {
  n <- sum(x)
  pi_ci <- binom_wilson_ci(x[1], n, conf.level)
  transform_pi_to_rate_ratio(pi_ci, T)
}

ci_rate_ratio_exact <- function(x, T, conf.level, ...) {
  pi_ci <- stats::binom.test(x[1], sum(x), conf.level = conf.level)$conf.int
  transform_pi_to_rate_ratio(pi_ci, T)
}

binom_wilson_ci <- function(x, n, conf.level) {
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha / 2)
  p <- x / n
  denominator <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denominator

  c(max(0, center - half_width), min(1, center + half_width))
}

transform_pi_to_rate_ratio <- function(pi_ci, T) {
  lower <- if (pi_ci[1] <= 0) {
    0
  } else {
    pi_ci[1] * T[2] / ((1 - pi_ci[1]) * T[1])
  }

  upper <- if (pi_ci[2] >= 1) {
    Inf
  } else {
    pi_ci[2] * T[2] / ((1 - pi_ci[2]) * T[1])
  }

  c(lower, upper)
}
