#' Test of Poisson Rates Using Normal Approximation
#'
#' Performs a large-sample (normal) test for one or two Poisson rates with known
#' exposures. The test is carried out on the observed count scale and then
#' translated to rates for reporting.
#'
#' @param x a vector of event counts. A single value specifies a one-sample
#'   test; a vector of length two specifies a two-sample comparison.
#' @param exposure a vector of exposures corresponding to \code{x}
#'   (e.g., person-time).
#'   Must have the same length as \code{x}.
#' @param r a positive number specifying the null rate per unit exposure,
#'   \eqn{\lambda_0}, for a one-sample test. Default to 1.0. Ignored for two-sample tests.
#' @param alternative a character string specifying the alternative hypothesis,
#'   one of \code{"two.sided"}, \code{"less"}, or \code{"greater"}.
#' @param conf.level confidence level for the confidence interval. For
#'   one-sample tests, the interval is the score confidence interval for the
#'   rate. For two-sample tests, the interval is the log-Wald confidence
#'   interval for the rate ratio.
#' @param correct logical; if \code{TRUE}, applies a continuity correction on
#'   the count scale. For one-sample tests, a \eqn{\pm 0.5} correction is applied
#'   to \eqn{x} in the direction determined by \code{alternative}. For
#'   \code{"two.sided"}, the correction uses the sign of \eqn{x-\mu_0}.
#'
#' @details
#' One-sample test:
#' Assume \eqn{X \sim \mathrm{Pois}(\mu)} with \eqn{\mu = \lambda T}. The null
#' hypothesis is \eqn{H_0:\lambda=\lambda_0}, where \eqn{\lambda_0=r}. Let
#' \eqn{\mu_0=\lambda_0 T}. The normal approximation gives
#' \deqn{Z = \frac{(x - c) - \mu_0}{\sqrt{\mu_0}} \approx N(0,1),}
#' where \eqn{c} is 0 if \code{correct = FALSE} and is a continuity correction on
#' the count scale if \code{correct = TRUE}. Two-sided p-values use
#' \eqn{2\{1-\Phi(|Z|)\}}.
#'
#' Two-sample test:
#' Assume independent \eqn{X_1 \sim \mathrm{Pois}(\lambda_1 T_1)} and
#' \eqn{X_2 \sim \mathrm{Pois}(\lambda_2 T_2)} with known exposures
#' \eqn{T_1} and \eqn{T_2}. The null hypothesis is
#' \eqn{H_0:\lambda_1=\lambda_2}.
#'
#' For hypothesis testing, the function uses the exact conditional
#' representation under \eqn{H_0}:
#' \deqn{
#' X_1 \mid (X_1+X_2=n)
#' \sim
#' \mathrm{Binom}\!\left(n,\; \frac{T_1}{T_1+T_2}\right),
#' }
#' and applies the same large-sample normal approximation as
#' \code{\link[stats]{prop.test}} (with optional Yates continuity correction)
#' to obtain the test statistic and p-value.
#'
#' For estimation and confidence intervals, the inference target is the
#' rate ratio,
#' \eqn{\rho = \lambda_1 / \lambda_2}.
#' The point estimate is
#' \eqn{\hat\rho = (X_1/T_1) / (X_2/T_2)}, and the confidence interval is
#' constructed on the log scale using estimated standard error
#' \deqn{
#' \sqrt{1/X_1 + 1/X_2}.
#' }
#' The continuity correction affects the hypothesis test but not the
#' confidence interval for \eqn{\rho}.
#'
#' Confidence intervals:
#' For one-sample tests, the confidence interval is the score interval for
#' \eqn{\lambda}, as computed by \code{\link{rate.1s.ci}} with
#' \code{method = "score"}.
#' For two-sample tests, the confidence interval is for the rate ratio
#' \eqn{\lambda_1/\lambda_2}. The log-Wald interval requires both event counts
#' to be positive.
#'
#' @return
#' An object of class \code{"htest"} containing:
#' \item{statistic}{the standardized normal statistic \code{z} for one-sample
#'   tests or the chi-square statistic from \code{\link[stats]{prop.test}} for
#'   two-sample tests.}
#' \item{parameter}{degrees of freedom (\code{df = 1}).}
#' \item{p.value}{the p-value.}
#' \item{conf.int}{a confidence interval for the rate in one-sample tests or
#'   for the rate ratio, \eqn{\lambda_1/\lambda_2}, in two-sample tests.}
#' \item{estimate}{estimated rate (one-sample) or estimated rates and rate ratio
#'   (two-sample).}
#' \item{null.value}{the null rate (one-sample) or the null rate ratio
#'   \eqn{\lambda_1/\lambda_2=1} (two-sample).}
#' \item{alternative}{the alternative hypothesis.}
#' \item{method}{a character string describing the test.}
#' \item{data.name}{a character string describing the data.}
#'
#' @seealso
#' \code{\link{rate.2s.ci}}, \code{\link[stats]{prop.test}},
#' \code{\link[stats]{poisson.test}}
#'
#' @examples
#' ## One-sample test: compare observed rate to a reference unit rate
#' rate.test(x = 411, exposure = 25800, r = 0.0119)
#' rate.test(x = 411, exposure = 25800, r = 0.0119, correct = FALSE)
#'
#' ## Two-sample test: compare two Poisson rates
#' rate.test(x = c(12, 5), exposure = c(100, 80))
#' rate.test(x = c(12, 5), exposure = c(100, 80), correct = FALSE)
#'
#' ## One-sided alternative
#' rate.test(x = 411, exposure = 25800, r = 0.0119,
#'           alternative = "greater")
#'
#' @export
rate.test <- function(x, exposure = 1.0, r = 1.0,
                      alternative = c("two.sided", "less", "greater"),
                      conf.level = 0.95,
                      correct = TRUE)
{
  alternative <- match.arg(alternative)

  if (length(x) != length(exposure)) {
    stop("'x' and 'exposure' must have the same length")
  }
  if (!length(x) %in% c(1L, 2L)) stop("'x' must have length 1 or 2")
  if (any(!is.finite(x)) || any(!is.finite(exposure))) {
    stop("'x' and 'exposure' must be finite")
  }
  if (any(x < 0) || any(abs(x - round(x)) > 0)) stop("'x' must be nonnegative integers")
  if (any(exposure <= 0)) stop("'exposure' must be positive")
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number in (0, 1)")
  }

  if (length(x) == 1L) {
      ## --- One-sample: H0: lambda = r (unit rate)
      if (is.null(r) || length(r) != 1L || !is.finite(r) || r <= 0) {
          stop("a single positive null unit rate 'r' must be specified")
      }
      
      mu0 <- r * exposure
      
      ## CC on count scale
      cc <- 0
      if (isTRUE(correct)) {
          cc <- switch(alternative,
                       "greater"   = 0.5,
                       "less"      = -0.5,
                       "two.sided" = if (x == mu0) 0 else 0.5 * sign(x - mu0)
                       )
      }
      
      ## Standardize observed count against null mean (score-style)
      z <- (x - cc - mu0) / sqrt(mu0)
      
      p.value <- switch(alternative,
                        "two.sided" = 2 * pnorm(-abs(z)),
                        "less"      = pnorm(z),
                        "greater"   = pnorm(z, lower.tail = FALSE)
                        )
      
      conf.int <- rate.1s.ci(
          x = x,
          exposure = exposure,
          conf.level = conf.level,
          method = "score",
          correct = correct
      )$conf.int
      
      estimate <- c(rate = x / exposure)
      null.value <- c(rate = r)
      
      ## Ingredients (shared names)
      statistic <- c(z = as.numeric(z))
      parameter <- c(df = 1)
      method <- "Normal approximation test for a Poisson rate"
      data.name <- paste0("x = ", x, ", exposure = ", exposure)
      
  } else {
      
      ## --- Two-sample: H0: lambda1 = lambda2
      x1 <- x[1]; x2 <- x[2]
      exposure1 <- exposure[1]; exposure2 <- exposure[2]
      
      if (x1 + x2 == 0) {
          stop("at least one event is required for the 2-sample test")
      }

      ci <- rate.2s.ci(
          x, exposure, conf.level = conf.level, method = "log"
      )
      
      ## Conditional binomial test under H0 via prop.test()
      n  <- x1 + x2
      p0 <- exposure1 / (exposure1 + exposure2)
      
      pt <- stats::prop.test(x = x1, n = n,p = p0, alternative = alternative,
                              conf.level = conf.level, correct = correct)
      
      ## Estimation target: rate ratio
      conf.int <- ci$conf.int
      
      estimate <- ci$estimate
      null.value <- c("rate ratio" = 1)
      
      ## Ingredients (shared names)
      statistic <- pt$statistic
      parameter <- pt$parameter
      p.value   <- pt$p.value
      method <- paste(
          "2-sample test for equality of Poisson rates",
          "(conditional binomial normal approximation);",
          "CI for rate ratio by log-Wald approximation"
      )
      data.name <- paste0(
          "x = c(", x1, ", ", x2, "), ",
        "exposure = c(", exposure1, ", ", exposure2, ")"
      )
  }
  
  ## --- Shared return block (single place)
  rval <- list(
      statistic   = statistic,
      parameter   = parameter,
      p.value     = as.numeric(p.value),
      conf.int    = structure(conf.int, conf.level = conf.level),
      estimate    = estimate,
      null.value  = null.value,
      alternative = alternative,
      method      = method,
      data.name   = data.name
  )
  class(rval) <- "htest"
  rval
}
