#' Confidence Intervals for a 2 x 2 Odds Ratio
#'
#' Computes confidence intervals for the odds ratio in a 2 x 2 table.
#'
#' @param x A 2 x 2 table of non-negative integer counts.
#' @param conf.level Confidence level for the interval.
#' @param method Method for confidence interval. One of:
#'   \code{"wald"}, \code{"adjusted"}, \code{"baptista-pike"}, or \code{"bp"}.
#' @param ... Reserved for future extensions.
#'
#' @details
#' For a table with entries \eqn{n_{11}}, \eqn{n_{12}}, \eqn{n_{21}}, and
#' \eqn{n_{22}}, the sample odds ratio is
#' \eqn{n_{11} n_{22} / (n_{12} n_{21})}.
#'
#' The \code{"wald"} method uses the usual large-sample interval on the
#' log odds-ratio scale. The \code{"adjusted"} method adds 0.5 to every
#' cell before applying the same calculation. The \code{"baptista-pike"}
#' and \code{"bp"} methods invert the conditional mid-p test based on the
#' noncentral hypergeometric distribution.
#'
#' @return An object of class \code{"ci"} containing the estimate and
#'   confidence limits.
#'
#' @references
#' Baptista, J., and Pike, M. C. (1977). Algorithm AS 115: Exact two-sided
#' confidence limits for the odds ratio in a 2 x 2 table.
#' \emph{Journal of the Royal Statistical Society. Series C}, 26(2), 214--220.
#'
#' Fagerland, M. W., Lydersen, S., and Laake, P. (2015). Recommended confidence
#' intervals for two independent binomial proportions.
#' \emph{Statistical Methods in Medical Research}, 24(2), 224--254.
#'
#' @examples
#' tab <- matrix(c(7, 27, 1, 33), nrow = 2, byrow = TRUE)
#' or.2x2.ci(tab)
#' or.2x2.ci(tab, method = "adjusted")
#' or.2x2.ci(tab, method = "bp")
#'
#' @importFrom stats qnorm uniroot
#' @export
or.2x2.ci <- function(
  x,
  conf.level = 0.95,
  method = c("wald", "adjusted", "baptista-pike", "bp"),
  ...
) {
  tab <- validate_2x2_table(x)
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
        !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1.")
  }

  method <- match.arg(method)
  if (method == "bp")
    method <- "baptista-pike"

  ci_methods <- list(
    wald = or_ci_wald,
    adjusted = or_ci_adjusted,
    "baptista-pike" = or_ci_bp_midp
  )

  estimate <- odds_ratio_estimate(tab)
  ci <- ci_methods[[method]](tab, conf.level, ...)

  structure(
    list(
      statistic = NULL,
      parameter = NULL,
      p.value = NULL,
      conf.int = structure(ci, conf.level = conf.level),
      estimate = c("odds ratio" = estimate),
      null.value = c("odds ratio" = 1),
      alternative = "two.sided",
      method = paste(
        switch(
          method,
          wald = "Wald",
          adjusted = "Adjusted Wald",
          "baptista-pike" = "Baptista-Pike mid-p"
        ),
        "CI for odds ratio"
      ),
      data.name = deparse(substitute(x))
    ),
    class = "ci"
  )
}

validate_2x2_table <- function(x) {
  tab <- as.matrix(x)
  if (!is.numeric(tab) || !identical(dim(tab), c(2L, 2L))) {
    stop("'x' must be a 2 x 2 numeric table.")
  }
  if (any(!is.finite(tab)) || any(tab < 0) || any(tab != floor(tab))) {
    stop("'x' must contain non-negative integer counts.")
  }
  storage.mode(tab) <- "integer"
  tab
}

odds_ratio_estimate <- function(tab) {
  a <- tab[1L, 1L]
  b <- tab[1L, 2L]
  c <- tab[2L, 1L]
  d <- tab[2L, 2L]
  (a * d) / (b * c)
}

or_ci_wald <- function(tab, conf.level, ...) {
  if (any(tab == 0)) {
    stop("Wald confidence interval requires all cell counts to be positive.")
  }

  alpha <- 1 - conf.level
  z <- stats::qnorm(1 - alpha / 2)
  estimate <- odds_ratio_estimate(tab)
  se <- sqrt(sum(1 / tab))

  exp(log(estimate) + c(-1, 1) * z * se)
}

or_ci_adjusted <- function(tab, conf.level, ...) {
  tab <- tab + 0.5
  alpha <- 1 - conf.level
  z <- stats::qnorm(1 - alpha / 2)
  estimate <- odds_ratio_estimate(tab)
  se <- sqrt(sum(1 / tab))

  exp(log(estimate) + c(-1, 1) * z * se)
}

or_ci_bp_midp <- function(tab, conf.level, tol = 1e-8, ...) {
  alpha <- 1 - conf.level
  estimate <- odds_ratio_estimate(tab)

  lower <- if (estimate == 0) {
    0
  } else {
    find_bp_limit(tab, alpha, lower = TRUE, estimate = estimate, tol = tol)
  }

  upper <- if (is.infinite(estimate)) {
    Inf
  } else {
    find_bp_limit(tab, alpha, lower = FALSE, estimate = estimate, tol = tol)
  }

  c(lower, upper)
}

find_bp_limit <- function(tab, alpha, lower, estimate, tol) {
  objective <- function(theta) bp_midp_value(tab, theta) - alpha

  if (lower) {
    high <- max(estimate, .Machine$double.eps)
    low <- high
    while (objective(low) > 0 && low > .Machine$double.xmin * 10) {
      low <- low / 10
    }
    if (objective(low) > 0)
      return(0)
    interval <- c(low, high)
  } else {
    low <- max(estimate, .Machine$double.eps)
    high <- low
    while (objective(high) > 0 && high < .Machine$double.xmax / 10) {
      high <- high * 10
    }
    if (objective(high) > 0)
      return(Inf)
    interval <- c(low, high)
  }

  stats::uniroot(objective, interval = interval, tol = tol)$root
}

bp_midp_value <- function(tab, theta) {
  x_obs <- tab[1L, 1L]
  row1 <- sum(tab[1L, ])
  row2 <- sum(tab[2L, ])
  col1 <- sum(tab[, 1L])

  support <- seq.int(max(0L, col1 - row2), min(row1, col1))
  prob <- noncentral_hypergeom_prob(support, row1, row2, col1, theta)
  prob_obs <- noncentral_hypergeom_prob(x_obs, row1, row2, col1, theta)

  sum(prob[prob < prob_obs]) + 0.5 * sum(prob[prob == prob_obs])
}

noncentral_hypergeom_prob <- function(x, row1, row2, col1, theta) {
  support <- seq.int(max(0L, col1 - row2), min(row1, col1))

  log_weights <- lchoose(row1, support) +
    lchoose(row2, col1 - support) +
    support * log(theta)
  log_weights <- log_weights - max(log_weights)
  weights <- exp(log_weights)
  probs <- weights / sum(weights)

  probs[match(x, support)]
}
