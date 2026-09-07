#' Confidence Intervals for Paired Proportion Difference
#'
#' Computes confidence intervals for the difference of two paired proportions.
#'
#' @param b Number of discordant pairs favorable to group 1.
#' @param c Number of discordant pairs favorable to group 2.
#' @param n Total number of matched pairs.
#' @param conf.level Confidence level for the interval.
#' @param method Method for confidence interval. One of \code{"score"},
#'   \code{"wald"}, \code{"waldcc"}, \code{"agresti-min"}, or
#'   \code{"wang"}.
#' @param ... Additional arguments for the selected method. For
#'   \code{"score"}, these include \code{tol}; for \code{"wang"}, they
#'   include \code{precision}, \code{grid.one}, and \code{grid.two}.
#'
#' @section Method-specific arguments:
#' For \code{method = "score"}, \code{tol} is the numerical tolerance used
#' when solving for the confidence limits; its default is \code{1e-8}. For
#' \code{method = "wang"}, \code{precision} controls the numerical precision
#' of the confidence limits (default \code{1e-5}), while \code{grid.one} and
#' \code{grid.two} control the sizes of the two nuisance-parameter search
#' grids (defaults \code{30} and \code{20}, respectively). Increasing these
#' values can improve numerical accuracy at the cost of longer computation.
#'
#' @details
#' For a matched-pair table with discordant counts \eqn{b} and \eqn{c}, the
#' estimated paired difference is \eqn{(b - c) / n}. The \code{"wald"} and
#' \code{"waldcc"} methods use the paired Wald standard error, with
#' \code{"waldcc"} adding the continuity-correction term \eqn{1 / n}. The
#' \code{"agresti-min"} method applies the Wald calculation to
#' \eqn{b + 1/2}, \eqn{c + 1/2}, and \eqn{n + 2}. The \code{"score"} method
#' inverts Tango's score statistic by numerically maximizing the multinomial
#' likelihood subject to the candidate difference. The \code{"wang"} method
#' computes Wang's exact interval based on an inductive ordering of the paired
#' sample space. This method can be time-consuming when the sample size is
#' moderate because it searches over nuisance parameters and refines the
#' confidence limits numerically.
#'
#' @return An object of class \code{"ci"} containing the estimate and
#'   confidence limits.
#'
#' @examples
#' prop.paired.ci(b = 8, c = 25, n = 180)
#' prop.paired.ci(b = 8, c = 25, n = 180, method = "waldcc")
#' prop.paired.ci(b = 8, c = 25, n = 180, method = "agresti-min")
#' prop.paired.ci(b = 3, c = 0, n = 4, method = "wang",
#'                 precision = 0.0001)
#'
#' @references
#' Tango, T. (1998). Equivalence test and confidence interval for the difference
#' in proportions for the paired-sample design. \emph{Statistics in Medicine},
#' 17, 891--908.
#'
#' Agresti, A., and Min, Y. (2005). Simple improved confidence intervals for
#' comparing matched proportions. \emph{Statistics in Medicine}, 24, 729--740.
#'
#' Wang, W. (2012). An inductive order construction for the difference of two
#' dependent proportions. \emph{Statistics & Probability Letters}, 82,
#' 1623--1628.
#'
#' @importFrom stats optimize qnorm uniroot
#' @export
prop.paired.ci <- function(
  b,
  c,
  n,
  conf.level = 0.95,
  method = c("score", "wald", "waldcc", "agresti-min", "wang"),
  ...
) {
  counts <- validate_paired_counts(b, c, n)
  b <- counts[["b"]]
  c <- counts[["c"]]
  n <- counts[["n"]]

  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
        !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("'conf.level' must be a single number between 0 and 1.")
  }

  method <- match.arg(method)
  estimate <- (b - c) / n
  ci <- switch(
    method,
    score = paired_ci_score(b, c, n, conf.level, ...),
    wald = paired_ci_wald(b, c, n, conf.level, correct = FALSE),
    waldcc = paired_ci_wald(b, c, n, conf.level, correct = TRUE),
    "agresti-min" = paired_ci_agresti_min(b, c, n, conf.level),
    wang = paired_ci_wang(b, c, n, conf.level, ...)
  )

  structure(
    list(
      statistic = NULL,
      parameter = NULL,
      p.value = NULL,
      conf.int = structure(ci, conf.level = conf.level),
      estimate = c("proportion difference" = estimate),
      null.value = c("proportion difference" = 0),
      alternative = "two.sided",
      method = paste(
        switch(
          method,
          score = "Score",
          wald = "Wald",
          waldcc = "Continuity-corrected Wald",
          "agresti-min" = "Agresti-Min",
          wang = "Wang exact"
        ),
        "CI for paired proportion difference"
      ),
      data.name = paste0("b = ", b, ", c = ", c, ", n = ", n)
    ),
    class = "ci"
  )
}

paired_ci_wang <- function(b, c, n, conf.level, precision = 0.00001,
                           grid.one = 30, grid.two = 20) {
  result <- wang.paired.ci(
    n10 = b,
    t = n - b - c,
    n01 = c,
    conf.level = conf.level,
    precision = precision,
    grid.one = grid.one,
    grid.two = grid.two
  )
  result$ExactCI
}

validate_paired_counts <- function(b, c, n) {
  vals <- c(b = b, c = c, n = n)
  if (any(lengths(list(b, c, n)) != 1L) ||
        any(!is.finite(vals)) || any(vals < 0) || any(vals != floor(vals))) {
    stop("'b', 'c', and 'n' must be single non-negative integer counts.")
  }
  if (b + c > n)
    stop("'b + c' must not exceed 'n'.")
  vals <- as.integer(vals)
  names(vals) <- c("b", "c", "n")
  vals
}

paired_ci_wald <- function(b, c, n, conf.level, correct) {
  alpha <- 1 - conf.level
  z <- stats::qnorm(1 - alpha / 2)
  estimate <- (b - c) / n
  se <- sqrt((b + c) - (b - c)^2 / n) / n
  correction <- if (correct) 1 / n else 0

  truncate_unit(estimate + c(-1, 1) * (z * se + correction))
}

paired_ci_agresti_min <- function(b, c, n, conf.level) {
  alpha <- 1 - conf.level
  z <- stats::qnorm(1 - alpha / 2)
  b_star <- b + 0.5
  c_star <- c + 0.5
  n_star <- n + 2
  estimate <- (b_star - c_star) / n_star
  se <- sqrt(n_star * (b_star + c_star) - (b_star - c_star)^2) /
    (n_star * sqrt(n_star))

  truncate_unit(estimate + c(-1, 1) * z * se)
}

paired_ci_score <- function(b, c, n, conf.level, tol = 1e-8, ...) {
  estimate <- (b - c) / n
  alpha <- 1 - conf.level
  z <- stats::qnorm(1 - alpha / 2)
  statistic <- function(delta) paired_score_stat(b, c, n, delta)

  lower <- if (estimate <= -1) {
    -1
  } else {
    find_paired_score_limit(
      statistic = statistic,
      target = z,
      interval = c(-1 + tol, estimate),
      boundary = -1,
      tol = tol
    )
  }

  upper <- if (estimate >= 1) {
    1
  } else {
    find_paired_score_limit(
      statistic = statistic,
      target = z,
      interval = c(estimate, 1 - tol),
      boundary = 1,
      tol = tol
    )
  }

  c(lower, upper)
}

paired_score_stat <- function(b, c, n, delta) {
  p <- constrained_paired_mle(b, c, n, delta)
  variance <- (p[["p10"]] + p[["p01"]] - delta^2) / n
  if (variance <= 0)
    return(Inf)

  abs((b - c) / n - delta) / sqrt(variance)
}

constrained_paired_mle <- function(b, c, n, delta) {
  a <- n - b - c
  d <- 0
  q <- c(a = a, b = b, c = c, d = d)
  lower_x <- max(0, -delta)
  upper_x <- (1 - delta) / 2

  neg_loglik <- function(x) {
    p10 <- x + delta
    p01 <- x
    concordant <- 1 - p10 - p01
    if (concordant < 0)
      return(Inf)

    if (a == 0) {
      p11 <- 0
      p00 <- concordant
    } else if (d == 0) {
      p11 <- concordant
      p00 <- 0
    } else {
      p11 <- concordant * a / (a + d)
      p00 <- concordant * d / (a + d)
    }

    probs <- c(p11, p10, p01, p00)
    if (any(probs < 0) || any(q > 0 & probs <= 0))
      return(Inf)
    -sum(q[q > 0] * log(probs[q > 0]))
  }

  opt <- stats::optimize(neg_loglik, c(lower_x, upper_x))
  x <- opt$minimum
  p10 <- x + delta
  p01 <- x
  concordant <- 1 - p10 - p01

  c(p10 = p10, p01 = p01, p11 = concordant, p00 = 0)
}

find_paired_score_limit <- function(statistic, target, interval, boundary,
                                    tol) {
  objective <- function(delta) {
    vapply(delta, function(x) statistic(x) - target, numeric(1))
  }
  values <- objective(interval)

  if (any(!is.finite(values)) || prod(values) > 0)
    return(boundary)

  stats::uniroot(objective, interval = interval, tol = tol)$root
}

truncate_unit <- function(x) {
  pmax(-1, pmin(1, x))
}
