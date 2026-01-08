#' One-sample proportion power and sample size calculation
#'
#' Computes power, sample size, or detectable effect size for a
#' one-sample binomial test. The interface mirrors
#' \code{stats::power.prop.test()}, but for the one-sample setting.
#'
#' Power can be computed using a normal approximation or exact binomial
#' methods. When \code{exact = TRUE}, the exact test is determined by
#' \code{exact.method}.
#'
#' Exactly one of \code{n}, \code{p0}, \code{p1}, \code{power}, or
#' \code{sig.level} must be \code{NULL}; the missing quantity is solved
#' numerically.
#'
#' @param n Sample size for the single group.
#' @param p0 Null hypothesis proportion.
#' @param p1 Alternative hypothesis proportion.
#' @param power Desired power.
#' @param sig.level Significance level.
#' @param alternative Character string specifying the alternative
#'   hypothesis; one of \code{"two.sided"}, \code{"less"}, or
#'   \code{"greater"}.
#' @param cc Logical; if \code{TRUE}, apply continuity correction in the
#'   normal approximation.
#' @param exact Logical; if \code{TRUE}, use an exact binomial method.
#' @param exact.method Method used for exact binomial power calculation.
#'   \code{"quantile"} uses fixed binomial rejection regions for stable
#'   power and sample-size inversion; \code{"midp"} applies the mid-p
#'   adjustment; \code{"cp"} inverts Clopper--Pearson acceptance regions
#'   to guarantee size not exceeding \code{sig.level}.
#' @param tol Numerical tolerance used in root finding.
#' @param max_n Maximum allowable sample size when solving for \code{n}.
#'
#' @return An object of class \code{"power.htest"} containing the
#'   computed quantity and test specifications.
#'
#' @examples
#' ## Normal approximation (default)
#' power.p1s.test(n = 50, p0 = 0.1, p1 = 0.25)
#'
#' ## Exact binomial power (quantile-based, default exact method)
#' power.p1s.test(n = 50, p0 = 0.1, p1 = 0.25, exact = TRUE)
#'
#' ## Exact mid-p power
#' power.p1s.test(
#'   n = 50, p0 = 0.1, p1 = 0.25,
#'   exact = TRUE, exact.method = "midp"
#' )
#'
#' ## Exact Clopper--Pearson power (guaranteed size control)
#' power.p1s.test(
#'   n = 50, p0 = 0.1, p1 = 0.25,
#'   exact = TRUE, exact.method = "cp"
#' )
#'
#' @export
power.p1s.test <- function(
  n = NULL,
  p0 = NULL,
  p1 = NULL,
  power = NULL,
  sig.level = 0.05,
  alternative = c("two.sided", "less", "greater"),
  cc = FALSE,
  exact = FALSE,
  exact.method = c("quantile", "midp", "cp"),
  tol = .Machine$double.eps^0.5,
  max_n = 1e7
) {
  alternative <- match.arg(alternative)
  exact.method <- match.arg(exact.method)

  ## ---- argument geometry ----
  is_null <- c(is.null(n), is.null(p0), is.null(p1),
               is.null(power), is.null(sig.level))
  if (sum(is_null) != 1L)
    stop("Exactly one of 'n', 'p0', 'p1', 'power', 'sig.level' must be NULL",
         call. = FALSE)

  ## ---- basic validation ----
  if (!is.null(n)) {
    if (!is.finite(n) || n <= 0)
      stop("'n' must be positive", call. = FALSE)
    n <- as.integer(ceiling(n))
  }
  if (!is.null(p0) && !(p0 > 0 && p0 < 1))
    stop("'p0' must be in (0,1)", call. = FALSE)
  if (!is.null(p1) && !(p1 > 0 && p1 < 1))
    stop("'p1' must be in (0,1)", call. = FALSE)
  if (!is.null(sig.level) && !(sig.level > 0 && sig.level < 1))
    stop("'sig.level' must be in (0,1)", call. = FALSE)
  if (!is.null(power) && !(power > 0 && power < 1))
    stop("'power' must be in (0,1)", call. = FALSE)

  if (!exact && exact.method != "quantile") {
    warning(
      "'exact.method' ignored when exact = FALSE",
      call. = FALSE
    )
  }

  ## ---- exact quantile rejection region cache ----
  rr_cache <- new.env(parent = emptyenv())

  get_rr <- function(n, p0, alpha) {
    key <- paste(n, p0, alpha, alternative, sep = "|")
    if (exists(key, envir = rr_cache, inherits = FALSE))
      return(get(key, envir = rr_cache, inherits = FALSE))

    if (alternative == "greater") {
      rr <- list(type = "greater",
                 k = qbinom(1 - alpha, n, p0))
    } else if (alternative == "less") {
      rr <- list(type = "less",
                 k = qbinom(alpha, n, p0))
    } else {
      rr <- list(
        type = "two.sided",
        k_lo = qbinom(alpha / 2, n, p0),
        k_hi = qbinom(1 - alpha / 2, n, p0)
      )
    }

    assign(key, rr, envir = rr_cache)
    rr
  }

  ## ---- power bodies ----

  approx_power_body <- function(n, p0, p1, alpha) {
    se0 <- sqrt(p0 * (1 - p0) / n)
    se1 <- sqrt(p1 * (1 - p1) / n)
    delta <- if (cc) 0.5 / n else 0

    if (alternative == "greater") {
      z <- qnorm(1 - alpha)
      1 - pnorm((p0 + z * se0 - delta - p1) / se1)
    } else if (alternative == "less") {
      z <- qnorm(1 - alpha)
      pnorm((p0 - z * se0 + delta - p1) / se1)
    } else {
      z <- qnorm(1 - alpha / 2)
      pnorm((p0 - z * se0 + delta - p1) / se1) +
        1 - pnorm((p0 + z * se0 - delta - p1) / se1)
    }
  }

  exact_quantile_power_body <- function(n, p0, p1, alpha, midp) {
    rr <- get_rr(n, p0, alpha)

    if (rr$type == "greater") {
      k <- rr$k
      if (midp)
        pbinom(k - 1, n, p1, lower.tail = FALSE) +
          0.5 * dbinom(k, n, p1)
      else
        pbinom(k - 1, n, p1, lower.tail = FALSE)
    } else if (rr$type == "less") {
      k <- rr$k
      if (midp)
        pbinom(k - 1, n, p1) +
          0.5 * dbinom(k, n, p1)
      else
        pbinom(k, n, p1)
    } else {
      k_lo <- rr$k_lo
      k_hi <- rr$k_hi
      if (midp) {
        pbinom(k_lo - 1, n, p1) +
          0.5 * dbinom(k_lo, n, p1) +
          pbinom(k_hi - 1, n, p1, lower.tail = FALSE) +
          0.5 * dbinom(k_hi, n, p1)
      } else {
        pbinom(k_lo, n, p1) +
          pbinom(k_hi - 1, n, p1, lower.tail = FALSE)
      }
    }
  }

  cp_power_body <- function(n, p0, p1, alpha) {
    x <- 0:n
    reject <- if (alternative == "greater") {
      pbinom(x - 1, n, p0, lower.tail = FALSE) <= alpha
    } else if (alternative == "less") {
      pbinom(x, n, p0) <= alpha
    } else {
      2 * pmin(
        pbinom(x, n, p0),
        pbinom(x - 1, n, p0, lower.tail = FALSE)
      ) <= alpha
    }
    sum(dbinom(x[reject], n, p1))
  }

  ## ---- dispatch ----
  if (!exact) {
    power_body <- approx_power_body
    midp <- FALSE
  } else if (exact.method == "cp") {
    power_body <- cp_power_body
    midp <- FALSE
  } else if (exact.method == "midp") {
    power_body <- function(n, p0, p1, alpha)
      exact_quantile_power_body(n, p0, p1, alpha, midp = TRUE)
    midp <- TRUE
  } else {
    power_body <- function(n, p0, p1, alpha)
      exact_quantile_power_body(n, p0, p1, alpha, midp = FALSE)
    midp <- FALSE
  }

  ## ---- solve for missing parameter ----
  if (is.null(power)) {
    power <- power_body(n, p0, p1, sig.level)
  }

  if (is.null(n)) {
    f <- function(nn) power_body(nn, p0, p1, sig.level) - power
    n_lo <- 1L
    if (f(n_lo) >= 0) {
      n <- n_lo
    } else {
      n_hi <- 2L
      while (f(n_hi) < 0) {
        n_hi <- n_hi * 2L
        if (n_hi > max_n)
          stop("Required n exceeds 'max_n'", call. = FALSE)
      }
      n <- as.integer(ceiling(uniroot(f, c(n_lo, n_hi), tol = tol)$root))
    }
  }

  if (is.null(p1)) {
    f <- function(pp) power_body(n, p0, pp, sig.level) - power
    p1 <- uniroot(f, c(tol, 1 - tol), tol = tol)$root
  }

  if (is.null(p0)) {
    f <- function(pp) power_body(n, pp, p1, sig.level) - power
    p0 <- uniroot(f, c(tol, 1 - tol), tol = tol)$root
  }

  if (is.null(sig.level)) {
    f <- function(a) power_body(n, p0, p1, a) - power
    sig.level <- uniroot(f, c(tol, 1 - tol), tol = tol)$root
  }

  ## ---- return object ----
  method <- if (!exact)
    "One-sample proportion power calculation (normal approximation)"
  else if (exact.method == "cp")
    "One-sample proportion power calculation (exact Clopper--Pearson)"
  else if (exact.method == "midp")
    "One-sample proportion power calculation (exact binomial, mid-p)"
  else
    "One-sample proportion power calculation (exact binomial)"

  structure(
    list(
      n = n,
      p0 = p0,
      p1 = p1,
      sig.level = sig.level,
      power = power,
      alternative = alternative,
      note = "n is sample size for the single group",
      method = method
    ),
    class = "power.htest"
  )
}
