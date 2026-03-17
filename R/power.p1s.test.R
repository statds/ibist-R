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
#' @param two.sided.small.tail Logical; for two-sided tests under the
#'   normal approximation, include the smaller tail in the power
#'   calculation. Setting FALSE reproduces the common textbook
#'   approximation that keeps only the dominant rejection tail.
#' @param tol Numerical tolerance used in root finding.
#' @param max_n Maximum allowable sample size when solving for \code{n}.
#' @param size.rule Character string specifying the rule used to select
#'   the sample size when \code{exact = TRUE}. Must be one of
#'   \code{"close"} (default) or \code{"minimal"}.
#'
#'   \code{"close"} selects the smallest sample size \eqn{n} such that
#'   the achieved significance level satisfies
#'   \eqn{\alpha_{\mathrm{act}} \le \alpha}, the power meets the target,
#'   and the achieved level is not too far below the nominal level:
#'   \deqn{\alpha_{\mathrm{act}} \ge c \alpha,}
#'   where \eqn{c =} \code{alpha.min.frac}. This rule avoids overly
#'   conservative designs and makes more efficient use of the type I
#'   error budget.
#'
#'   \code{"minimal"} selects the smallest sample size \eqn{n} such that
#'   the achieved significance level satisfies
#'   \eqn{\alpha_{\mathrm{act}} \le \alpha} and the power is at least the
#'   target value. This is the classical exact test design but may be
#'   conservative, with \eqn{\alpha_{\mathrm{act}} \ll \alpha}.
#' @param alpha.min.frac Numeric value in (0, 1) specifying the minimum
#'   fraction of the nominal significance level required for the achieved
#'   level when \code{size.rule = "close"}.
#'
#'   The achieved significance level must satisfy
#'   \deqn{\alpha_{\mathrm{act}} \ge
#'   \text{alpha.min.frac} \times \alpha.}
#'
#'   The default value \code{0.9} requires the achieved level to be at
#'   least 90\% of the nominal level. Smaller values allow more
#'   conservative designs, while larger values enforce closer agreement
#'   with the nominal level.
#' 
#'   See also \code{achieved.sig.level} in the returned value.
#' 
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
#' @importFrom stats dbinom pbinom qbinom
#' @importFrom stats dnorm pnorm qnorm
#' @importFrom stats uniroot density
#
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
  two.sided.small.tail = TRUE,
  tol = .Machine$double.eps^0.5,
  max_n = 1e7,
  size.rule = c("close", "minimal"),
  alpha.min.frac = 0.9                
) {
  alternative <- match.arg(alternative)
  exact.method <- match.arg(exact.method)
  size.rule <- match.arg(size.rule)

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

  midp = exact.method == "midp"

  ## ---- exact quantile rejection region cache ----
  rr_cache <- new.env(parent = emptyenv())

  get_rr <- function(n, p0, sig.level) {
    key <- paste(n, p0, sig.level, alternative, sep = "|")
    if (exists(key, envir = rr_cache, inherits = FALSE))
      return(get(key, envir = rr_cache, inherits = FALSE))

    if (alternative == "greater") {
      rr <- list(type = "greater",
                 k = qbinom(1 - sig.level, n, p0) + 1L)
    } else if (alternative == "less") {
      rr <- list(type = "less",
                 k = qbinom(sig.level, n, p0) - 1L)
    } else {
      rr <- list(
        type = "two.sided",
        k_lo = qbinom(sig.level / 2, n, p0) - 1L,
        k_hi = qbinom(1 - sig.level / 2, n, p0) + 1L
      )
    }

    assign(key, rr, envir = rr_cache)
    rr
  }

  ## ---- power bodies ----

  approx_power_body <- quote({
    se0 <- sqrt(p0 * (1 - p0) / n)
    se1 <- sqrt(p1 * (1 - p1) / n)
    delta <- if (cc) 0.5 / n else 0

    if (alternative == "greater") {
      z <- qnorm(1 - sig.level)
      1 - pnorm((p0 + z * se0 + delta - p1) / se1)
    } else if (alternative == "less") {
      z <- qnorm(1 - sig.level)
      pnorm((p0 - z * se0 - delta - p1) / se1)
    } else {
      z <- qnorm(1 - sig.level / 2)
      lower_tail <- pnorm((p0 - z * se0 - delta - p1) / se1)
      upper_tail <- 1 - pnorm((p0 + z * se0 + delta - p1) / se1)
      if (two.sided.small.tail) {
        lower_tail + upper_tail
      } else {
        dominant_upper <- p1 >= p0
        if (dominant_upper) upper_tail else lower_tail
      }
    }
  })

  exact_quantile_power_body <- quote({
    rr <- get_rr(n, p0, sig.level)

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
        pbinom(k - 1, n, p1) + 0.5 * dbinom(k, n, p1)
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
  })

  cp_power_body <- quote({
    x <- 0:n
    reject <- if (alternative == "greater") {
      pbinom(x - 1, n, p0, lower.tail = FALSE) <= sig.level
    } else if (alternative == "less") {
      pbinom(x, n, p0) <= sig.level
    } else {
      2 * pmin(
        pbinom(x, n, p0),
        pbinom(x - 1, n, p0, lower.tail = FALSE)
      ) <= sig.level
    }
    sum(dbinom(x[reject], n, p1))
  })

  ## ---- dispatch ----
  if (!exact) {
    power_body <- approx_power_body
    midp <- NA
  } else if (exact.method == "cp") {
    power_body <- cp_power_body
    midp <- NA
  } else {
    power_body <- exact_quantile_power_body
  }

  ## ---- solve for missing parameter ----
  if (is.null(power)) {
    power <- eval(power_body)
  }

  if (is.null(n)) {
      if (!exact) {
          n <- uniroot(function(n) eval(power_body) - power,
                       c(1, max_n), tol = tol,
                       extendInt = "upX")$root
      } else {
          power_at_n <- function(nn) {
              n <- as.integer(nn)
              eval(power_body)
          }
          
          alpha_at_n <- function(nn) {
              n <- as.integer(nn)
              p1 <- p0
              eval(power_body)
          }
          
          if (exact.method == "cp") {
              feasible <- function(nn) {
                  a <- alpha_at_n(nn)
                  pw <- power_at_n(nn)
                  if (size.rule == "minimal") {
                      pw >= power
                  } else {
                      (a <= sig.level) &&
                          (a >= alpha.min.frac * sig.level) &&
                          (pw >= power)
                  }
              }              
          } else {
              feasible <- function(nn) {
                  a <- alpha_at_n(nn)
                  pw <- power_at_n(nn)
                  if (size.rule == "minimal") {
                      pw >= power
                  } else {
                      (a <= sig.level) &&
                          (a >= alpha.min.frac * sig.level) &&
                          (pw >= power)
                  }
              }
          }
          
          ## ---- exponential bracketing + integer bisection ----
          ## n_lo <- 1L
          ## if (feasible(n_lo)) {
          ##     n <- n_lo
          ## } else {
          ##     n_hi <- 2L
          ##     while (!feasible(n_hi)) {
          ##         n_lo <- n_hi
          ##         n_hi <- n_hi * 2L
          ##         if (n_hi > max_n)
          ##             stop("Required n exceeds 'max_n'", call. = FALSE)
          ##     }
          ##     while (n_hi - n_lo > 1L) {
          ##         n_mid <- (n_lo + n_hi) %/% 2L
          ##         if (feasible(n_mid)) n_hi <- n_mid else n_lo <- n_mid
          ##     }
          ##     n <- n_hi
          ## }
          
          ## Step 1: monotone search
          hi <- 10
          while (power_at_n(hi) < power) hi <- hi * 2
          
          ## Step 2: local scan
          n_seq <- seq.int(max(1, hi/2), hi * 2)
          feas <- vapply(n_seq, feasible, logical(1))
          n <- n_seq[min(which(feas))]
          power <- power_at_n(n)
      }
  }

  
  if (is.null(p1)) {
      p1 <- uniroot(function(pp) eval(power_body) - power,
                    c(tol, 1 - tol), tol = tol)$root
  }

  if (is.null(p0)) {
      p0 <- uniroot(function(pp) eval(power_body) - power,
                    c(tol, 1 - tol), tol = tol)$root
  }
  
  if (is.null(sig.level)) {
      sig.level <- uniroot(function(a) eval(power_body) - power,
                           c(tol, 1 - tol), tol = tol)$root
  }

  achieved.sig.level <- if (exact) {
                            p1_save <- p1
                            p1 <- p0
                            out <- eval(power_body)
                            p1 <- p1_save
                            out
                        } else {
                            sig.level
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
  
  out <- list(
      n = n,
      p0 = p0,
      p1 = p1,
      sig.level = sig.level,
      power = power,
      alternative = alternative,
      method = method
  )
  if (exact) {
  out$achieved.sig.level <- achieved.sig.level
}

  structure(out, class = "power.htest")
}



## Two cautions
## First, if you adopt the “close” rule, there may be cases where no
## sample size in a modest range satisfies both the power target and the
## closeness requirement. In that case, you should either keep searching or
## fall back to "minimal" with a warning.
## Second, for exact.method = "midp", your current use of
## alpha_at_n(nn) <= sig.level treats the mid-plevel as if it were a
## strict size bound. That is a design choice, but mid-p p is usually
## motivated precisely because it relaxes conservativeness. So you should
## decide explicitly whether you want "midp" to be planned under nominal
## mid-p level or under conservative actual size.
