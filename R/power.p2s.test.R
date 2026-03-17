#' Power Calculation for Two-Sample Proportion Test with Unequal Group Sizes
#'
#' Computes power, sample size, detectable proportions, or significance
#' level for a two-sample test comparing proportions, allowing unequal
#' group sizes through a user-specified allocation ratio.
#'
#' Exactly one of \code{n}, \code{p1}, \code{p2}, \code{power}, or
#' \code{sig.level} must be \code{NULL}, and that parameter is solved
#' from the others.
#'
#' @param n Sample size in the first group. The second group size is
#'   \code{n * group.rate}.
#' @param p1 Proportion in the first group.
#' @param p2 Proportion in the second group.
#' @param sig.level Significance level (type I error rate), must be in
#'   \code{[0, 1]}.
#' @param power Target power of the test.
#' @param group.rate Ratio of sample sizes between the second and first
#'   groups, defined as \eqn{n_2 / n_1}. Defaults to 1 (equal group sizes).
#' @param alternative Character string specifying the alternative
#'   hypothesis, either \code{"two.sided"} or \code{"one.sided"}.
#' @param strict Logical; if \code{TRUE} and \code{alternative = "two.sided"},
#'   computes power accounting for both tails explicitly. Otherwise uses
#'   the standard one-sided approximation.
#' @param tol Tolerance for numerical root-finding when solving for an
#'   unknown parameter.
#'
#' @details
#' The calculation is based on the normal approximation to the binomial
#' distribution. Let \eqn{n_1 = n} and \eqn{n_2 = n \times group.rate}.
#' The test statistic uses a pooled variance under the null and an
#' unpooled variance under the alternative.
#'
#' When \code{strict = TRUE} for a two-sided test, power is computed as
#' the sum of tail probabilities. Otherwise, a one-sided approximation
#' is used.
#'
#' @return
#' An object of class \code{"power.htest"} with components:
#' \item{n}{Sample size in the first group.}
#' \item{p1}{Proportion in the first group.}
#' \item{p2}{Proportion in the second group.}
#' \item{sig.level}{Significance level.}
#' \item{power}{Power of the test.}
#' \item{alternative}{Type of alternative hypothesis.}
#' \item{note}{Clarifies that \code{n} refers to the first group size.}
#' \item{method}{Description of the method.}
#'
#' @note
#' This function generalizes \code{stats::power.prop.test()} by allowing
#' unequal group sizes through \code{group.rate}. When
#' \code{group.rate = 1}, the two functions are equivalent.
#'
#' @examples
#' # Power with unequal group sizes (n2 = 2 * n1)
#' power.p2s.test(n = 100, p1 = 0.3, p2 = 0.5, group.rate = 2)
#'
#' # Required sample size in group 1
#' power.p2s.test(p1 = 0.3, p2 = 0.5, power = 0.8, group.rate = 1.5)
#'
#' @export
power.p2s.test <- function(n = NULL, p1 = NULL, p2 = NULL,
                           sig.level = 0.05, power = NULL,
                           group.rate = 1,
                           alternative = c("two.sided", "one.sided"),
                           strict = FALSE, 
                           tol = .Machine$double.eps^0.25) {
    if (sum(sapply(list(n, p1, p2, power, sig.level), is.null)) != 1) 
        stop("exactly one of 'n', 'p1', 'p2', 'power', and 'sig.level' must be NULL")
    if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
        sig.level | sig.level > 1)) 
        stop("'sig.level' must be numeric in [0, 1]")
    alternative <- match.arg(alternative)
    tside <- switch(alternative, one.sided = 1, two.sided = 2)
    p.body <- if (strict && tside == 2) quote({
        n1 <- n; n2 <- n * grate
        qu <- qnorm(sig.level/tside, lower.tail = FALSE)
        d <- abs(p1 - p2)
        q1 <- 1 - p1; q2 <- 1 - p2
        pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
        qbar <- 1 - pbar
        v1 <- p1 * q1; v2 <- p2 * q2 
        se.num <- sqrt(pbar * qbar * (1 / n1 + 1 / n2))
        se.den <- sqrt(v1 / n1 + v2 / n2)
        pnorm( (d - qu * se.num) / se.den ) +
            pnorm((d + qu * se.num) / se.den, lower.tail=FALSE) 
    }) else quote({
        n1 <- n; n2 <- n * grate
        qu <- qnorm(sig.level/tside, lower.tail = FALSE)
        d <- abs(p1 - p2)
        q1 <- 1 - p1; q2 <- 1 - p2
        pbar <- (n1 * p1 + n2 * p2) / (n1 + n2)
        qbar <- 1 - pbar
        v1 <- p1 * q1; v2 <- p2 * q2 
        se.num <- sqrt(pbar * qbar * (1 / n1 + 1 / n2))
        se.den <- sqrt(v1 / n1 + v2 / n2)
        pnorm( (d - qu * se.num) / se.den )
    })
    if (is.null(power)) 
        power <- eval(p.body)
    else if (is.null(n)) 
        n <- uniroot(function(n) eval(p.body) - power, c(1, 1e+07), 
            tol = tol, extendInt = "upX")$root
    else if (is.null(p1)) {
        p1 <- uniroot(function(p1) eval(p.body) - power, c(0, 
            p2), tol = tol, extendInt = "yes")$root
        if (p1 < 0) 
            warning("No p1 in [0, p2] can be found to achieve the desired power")
    }
    else if (is.null(p2)) {
        p2 <- uniroot(function(p2) eval(p.body) - power, c(p1, 
            1), tol = tol, extendInt = "yes")$root
        if (p2 > 1) 
            warning("No p2 in [p1, 1] can be found to achieve the desired power")
    }
    else if (is.null(sig.level)) {
        sig.level <- uniroot(function(sig.level) eval(p.body) - 
            power, c(1e-10, 1 - 1e-10), tol = tol, extendInt = "upX")$root
        if (sig.level < 0 || sig.level > 1) 
            warning("No significance level [0, 1] can be found to achieve the desired power")
    }
    else stop("internal error", domain = NA)
    structure(list(n = n, p1 = p1, p2 = p2, sig.level = sig.level, 
        power = power, alternative = alternative, note = "n is number in the 1st group", 
        method = "Two-sample comparison of proportions power calculation"), 
        class = "power.htest")
}
