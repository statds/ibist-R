#' Print an ibist Confidence Interval Object
#'
#' @param x An object of class \code{"ibist_ci"}.
#' @param digits Number of significant digits to print.
#' @param ... Reserved for future extensions.
#'
#' @export
print.ibist_ci <- function(x, digits = getOption("digits"), ...) {
  cat("\n", x$method, "\n", sep = "")
  cat("\n", 100 * x$conf.level, "% confidence interval:\n", sep = "")
  print(x$conf.int, digits = digits)
  invisible(x)
}

ci_table <- function(method, estimate, intervals, conf.level) {
  data.frame(
    method = method,
    estimate = estimate,
    lower = vapply(intervals, `[`, numeric(1), 1L),
    upper = vapply(intervals, `[`, numeric(1), 2L),
    conf.level = conf.level
  )
}
