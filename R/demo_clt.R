#' Demonstrate the Central Limit Theorem
#'
#' The \code{demo_clt()} function generates plots to illustrate the
#' Central Limit Theorem (CLT) using a specified random number generator.
#' The function displays standardized sampling distributions for
#' different sample sizes and overlays the standard normal density.
#'
#' @param rng A random number generator function taking the sample size
#'   as its first argument (e.g., \code{runif}, \code{rnorm},
#'   \code{rgamma}).
#' @param n A numeric vector of sample sizes (e.g., \code{c(5, 10, 20,
#'   40)}).
#' @param nrep The number of repetitions for generating sample means
#'   (default is 10000).
#' @param ... Additional arguments passed to the random number generator
#'   (e.g., \code{shape} and \code{rate} for \code{rgamma}).
#' @param mean The theoretical mean of the distribution. If \code{NULL},
#'   it is estimated from a large Monte Carlo sample.
#' @param sd The theoretical standard deviation of the distribution.
#'   If \code{NULL}, it is estimated from a large Monte Carlo sample.
#'
#' @return A \code{ggplot2} object showing the standardized sampling
#'   distributions for different sample sizes, compared against the
#'   standard normal curve.
#'
#' @examples
#' set.seed(123)
#' demo_clt(runif, n = c(5, 10, 20, 40), min = 0, max = 1)
#'
#' demo_clt(
#'   rgamma,
#'   n = c(5, 10, 20, 40),
#'   shape = 2, rate = 1,
#'   mean = 2, sd = sqrt(2)
#' )
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_histogram geom_line aes
#' @importFrom ggplot2 facet_wrap labs theme_minimal after_stat
#' @export
demo_clt <- function(
  rng,
  n,
  nrep = 10000,
  ...,
  mean = NULL,
  sd = NULL
) {
  ## ---- basic validation ----
  if (!is.function(rng)) {
    stop("The argument 'rng' must be a function.", call. = FALSE)
  }

  if (!is.numeric(n) || any(n <= 0)) {
    stop(
      "The argument 'n' must be a numeric vector of positive values.",
      call. = FALSE
    )
  }

  if (!is.numeric(nrep) || length(nrep) != 1L || nrep <= 0) {
    stop(
      "The argument 'nrep' must be a positive integer.",
      call. = FALSE
    )
  }

  ## ---- estimate mean and sd if needed ----
  if (is.null(mean) || is.null(sd)) {
    sample_data <- rng(100000, ...)
    mean <- base::mean(sample_data)
    sd <- stats::sd(sample_data)
  }

  ## ---- generate standardized sample means ----
  results <- vector("list", length(n))
  names(results) <- as.character(n)

  for (size in n) {
    sample_means <- replicate(
      nrep,
      base::mean(rng(size, ...))
    )

    standardized_means <-
      (sample_means - mean) / (sd / sqrt(size))

    results[[as.character(size)]] <- data.frame(
      StandardizedMean = standardized_means,
      SampleSize = factor(size)
    )
  }

  data <- do.call(rbind, results)

  ## ---- standard normal reference ----
  x_vals <- seq(-4, 4, length.out = 200)
  normal_data <- data.frame(
    x = x_vals,
    y = stats::dnorm(x_vals)
  )

  ## ---- plot ----
  ggplot2::ggplot(
  data,
  ggplot2::aes(x = .data$StandardizedMean)
  ) +
  ggplot2::geom_histogram(
    ggplot2::aes(y = ggplot2::after_stat(density)),
    bins = 30,
    color = "black",
    fill = "skyblue"
  ) +
  ggplot2::geom_line(
    data = normal_data,
    ggplot2::aes(x = .data$x, y = .data$y),
    color = "red",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  ggplot2::facet_wrap(~ SampleSize) +
  ggplot2::labs(
    title = "Demonstrating the Central Limit Theorem",
    x = "Standardized Sample Mean",
    y = "Density"
  ) +
  ggplot2::theme_minimal()
}
