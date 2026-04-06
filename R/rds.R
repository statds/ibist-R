#' Neonatal Respiratory Distress Syndrome Data
#'
#' Aggregated counts of neonatal outcomes by birth weight,
#' surfactant use, and survival status.
#'
#' @format A data frame with 16 rows and 4 variables:
#' \describe{
#'   \item{bwt}{Birth weight category (factor with 4 levels).}
#'   \item{surf}{Surfactant use (Yes/No).}
#'   \item{death}{Outcome (1 = death; 0 = alive).}
#'   \item{count}{Number of infants in each group.}
#' }
#'
#' @details
#' The data are presented in aggregated form. Each row corresponds
#' to a combination of birth weight category, surfactant use, and
#' outcome, with \code{count} indicating the number of observations.
#'
#' Analyses should account for the grouped structure, e.g.,
#' via weighted models or binomial responses.
#'
#' @source
#' Classical neonatal RDS dataset (exact source to be specified).
"rds"
