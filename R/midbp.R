#' Blood Pressure Cohort Data (midbp)
#'
#' Aggregated data of cases and person-years by diastolic blood pressure
#' category and gender.
#'
#' @format A data frame with 14 rows and 4 variables:
#' \describe{
#'   \item{dbp}{Diastolic blood pressure category (1--7).}
#'   \item{gender}{Gender (Men, Women).}
#'   \item{cases}{Number of observed cases.}
#'   \item{pyears}{Person-years of follow-up.}
#' }
#'
#' @details
#' The dataset is aggregated by DBP category and gender, suitable for
#' rate modeling (e.g., Poisson regression with offset).
#'
#' @examples
#' data(midbp)
#' head(midbp)
#'
"midbp"
