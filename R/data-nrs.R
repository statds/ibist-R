#' Non-restorative sleep and physical activity (Japan cohort study)
#'
#' A large observational dataset from a cohort study conducted in Japan
#' to examine the association between non-restorative sleep (NRS) and
#' physical activity, gender, and age. The data are used to illustrate
#' logistic regression modeling for a binary outcome in a large-sample
#' setting.
#'
#' @format
#' A data frame with 90,122 observations on the following variables:
#' \describe{
#'   \item{id}{Subject identifier.}
#'   \item{gender}{Gender of the subject (integer-coded).}
#'   \item{age}{Age in years in 2013.}
#'   \item{ex}{Indicator of regular exercise in 2013
#'   (integer-coded).}
#'   \item{pa}{Physical activity measure in 2013
#'   (integer-coded).}
#'   \item{nrs}{Indicator of non-restorative sleep in 2013
#'   (1 = presence, 0 = absence).}
#' }
#'
#' @details
#' Non-restorative sleep (NRS) is defined as a subjective feeling of lack
#' of refreshment on awakening and reflects qualitative aspects of sleep.
#' Hidaka et al. (2019) analyzed these data using logistic regression to
#' assess whether the probability of NRS is associated with physical
#' activity, gender, and age in a large cohort of adult subjects in
#' Japan. Within this package, the dataset is provided for methodological
#' illustration of binary regression models rather than for substantive
#' epidemiological inference.
#'
#' All variables are stored as integer codes. Missing values are
#' represented as \code{NA}.
#'
#' @source
#' Hidaka et al. (2019).
#'
"nrs"
