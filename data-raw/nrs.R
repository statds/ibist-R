## data-raw/nrs.R

nrs <- read.csv("data-raw/nrs.csv", header = TRUE)

usethis::use_data(nrs, overwrite = TRUE)
