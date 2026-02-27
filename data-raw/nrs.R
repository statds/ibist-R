## data-raw/nrs.R
nrs_raw <- read.csv("data-raw/nrs.csv", header = TRUE)

# Keep only required columns
nrs <- subset(
  nrs_raw,
  select = c(id, Gender, Age_2013, EX_2013, PA_2013, NRS_2013)
)

# Standardize names to lower case
names(nrs) <- tolower(names(nrs))

# Remove "_2013" suffix
names(nrs) <- sub("_2013$", "", names(nrs))

usethis::use_data(nrs, overwrite = TRUE)
