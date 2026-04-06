## data-raw/midbp.R

midbp <- read.csv("data-raw/midbp.csv")

## Clean and standardize
midbp$dbp <- as.integer(midbp$dbp)
midbp$gender <- factor(midbp$gender, levels = c("Men", "Women"))
midbp$cases <- as.integer(midbp$cases)
midbp$pyears <- as.numeric(midbp$pyears)

## Optional: ensure ordering
midbp <- midbp[order(midbp$dbp, midbp$gender), ]

## Save
usethis::use_data(midbp, overwrite = TRUE)
