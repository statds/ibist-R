# data-raw/rds.R

# Read raw data
rds <- read.table("data-raw/rds.tab", header = TRUE)

# Clean column names (lowercase already OK)
names(rds) <- tolower(names(rds))

# Ensure types
rds$bwt   <- as.integer(rds$bwt)
rds$surf  <- as.integer(rds$surf)
rds$death <- as.integer(rds$death)
rds$count <- as.integer(rds$count)

# Optional: convert to factor with labels (recommended for users)
rds$bwt  <- factor(rds$bwt,
                   levels = 1:4,
                   labels = c("500-749g", "750-999g",
                              "1000-1249", "1250-1500g"))

rds$surf <- factor(rds$surf,
                   levels = c(0, 1),
                   labels = c("No", "Yes"))

# Save dataset
usethis::use_data(rds, overwrite = TRUE)
