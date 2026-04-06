# data-raw/depress.R
# Purpose: Process raw Kalmbach et al. (2018) data into analysis dataset

library(haven)

# -----------------------------
# READ RAW DATA
# -----------------------------
unzip("data-raw/kalmbach2018.zip", exdir = "data-raw")
data <- read_sav("data-raw/S1 dataset.sav")

# -----------------------------
# EXTRACT RELEVANT VARIABLES
# -----------------------------
dat <- data[,c(5,6, 7, 9, 12, 68, 70, 159,
               10741, 10814, 10831, 10834,
               10850, 10851, 10852,
               10825, 10826, 10847, 10820)]

names(dat) <- c(
  "age", "gender", "marstat", "edu", "race",
  "sol", "waso", "bstresscount",
  "stressorA", "ies", "depressbase",
  "noctwake", "dummyinsom", "dummynp",
  "dummynpinsom", "depress1yr",
  "depress2yr", "depress1or2yr", "IEStert"
)

## -----------------------------
## APPLY INCLUSION CRITERIA
## -----------------------------
newdat2 <- dat[dat$depressbase < 10,] # 2665 obs (non-depressed at baseline)
newdat2 <- newdat2[newdat2$age > 17,] # 2664 obs (age 18+)
newdat2 <- newdat2[!is.na(newdat2$gender),] # 2662 obs (non-missing gender)
newdat2 <- newdat2[!is.na(newdat2$ies),] # 2083 obs (non-missing cognitive intrusions)
newdat2 <- newdat2[!is.na(newdat2$waso),] # 2078 obs (non-missing insomnia)
newdat2 <- newdat2[newdat2$stressorA >= 1,] # 2078 obs (had at least 1 stressor in past year)
newdat2 <- newdat2[newdat2$noctwake <= 90,] # 1730 obs (no extreme outliers in noctural wakefulness)

# if want to get same analytic sample as paper, additionally run this line:
# newdat2<- newdat2[!is.na(newdat2$IEStert),] #1126 # removes middle cognitive tertial 
#                                                     equiv. to removing IESgrp==1 below 

# -----------------------------
# DERIVED VARIABLES
# -----------------------------

# create cognitive intrusion variable
newdat2$IESgrp <- rep(0, nrow(newdat2))      # 0=lowest IES tertial
newdat2$IESgrp[newdat2$IEStert == 1] <- 2    # 2=highest IES tertial
newdat2$IESgrp[is.na(newdat2$IEStert)] <- 1  # 1=middle IES tertial

# create insomnia variable
newdat2$insom <- rep(0, nrow(newdat2))
newdat2$insom[newdat2$sol > 30] <- 1
newdat2$insom[newdat2$waso > 30] <- 1

# -----------------------------
# FINAL DATASET
# -----------------------------
depress <- newdat2[,c(
  "age", "gender",
  "stressorA",
  "depress1yr", "depress2yr", "depress1or2yr",
  "IEStert", "IESgrp", "insom"
)]

# -----------------------------
# SAVE
# -----------------------------
depress <- as.data.frame(depress)
usethis::use_data(depress, overwrite = TRUE)
