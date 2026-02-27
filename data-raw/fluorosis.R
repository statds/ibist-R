fluorosis <- data.frame(
  om = rep(c(0, 1), each = 12),
  fl_lev = rep(rep(1:3, each = 4), times = 2),
  amox = rep(rep(c(0, 1), each = 2), times = 6),
  fluorosis = rep(c(1, 0), times = 12),
  count = c(
    2,30,6,28,7,27,10,18,5,24,8,13,
    0,6,4,35,2,3,11,20,1,1,35,12
  )
)

# Convert to factors for modeling use
fluorosis$om        <- factor(fluorosis$om)
fluorosis$fl_lev    <- factor(fluorosis$fl_lev)
fluorosis$amox      <- factor(fluorosis$amox)
fluorosis$fluorosis <- factor(fluorosis$fluorosis)

usethis::use_data(fluorosis, overwrite = TRUE)
