# data-raw/mihdl.R

mihdl <- read.csv("data-raw/mihdl.csv", stringsAsFactors = FALSE)

# Type enforcement
mihdl$hdl <- as.integer(mihdl$hdl)
mihdl$gender <- factor(mihdl$gender, levels = c("Men", "Women"))
mihdl$cases <- as.integer(mihdl$cases)
mihdl$pyears <- as.numeric(mihdl$pyears)

# Optional: sort for consistency
mihdl <- mihdl[order(mihdl$hdl, mihdl$gender), ]

# Save to package
usethis::use_data(mihdl, overwrite = TRUE)
