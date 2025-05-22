##### Converts monthly ONI ENSO values to season-year values #####

setwd("F:/PDKE/git/PDKE/CCVD/CCVD_INPUTS/")

# Read in ONI data
oni <- read.csv("enso_oni_raw.csv")

# Custom function: maximum absolute value with sign-check
season_max_value <- function(values) {
  values <- values[!is.na(values)]
  if (length(values) == 0) return(NA)
  
  abs_vals <- abs(values)
  max_val <- max(abs_vals)
  candidates <- values[abs_vals == max_val]
  
  if (length(candidates) > 1 && any(candidates > 0) && any(candidates < 0)) {
    return(mean(values))
  } else {
    return(candidates[1])
  }
}

# Initialize results
results <- data.frame(Year = oni$Year, MEI_W = NA, MEI_D = NA)

# Loop over each year
for (i in 1:nrow(oni)) {
  year <- oni$Year[i]
  
  ## --- Wet Season (Nov–Apr of next year) ---
  if (i < nrow(oni)) {
    wet_values <- c(
      oni$OND[i],     # Nov–Dec of this year
      oni$NDJ[i],     # Dec–Jan
      oni$DJF[i + 1], # Jan–Feb
      oni$JFM[i + 1], # Feb–Mar
      oni$FMA[i + 1], # Mar–Apr
      oni$MAM[i + 1]  # Apr–May
    )
    results$MEI_W[i] <- round(season_max_value(wet_values), 1)
  }
  
  ## --- Dry Season (May–Oct of same year) ---
  dry_values <- c(
    oni$AMJ[i],
    oni$MJJ[i],
    oni$JJA[i],
    oni$JAS[i],
    oni$ASO[i],
    oni$SON[i]
  )
  results$MEI_D[i] <- round(season_max_value(dry_values), 1)
}

# View the first few rows
head(results, 10)
tail(results)

write.csv(results, "ONI_Season2.csv")
