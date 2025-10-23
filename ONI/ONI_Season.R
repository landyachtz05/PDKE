
library(here)
library(dplyr)
library(lubridate)
library(here)

##### Converts monthly ONI ENSO values to wet/dry season values #####

BASE_DIR <- here()
CCVD_INPUTS_DIR <- file.path(BASE_DIR, "CCVD", "CCVD_INPUTS")
oni_file <- file.path(BASE_DIR, "ONI", "enso_oni_raw.csv")

cat("CCVD_INPUTS_DIR:", CCVD_INPUTS_DIR, "\n")
cat("oni_file:", oni_file, "\n")

# --- Read monthly ONI data ---
oni <- read.csv(oni_file)
oni$Date <- as.Date(oni$Date)
oni <- oni %>%
  mutate(
    Year = year(Date),
    Month = month(Date)
  )

# --- Custom helper for max absolute value with sign check ---
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

# --- Determine all relevant wet-season years ---
# We’ll start at the second year (since wet 1950 needs 1949 Oct–Dec)
years <- sort(unique(oni$Year))
results <- data.frame(Year = years, MEI_W = NA, MEI_D = NA)

for (i in seq_along(years)) {
  yr <- years[i]
  
  # Dry season: May–Sep of the same year
  dry_values <- oni %>%
    filter(Year == yr, Month >= 5, Month <= 9) %>%
    pull(ONI)
  results$MEI_D[i] <- round(season_max_value(dry_values), 1)
  
  # Wet season: Oct–Dec of previous year + Jan–Apr of this year
  wet_values <- oni %>%
    filter(
      (Year == yr - 1 & Month >= 10) |  # Oct–Dec of previous year
        (Year == yr & Month <= 4)         # Jan–Apr of current year
    ) %>%
    pull(ONI)
  results$MEI_W[i] <- round(season_max_value(wet_values), 1)
}

# --- Output ---
print(head(results, 10))
print(tail(results))

output_file <- file.path(CCVD_INPUTS_DIR, "ONI_Season.csv")
current_date <- format(Sys.Date(), "%Y%m%d")
backup_file <- file.path(CCVD_INPUTS_DIR, paste0("ONI_Season_backup_", current_date, ".csv"))

if (file.exists(output_file)) {
  file.copy(output_file, backup_file)
  cat("Existing ONI_Season.csv backed up to:", backup_file, "\n")
} else {
  cat("ONI_Season.csv does not exist, no backup created.\n")
}

cat("Saving to:", output_file, "\n")
write.csv(results, output_file, row.names = FALSE)
