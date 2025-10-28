# To run this, either run PDKESite/server.R to get a GUI which allows you to
# select an area on which to run this script, or, using one of the examples
# below, copy and paste one of the lines into the command line/terminal window.
# either one will run this script followed by the CCVD_portfolio_ppt.R script.
# A run generally takes about 40 minutes.
#
# No matter which method you choose, you need to edit the credentials.json file
# to set the variables correctly for your environment.
# Set ENV_TYPE to either "linux" or "windows", depending on the box you will be running this on.
# For the PROJ_LIB_VAL field search your system for a file called proj.db and use that path as the value.
# If you have none, install proj: https://proj.org/en/stable/index.html
# If that us set incorrectly, you will get the error message:
# “proj_create: Cannot find proj.db”
#
# to test via command line, no interface needed
# /Library/Frameworks/R.framework/Resources/Rscript /Users/jgeis/Work/PDKE/CCVD_portfolio_content.R 'jgeis@hawaii.edu' '/Users/jgeis/Work/PDKE/PDKESite/Shapefiles/SelectedPolygon/Molokai_2025_02_04_09_31_29.shp' 'Molokai' 'Molokai' 'Molokai' 'MO'"
# /Library/Frameworks/R.framework/Resources/bin/Rscript /Users/jgeis/Work/PDKE/CCVD_portfolio_content.R 'jgeis26@gmail.com' '/Users/jgeis/Work/PDKE/PDKESite/Shapefiles/SelectedPolygon/Kawela Ahupuaa_2025_04_21_11_48_04.shp' 'Kawela Ahupuaa' 'Kawela' 'Oahu' 'OA'

# for this to run, user must manually install proj: https://proj.org/en/9.3/about.html
start_time <- Sys.time()

library(gstat)
library(raster)
library(sp) 
#library(maptools)
library(RColorBrewer)
library(gridExtra)
library(ggplot2)
library(grid)
library(data.table)
library(devtools)
library(DescTools)
library(lubridate)
library(latticeExtra)
library(rasterVis)
library(plotrix)
library(plyr)
library(dplyr)
library(xts)
library(timeSeries)
library(ggfortify)
library(changepoint)
library(scales)
library(reshape2)
library(hydroTSM)
library(tiff)
library(lmomco)
library(parallel)
library(SPEI)
library(sf)
library(ggpubr)
library(terrainr)
library(ggmap)
library(ggthemes)
library(zoo)
library(classInt)
library(jsonlite)
library(here)


##########   This code will the generate the inputs for a CCVD portfolio
##########   Required Inputs: 1. Folder destinations 2. Measurement units 3. Shapefiles 4. List of files and output names
##########   Searchable sections of the Code:
##########   Maps - Elevation - Mean Climate - Downscaling - Rainfall Extract - SPI 1 - SPI 2 - MEI

##########################################################################################################################

print("In CCVD_portfolio_content.R")

read_credentials <- function(filepath) {
  tryCatch({
    credentials <- fromJSON(filepath)
    return(credentials)
  }, error = function(e) {
    print(paste("Error reading credentials file:", e$message))
    return(NULL) # Or handle the error as you see fit
  })
}

# default values for prod
ENV_TYPE = "windows"
PROJ_LIB_IN <- "/usr/share/proj/"
RSCRIPT_PATH <- "/usr/bin/Rscript"
BASE_DIR <- here() # Gets the project root
print(paste("RSCRIPT_PATH:", RSCRIPT_PATH))
print(paste("BASE_DIR:", BASE_DIR))

credentials_file <- paste0(BASE_DIR, "/credentials.json")
creds <- read_credentials(credentials_file)
if (!is.null(creds)) {
  ENV_TYPE <- creds$ENV_TYPE
  PROJ_LIB_IN <- creds$PROJ_LIB_VAL
  RSCRIPT_PATH <- creds$RSCRIPT_PATH
} else {
  print("Credentials could not be loaded")
}
print(paste("PROJ_LIB_IN:", PROJ_LIB_IN))
print(paste("BASE_DIR:", BASE_DIR))
print(paste("RSCRIPT_PATH:", RSCRIPT_PATH))
Sys.setenv(PROJ_LIB = PROJ_LIB_IN)


# non-user provided values
WORKING_DIR <- paste0(BASE_DIR, "/CCVD/MINI_Phase2/")
print(WORKING_DIR)
setwd(WORKING_DIR)               # WORKING DIRECTORY
INPUTS_FOLDER <- paste0(BASE_DIR, "/CCVD/CCVD_INPUTS/")       # INPUT FOLDER
OUTPUTS_FOLDER <- paste0(BASE_DIR, "/CCVD/CCVD_OUTPUTS/")       # OUTPUT FOLDER
MYSCRIPT_PATH = paste0(BASE_DIR, "/CCVD_portfolio_ppt.R")

#datetime_str <- format(Sys.time(), "%Y_%m_%d_%H_%M_%S")

# default values for testing, please don't check in modified values
email <- 'djford@hawaii.edu';
NP_FILE <- 'D:/PDKE/git/PDKE/PDKESite/Shapefiles/UserDefinedPolygon/Moaka Property_2025_07_08_14_20_29.shp';
NM <- 'Moaka Property';
NM_s <- 'Moaka';
ILE <- 'Oahu';
ILE_s <- 'OA';

# Get the command-line arguments passed from the main script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  email <- args[1];
  NP_FILE <- args[2];
  NM <- args[3];
  NM_s <- args[4];
  ILE <- args[5];
  ILE_s <- args[6];
} else {
  print("No command line args provided.  Using default values.")
}

# Extract the date and time string from the passed-in shape file using a regular expression
datetime_str <- sub(".*_(\\d{4}_\\d{2}_\\d{2}_\\d{2}_\\d{2}_\\d{2})\\.shp$", "\\1", NP_FILE)
PROJECT_WITH_DATE = paste0(NM_s, "_", datetime_str)

options(warn = -1)  # Suppress all warnings
#options(warn = 0)  # Re-enable warnings
debug_print <- function(content) {
  # would normally do a print or cat w/o using stderr, but this the only thing
  # that gets the shiny server on prod to actually print something.
  cat(file = stderr(), PROJECT_WITH_DATE, " PDKE: ", content, "\n") # prod
  #print(paste0(PROJECT_WITH_DATE, ", PDKE: ", content)) # dev
}

debug_print(paste0("args:", length(args)))
debug_print(paste("1,BASE_DIR: ", BASE_DIR))
debug_print(paste("1,WORKING_DIR: ", WORKING_DIR))
debug_print(paste("1,INPUTS_FOLDER: ", INPUTS_FOLDER))
debug_print(paste("1,OUTPUTS_FOLDER: ", OUTPUTS_FOLDER))

# These are now provided via command line or gui, do not hardcode
#NP_DIR <- paste0(INPUTS_FOLDER, "waikiki_watershed/")
#NP_FILE <- paste0(NP_DIR, "waikiki_watershed.shp")
#NM <- "Waikiki Watershed"
#NM_s <- "Waikiki"
#ILE <- "Oahu"
#ILE_s <- "OA"

debug_print(paste0("email: ", email))
debug_print(paste0("NP_FILE: ", NP_FILE))
debug_print(paste0("NM: ", NM))
debug_print(paste0("NM_s: ", NM_s))
debug_print(paste0("ILE: ", ILE))
debug_print(paste0("ILE_s: ", ILE_s))

###########################################################################################################################

###########################################################################################################################

##########   SET MEASUREMENT UNITS

TUnit = "\u00B0C"
TUnit2 = " \u00B0C"
RFUnit = " mm"
RFUnit2 = "mm"
ELUnit = " m"
ELUnit2 = "m"

TUnit = "\u00B0F"
TUnit2 = " \u00B0F"
RFUnit = " in"
RFUnit2 = "in"
ELUnit = " ft"
ELUnit2 = "ft"

############################################################################################################################

##########   Mean Climate Data and EXAMPLE RASTER To Set Spatial Extent
debug_print("define Mean_CLIM")
Mean_CLIM = dir(
  paste0(INPUTS_FOLDER, "Mean Climate"),
  pattern = "*x.adf",
  recursive = TRUE,
  full.names = TRUE
)
EXAMP <- raster(Mean_CLIM[71])
Mean_CLIM[71]
plot(EXAMP)

############################################################################################################################
##########   Coastal Shape Files
debug_print("Coastal Shape Files")

COAST_PATH <- paste0(INPUTS_FOLDER, "COAST/Coast_2/coast_geo_shp.dbf")
debug_print(paste("1,COAST_PATH", COAST_PATH))
Coast <- st_read(COAST_PATH)
plot(Coast[which(Coast$COAST_GEO_ == 10), ])

# Define the cropping extent as a bounding box
bbox <- st_bbox(c(xmin = -160.0, xmax = -154.8066, ymin = 18.91069, ymax = 22.23), crs = st_crs(Coast))
# Crop the `sf` object
Coast_Crop <- st_crop(Coast, bbox)
plot(Coast_Crop)

# BIG ISLAND
Coast_BI <- st_transform(Coast_Crop[which(Coast_Crop$COAST_GEO_ == 14),], st_crs(EXAMP))
plot(Coast_BI)

# MAUI
Coast_MN <- st_transform(Coast_Crop[which(Coast_Crop$COAST_GEO_ == 10),], st_crs(EXAMP))
plot(Coast_MN)

#Molokai
Coast_MO <- st_transform(Coast_Crop[which(Coast_Crop$COAST_GEO_ == 9),], st_crs(EXAMP))

#Lanai
Coast_LA <- st_transform(Coast_Crop[which(Coast_Crop$COAST_GEO_ == 11),], st_crs(EXAMP))

#Oahu
Coast_OA <- st_transform(Coast_Crop[which(Coast_Crop$COAST_GEO_ == 5),], st_crs(EXAMP))

#KAUI
Coast_KA <- st_transform(Coast_Crop[which(Coast_Crop$COAST_GEO_ == 2),], st_crs(EXAMP))

#Kahoolawe 
Coast_KO <- st_transform(Coast_Crop[which(Coast_Crop$COAST_GEO_ == 13),], st_crs(EXAMP))

############################################################

debug_print(paste("1,NP_FILE", NP_FILE))
debug_print(file.exists(NP_FILE))

# this now gets passed in via command line or GUI, do not hard-code
#NP_ALL <- readOGR("D:/PDKE/CCVD/sites/honouliuli_national_historic_site.shp") 
NP_ALL <- st_read(NP_FILE)

# these all get passed in via command line or GUI, do not hard-code
#HALE <- NP_ALL
#NM <- "Honouliuli National Historic Site"
#NM_s <- "HONO"
#ILE <- "Oahu"
#ILE_s <- "OA"

debug_print("1")

HALE <- NP_ALL

# Match AOI CRS with Coastline Polygon
HALE <- st_transform(HALE, st_crs(EXAMP))

plot(HALE)

# Ensure HALE is a POLYGON and extract coordinates properly
HALE_df <- HALE %>%
  st_cast("POLYGON") %>%  # Ensure proper geometry type
  st_coordinates() %>%    # Extract coordinates
  as.data.frame() %>%
  dplyr::rename(long = X, lat = Y)  # Rename for clarity

# Add group information for ggplot
HALE_df$group <- HALE_df$L1 

##########   Calculate area of AOI and buffer if it's too small
debug_print("Calculate area of AOI and buffer if it's too small")

# Reproject to UTM Zone 4 (if HALE is an sf object)
HALE_planar <- st_transform(HALE, crs = "+proj=utm +zone=4 +datum=WGS84 +units=m +no_defs")
# Calculate area in square meters
HALE_planar$area_m2 <- st_area(HALE_planar)
# Convert area to acres
HALE_planar$area_acres <- as.numeric(HALE_planar$area_m2) / 4046.85642
HALE_planar$area_acres

if (as.numeric(HALE_planar$area_acres) < 100) {
  
  # Define the target area in square meters
  desired_area_m2 <- 100 * 4046.85642
  
  # Calculate the side length of a square with this area
  side_length <- sqrt(desired_area_m2)
  
  # Get the centroid of the HALE polygon
  centroid <- st_centroid(HALE_planar)
  
  # Extract centroid coordinates
  centroid_coords <- st_coordinates(centroid)
  
  # Create a square polygon around the centroid
  square_coords <- matrix(c(
    centroid_coords[1] - side_length / 2, centroid_coords[2] - side_length / 2,
    centroid_coords[1] + side_length / 2, centroid_coords[2] - side_length / 2,
    centroid_coords[1] + side_length / 2, centroid_coords[2] + side_length / 2,
    centroid_coords[1] - side_length / 2, centroid_coords[2] + side_length / 2,
    centroid_coords[1] - side_length / 2, centroid_coords[2] - side_length / 2 # Close the polygon
  ), ncol = 2, byrow = TRUE)
  
  # Convert to an `sf` polygon
  square_buffer <- st_sfc(st_polygon(list(square_coords)), crs = st_crs(HALE_planar))
  
  # Convert to sf object
  square_buffer <- st_sf(geometry = square_buffer)
  
  # Assign the square buffer to HALE
  HALE <- square_buffer
  
  # Plot the square buffer in red
  plot(st_geometry(square_buffer), col = 'red', border = 'darkred', main = "HALE and Square Buffer")
  
  # Then plot the HALE polygon on top in blue
  plot(st_geometry(HALE_planar), col = 'blue', border = 'darkblue', add = TRUE)
}


##########   Set island boundary
debug_print("Set island boundary")

# Set island boundary and plot the corresponding coast
if (ILE_s == "BI") { 
  plot(st_geometry(Coast_BI), main = ILE)  # Plot the island boundary
  plot(st_geometry(HALE), add = TRUE, col = "red")  # Overlay HALE geometry
} 
if (ILE_s == "MN") { 
  plot(st_geometry(Coast_MN), main = ILE)
  plot(st_geometry(HALE), add = TRUE, col = "red")  
} 
if (ILE_s == "MO") { 
  plot(st_geometry(Coast_MO), main = ILE)
  plot(st_geometry(HALE), add = TRUE, col = "red")  
} 
if (ILE_s == "LA") { 
  plot(st_geometry(Coast_LA), main = ILE)
  plot(st_geometry(HALE), add = TRUE, col = "red")  
} 
if (ILE_s == "KO") { 
  plot(st_geometry(Coast_KO), main = ILE)
  plot(st_geometry(HALE), add = TRUE, col = "red")  
} 
if (ILE_s == "OA") { 
  plot(st_geometry(Coast_OA), main = ILE)
  plot(st_geometry(HALE), add = TRUE, col = "red")  
} 
if (ILE_s == "KA") { 
  plot(st_geometry(Coast_KA), main = ILE)
  plot(st_geometry(HALE), add = TRUE, col = "red")  
}

##################################################
##########   Define Units
debug_print("Define Units")

##########   Assign AOI shapefile coordinates
HALE <- st_transform(HALE, crs = crs(EXAMP))
HALE <- as_Spatial(HALE)

##########   SHAPE FILES FOR ANALYSIS
UNIT_X <-   c(HALE, HALE)

##########   ISLAND - NEED TO CHANGE MANUALLY BASED ON ISLAND
UNIT_I <-  c(ILE_s, ILE_s)

##########  UNIT NAME
##########  UNIT NAME (For CCVD Narrative)
UNIT_N <- as.vector(as.character(c(NM, NM)))

##########  Short Name (For Figure Titles)
UNIT_Ns <- as.vector(as.character(c(NM_s, NM_s)))

#20,24,-26 34


##########   COUNT INPUTS
debug_print("COUNT INPUTS")

UNIT_Shape <- sum(!is.na(UNIT_X))
UNIT_C <- sum(!is.na(UNIT_N))  # This will mark the end of the Loop
UNIT_Island <- sum(!is.na(UNIT_I))
ShortNm <- sum(!is.na(UNIT_Ns))


##########   CONFIRM INPUTS ARE THE SAME

# UNIT_C %>%
# ShortNm
# UNIT_Shape
# UNIT_Island

##########  Mean Rainfall DATA
debug_print("Define MeanRF_ALL")

Mean_RF_Data_dir <- paste0(INPUTS_FOLDER, "Mean_RF_Data/StateMaps/")
MeanRF_ALL = dir(
  Mean_RF_Data_dir,
  pattern = "*x.adf",
  recursive = T,
  full.names = T
)

#
# ## Find state-wide mean monthly rainfall stats (min, mean, max, percentiles)
# # loop through month maps and calculate stats, also make average monthly rainfall map
# n<-1
# fq<-data.frame(quantile(c(0:0)))
#
# sm<-stack()
#
# for (i in MeanRF_ALL[1:12]) {
#   m<-raster(i)
#
#   # convert mm to inches
#   if(RFUnit == " in") {m  <- m  * 0.0393701}
#
#   # calculate quantiles and put in dataframe
#   mq<-data.frame(quantile(m))
#   colnames(mq)<-n
#   fq<-cbind(fq,mq)
#   n<-n+1
#
#   # add raster to stack
#   sm<-stack(sm, m)
# }
#
# sm
#
# rm(mean)
#
# fq<-round(fq[2:ncol(fq)], 2)
# fq
# write.csv(fq, paste0(INPUTS_FOLDER,"Mean_RF_Data/stateRFM_quantiles_input.csv"))

# TODO: Move this stuff to where it's actually getting used (around line 5816)
# import static data table of statewide monthly rainfall quantiles
Mean_RF_Data_file <- paste0(Mean_RF_Data_dir, "stateRFM_quantiles_input.csv")
debug_print(paste0("Mean_RF_Data_file: ", Mean_RF_Data_file))
fq <- data.frame(read.csv(
  Mean_RF_Data_file,
  header = TRUE
))
stateRF <- rowMeans(fq[2:ncol(fq)])
# stateRFM<-calc(sm, mean, na.rm=TRUE)
# stateRFM

# # make dataframe from average monthly rainfall map
# stateRFMd<-as.data.frame(stateRFM, na.rm=T)
# stateRFMd<-data.frame(stateRFMd[1])
# write.csv(stateRFMd, "stateRFMd_input.csv")

# read in static state-wide rainfall raster values table
stateRFMd_file <- paste0(Mean_RF_Data_dir, "stateRFMd_input.csv")
debug_print(paste0("stateRFMd_file: ", stateRFMd_file))
stateRFMd <- read.csv(stateRFMd_file)

u <- 1

##########   FIRE OCCURRENCE AND RISK
#Fire Occurrence Shape 2022 (From Clay)
#FIRE_Shape <- readOGR(paste0(INPUTS_FOLDER,"StateFire_1999/fires_1999_2022/fires_1999_2022.shp"))
FIRE_Shape_File <- paste0(INPUTS_FOLDER, "StateFire_1999/fires_1999_2022/fires_1999_2022.shp")
debug_print(paste("17,FIRE_Shape_File", FIRE_Shape_File))
FIRE_Shape <- st_read(FIRE_Shape_File)
FIRE_Shape_T <- st_transform(FIRE_Shape, st_crs(EXAMP))
# Fire_Mask <- mask(x = FIRE_Shape_T, mask = UNIT_X[[u]])
# Fire_Crop <- crop(x = FIRE_Shape_T, y = extent(UNIT_X[[u]]))

# #Fire risk categories 2023 (From Clay)
# FRISK <- raster(paste0(INPUTS_FOLDER,"FireRisk/Avg_Landscape_Fire_Risk_Hawaii_2023/Avg_Landscape_Fire_Risk_Hawaii_fire_risk_categories_2023.tif"))
# crs(FRISK) <- "+proj=longlat +datum=WGS84 +no_defs"
#
# FRISK_Sh
# FRISK_Shape_T <- spTransform(FRISK_Shape, crs(EXAMP))
# crs(FRISK_Shape_T) <- CRS(proj4string(EXAMP))

##########   Forest Roads
debug_print(paste0("22.2, INPUTS_FOLDER: ", INPUTS_FOLDER))

#F_Roads <- readOGR(paste0(INPUTS_FOLDER, "Forestry Roads/Forestry_Roads.shp"))
F_Roads_File <- paste0(INPUTS_FOLDER, "Forestry Roads/Forestry_Roads.shp")
debug_print(paste("23,F_Roads_File", F_Roads_File))
F_Roads <- st_read(F_Roads_File)
F_Roads <- st_transform(F_Roads, st_crs(EXAMP))
plot(F_Roads[1])

##########   Digital Elevation Model

ELEV2 = dir(
  paste0(INPUTS_FOLDER, "ned_dem/"),
  pattern = "*x.adf",
  recursive = T,
  full.names = T
)
# ### if doing portfolio for whole island (needed for Big Island at least)
#   ELEV2 = dir(paste0(INPUTS_FOLDER), pattern="*ned_dem_crs.tif", recursive=T, full.names=T)
ELEV <- raster(ELEV2[1])
if (ELUnit == " ft") { ELEV = ELEV * 3.28084 }
crs(ELEV) <- "+proj=longlat +datum=WGS84 +no_defs"

##########  Hillshade
HS <- raster(paste0(INPUTS_FOLDER, "hawaii_hillshade.tif"))
crs(HS) <- "+proj=longlat +datum=WGS84 +no_defs"

##########  Landcover

LC2 = dir(
  paste0(INPUTS_FOLDER, "Landcover/"),
  pattern = "*LCMAP_HI_2020_V10_LCPRI_crs.tif$",
  recursive = T,
  full.names = T
)
LC <- raster(LC2)
crs(LC) <- "+proj=longlat +datum=WGS84 +no_defs"

##########  Moku and Ahupuaa
debug_print("Moku and Ahupuaa")

MOKU_file <- paste0(INPUTS_FOLDER, "Moku.shp")
debug_print(paste0("MOKU_file: ", MOKU_file))
MOKU <- st_read(MOKU_file) %>%
  st_transform(st_crs(EXAMP))

AHU_file <- paste0(INPUTS_FOLDER, "Ahupuaa3.shp")
debug_print(paste0("AHU_file: ", AHU_file))
AHU <- st_read(AHU_file) %>%
  st_transform(st_crs(EXAMP))

##########  Streams and Aquifers
debug_print("Streams and Aquifers")

STRM <- st_read(paste0(INPUTS_FOLDER, "NHD_H_Hawaii_State_Shape/NHD_Flowlines2.shp")) %>%
  st_transform(st_crs(EXAMP))

AQU <- st_read(paste0(INPUTS_FOLDER, "Aquifers/DOH_Aquifers_type_status.shp")) %>%
  st_transform(st_crs(EXAMP))

##########  Downscaling
debug_print("Downscaling")

# should be 132 DyDs_files. If TIF's were opened in arcgis, extra files may have been created
DyDs_files = dir(
  paste0(INPUTS_FOLDER, "HI_Downscaling_Final/"),
  recursive = T,
  full.names = T,
  pattern = "\\.tif$"
)
D2_files = dir(
  paste0(INPUTS_FOLDER, "Lulin_Downscaling/"),
  recursive = T,
  full.names = T
)

##########   Month Year Rainfall Maps (Frazier et al., 2016)
debug_print("Month Year Rainfall Maps (Frazier et al., 2016)")

# RF_Map_Path_A <- ("D:/PDKE/CCVD/data/production/rainfall/legacy/month/statewide/data_map/")
RF_Map_Path_A_DIR <- paste0(BASE_DIR, "/CCVD/data/production/rainfall/legacy/month/statewide/data_map/")
debug_print(paste("27,RF_Map_Path_A_DIR: ", RF_Map_Path_A_DIR))
RF_Map_Path_A <- (RF_Map_Path_A_DIR)

##########   Month-Year Rainfall Maps New Lucas et al., In Review

# Need to download maps from HCDP
# Legacy maps 1920-1989
# Matty maps 1990-NRT
# Put in INput folder
#For now Read this path in seperate as not in input folder
#RF_Map_Path <- ("D:/PDKE/CCVD/NEW RF MAPS/statewide_rf_mm/rf_mm/")
RF_Map_Path_DIR <- paste0(BASE_DIR, "/CCVD/NEW_RF_MAPS/statewide_rf_mm/rf_mm/")
RF_Map_Path <- (RF_Map_Path_DIR)
debug_print(paste("28,RF_Map_Path: ", RF_Map_Path))

##########  Multi-Variate ENSO Index

#Replace this with the ONI Dataset at some point
MEI <- read.csv(paste0(INPUTS_FOLDER, "ONI_Season.csv"), sep = ",")
ONI_raw <- read.csv(paste0(BASE_DIR, "/ONI/enso_oni_raw.csv"))

########## Month Year Air Temperature Maps
# moved the definition of this to where it's actually used.
#AT_Map_Path_A <- ("D:/PDKE/CCVD/CCVD INPUTS/air_temp/data_map_newer/") 

########## Monthly Climatology Air Temperature Maps
#AT_CLIM_Path_A <- ("D:/PDKE/CCVD/CCVD INPUTS/air_temp/air_temp_climatology/")
AT_CLIM_Path_A <- paste0(INPUTS_FOLDER, "/air_temp/air_temp_climatology/")

#################################################################################################################################

##########   DataFrames For Output

#Set up Group matrix for all climate Variables
#Climate For all Ranches
Cell.DataCL <- data.frame(matrix(ncol = 15, nrow = UNIT_C))
colnames(Cell.DataCL) <- c("Unit","Elev Min","Elev Mean", "Elev Max","RF","Tavg","Tmax","Tmin","RH","SM","KD","ET","CF","Island","Short Name")
Cell.DataCL[1] <- UNIT_N
Cell.DataCL[16] <- UNIT_Ns

#Climate For Individual Ranches
Cell.DataCLR <- data.frame(matrix(ncol = 16, nrow = 3))
colnames(Cell.DataCLR) <- c("Unit","Elev","RF","Tavg","Tmax","Tmin","RH","SM","KD","ET","CF","WS","DryRF","WetRF","Island","Short Name")

#Create a matrix For a site Comparision of rainfall
Cell.RF_Year <- data.frame(matrix(ncol = 15, nrow = UNIT_C))
colnames(Cell.RF_Year) <- c("Unit","Mean ELEV","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC","ANN")
Cell.RF_Year[1] <- UNIT_N

dpi = 300

u <- 1

# for (u in 1:1) {

#debug_print("Analysis UNIT ")  # Unit Analyzed
#debug_print(u)                 # Number In Loop
#debug_print(UNIT_N[u])         # Name Of Unit
#debug_print("Maps")

debug_print(paste("OUTPUTS_FOLDER: ", OUTPUTS_FOLDER))
debug_print(paste("UNIT_N[u]: ", UNIT_N[u]))
debug_print(paste("datetime_str: ", datetime_str))

#Create A Directory for Output
PATH <- paste0(OUTPUTS_FOLDER, UNIT_N[u], "_", datetime_str)
PATH_WITH_PROJECT_NAME <- paste0(PATH, "/", UNIT_N[u], "_") 
debug_print(paste("29, PATH:", PATH))
dir.create(PATH, showWarnings = TRUE, recursive = FALSE)
debug_print(paste("29, PATH_WITH_PROJECT_NAME:", PATH_WITH_PROJECT_NAME))

IS <- UNIT_I[1]       # Island
HALE<-st_as_sf(HALE)
SH <- HALE     # Shape File

##########   Export AOI shapefile

# png(paste0(OUTPUTS_FOLDER,UNIT_N[u],"/",UNIT_N[u]," Mokupuni.png"),width=5*dpi,height=4*dpi,res=dpi)
#
# writeOGR(obj = HALE,
#          dsn = paste0(OUTPUTS_FOLDER,UNIT_N[u],"/AOI/"),
#          layer = "AOI",
#          driver = "ESRI Shapefile")

##########   Get Coast for Island
debug_print("30, Get Coast for Island")

if(IS == "BI") {CoastM <- Coast_BI; Iname <- "Hawaii"}
if(IS == "MN") {CoastM <- Coast_MN; Iname <- "Maui"}
if(IS == "MO") {CoastM <- Coast_MO; Iname <- "Molokai"}
if(IS == "LA") {CoastM <- Coast_LA; Iname <- "Lanai"}
if(IS == "OA") {CoastM <- Coast_OA; Iname <- "Oahu"}
if(IS == "KA") {CoastM <- Coast_KA; Iname <- "Kauai"}
if(IS == "KO") {CoastM <- Coast_KO; Iname <- "Kahoolawe"}
if(IS == "ALL") {CoastM <- Coast_Crop_T; Iname <- "Hawaiian Islands"}

##########   GET Central Latitude and Longitude for shapefile

# Ensure CoastM is a POLYGON and extract coordinates properly
CoastM_df <- CoastM %>%
  st_cast("POLYGON") %>%  # Ensure proper geometry
  st_coordinates() %>%    # Extract coordinates
  as.data.frame() %>% 
  dplyr::rename(long = X, lat = Y)  # Rename columns for plotting

# Add group information for ggplot
CoastM_df$group <- CoastM_df$L1  

SH_sf <- st_as_sf(SH)
SH_sf$geometry

SH_centroid <- st_centroid(SH_sf$geometry)
LAT <- st_coordinates(SH_centroid)[,2]  # Y coordinate (latitude)
LONG <- st_coordinates(SH_centroid)[,1] # X coordinate (longitude)


##########   Make a Plot of Shapefile with Island Coast and LAT and LON
debug_print("31, Make a Plot of Shapefile with Island Coast and LAT and LON")

TIT <- paste0(UNIT_Ns[u], " (", Iname, ")")  # Title For Figure
dpi = 300

# find width of map
width <- abs(extent(CoastM)[1]) - abs(extent(CoastM)[2])

# convert width units (DD to miles)
# set latitude in radians
lat <- mean(c(extent(CoastM)[3], extent(CoastM)[4]))
lat.r <- DegToRad(lat)

# calculate width (miles) of 1 degree at latitude
width.d <- abs(cos(lat) * 69.172)

# calculate width (miles) of map
width.m <- round(width * width.d, 1)

# set scalebar width and convert to kilometers (to match coordinate system)
sb <- width.m * .25

if(sb<=3) {sb2<-1}
if(sb>3 && sb<=5) {sb2<-5}
if(sb>5 && sb<=10) {sb2<-10}
if(sb>10 && sb<=20) {sb2<-20}

sb2 <- sb2 * 1.60934

############# plot Mokupuni (island)
debug_print("32, Mokupuni")

dpi = 300
png(
  paste0(PATH_WITH_PROJECT_NAME, "Mokupuni.png"),
  width = 5 * dpi,
  height = 4 * dpi,
  res = dpi
)

ggplot() +
  geom_sf(data = Coast, fill = "aliceblue", color = "black", linewidth = 0.3) +  # Thicker outline for Coast
  geom_sf(data = CoastM, fill = "aliceblue", color = "brown", linewidth = 0.3) + # Thicker outline for CoastM
  geom_sf(data = HALE, fill = "aliceblue", color = "blue", linewidth = 0.5) +       # HALE in blue
  theme_void() +  # Remove all background (gridlines, axes, etc.)
  theme(legend.position = "none")  # Hide the legend if not needed

dev.off()

############# plot Moku
debug_print("33, Plot Moku")

# Ensure MOKU and HALE have the same CRS
HALE <- st_transform(HALE, st_crs(MOKU))

# Perform the intersection using sf
MIt <- st_intersection(MOKU, HALE)

MI2<-data.frame(MIt$objectid)
colnames(MI2)<-"ObjectID"

# subset moku for study area
MI3<-MOKU[MI2$ObjectID,]

# get subsetted moku as a dataframe for text labels
MId <- st_drop_geometry(MI3)
write.csv(data.frame(MId),
  paste0(PATH_WITH_PROJECT_NAME, "Moku.csv"),
  row.names = F
)

bbox <- st_bbox(MI3)  # Get bounding box

xm <- bbox["xmin"]
xma <- bbox["xmax"]
ym <- bbox["ymin"] * 0.999
yma <- bbox["ymax"] * 1.0001

# make map
png(
  paste0(PATH_WITH_PROJECT_NAME, "Moku.png"),
  width = 5 * dpi,
  height = 4 * dpi,
  res = dpi
)

ggplot() +
  geom_sf(data = Coast, col = "black", fill = "aliceblue", alpha = 0.5, linewidth = 1) +
  geom_sf(data = MI3, col = "brown", fill = NA, linewidth = 1) +
  geom_sf(data = HALE, col = "blue", fill = NA, linewidth = 1) +
  geom_text(data = MId, aes(label = moku, x = cent_long, y = cent_lat), 
            colour = "brown", size = 4, fontface = "bold") +
  coord_sf(xlim = c(xm, xma), ylim = c(ym, yma)) +
  theme_void()

dev.off()

########## Plot Ahupuaa

debug_print("34, Plot Ahupuaa")

# Ensure MOKU and HALE have the same CRS
HALE <- st_transform(HALE, st_crs(AHU))

# Perform the intersection using sf
AIt <- st_intersection(AHU, HALE)
AI2<-data.frame(AIt$objectid)
AI2
colnames(AI2)<-"ObjectID"

# subset ahupuaa for study area
AI3<-AHU[AI2$ObjectID,]

# clip polygons and calculate areas
AI3_sf<-st_as_sf(AI3)
HALE_sf<-st_as_sf(HALE)
HALE_area_m <- st_area(HALE_sf)  # Total area of HALE in square meters

AI3_clipped <- st_intersection(AI3_sf, HALE_sf)
AI3_clipped$area_m <- st_area(AI3_clipped)

# get subsetted ahupuaa as a dataframe for text labels
AId<-as.data.frame(AI3_clipped)

AI3_clipped_df <- st_drop_geometry(AI3_clipped)

write.csv(data.frame(AI3_clipped_df), paste0(PATH_WITH_PROJECT_NAME, "Ahupuaa.csv"),row.names = F)

bbox <- st_bbox(AI3)  # Get bounding box

xm <- bbox["xmin"]
xma <- bbox["xmax"]
ym <- bbox["ymin"] * 0.999
yma <- bbox["ymax"] * 1.0001

dpi=300

# make map
png(
  paste0(PATH_WITH_PROJECT_NAME, "Ahupuaa.png"),
  width = 5 * dpi,
  height = 4 * dpi,
  res = dpi
)

ggplot() +
  geom_sf(data = Coast, col = "black", fill = "aliceblue", alpha = 0.5, linewidth = 1) +
  geom_sf(data = AI3, col = "brown", fill = NA, linewidth = 1) +
  geom_sf(data = HALE, col = "blue", fill = NA, linewidth = 1) +
  geom_text(data = AId, aes(label = ahupuaa, x = cent_long, y = cent_lat), 
            colour = "black", size = 4, fontface = "bold") +
  coord_sf(xlim = c(xm, xma), ylim = c(ym, yma)) +
  theme_void()

dev.off()

##########  Fire Maps

debug_print("35, Fire Maps")

# Fire history map
FIRE_CropI <- FIRE_Shape_T
TITF <- "Fire Occurrence 1999-2022"

dpi = 300
png(
  paste0(PATH_WITH_PROJECT_NAME, "Fire.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

par(mfrow=c(1,1))

plot(st_geometry(CoastM), lwd = 1, lty = "solid", axes = FALSE, las = 2)
title(TITF, line = 0.5, cex.main = 1.5)
plot(st_geometry(SH_sf), add = TRUE, lwd = 2, col = "cyan")
plot(FIRE_CropI,add=TRUE,col="darkred")

# legend("topright", 
#        legend = c("Fire Occurrence", UNIT_Ns[[u]]), 
#        pch = c(20, 20), 
#        col = c("darkred", "cyan"),
#        bty = "n", 
#        horiz = FALSE, 
#        inset = c(-0.08, -0.08), 
#        cex = 0.8,          # Adjust text size
#        pt.cex = 2)         # Adjust symbol size (increase to make bigger)

legend("topright", 
       legend = c("Fire Occurrence", UNIT_Ns[[u]]), 
       pch = c(20, 20), 
       col = c("darkred", "cyan"),
       bty = "n", 
       horiz = FALSE, 
       inset = c(0.02, 0.02),  # always 2% inward from the corner
       cex = 0.8,
       pt.cex = 2)

dev.off()

# Fire risk map

### prepare hillshade for plotting
#clip to island
HSmoku <- mask(x = HS, mask = CoastM)

HSmoku <- crop(x = HS, y = extent(CoastM))
HSmoku[HSmoku == 0] <- NA

gplot_data <- function(x, maxpixels = ncell(x))  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x)))
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]],
                            by = c("value" = "ID"))
  }
  dat
}

HSd <- gplot_data(HSmoku)

bbox <- st_bbox(HALE)  # Get bounding box

xm <- bbox["xmin"]
xma <- bbox["xmax"]
ym <- bbox["ymin"] * 0.999
yma <- bbox["ymax"] * 1.0001

########### Elevation Extract and Map

debug_print("36, Elevation")

ELEV_MaskI <- mask(x = ELEV, mask = CoastM)
ELEV_CropI <- crop(x = ELEV_MaskI, y = extent(CoastM))

ELEV_Mask <- mask(x = ELEV, mask = HALE)

# Crop the masked raster to the extent of the polygon UNIT_X[[u]]
ELEV_Crop <- crop(x = ELEV_Mask, y = extent(HALE))

# Calculate statistics (mean, max, and min) of the cropped raster
ELEV_Mean <- round(cellStats(ELEV_Crop, 'mean'), 1)
ELEV_Max <- round(cellStats(ELEV_Crop, 'max'), 1)
ELEV_Min <- round(cellStats(ELEV_Crop, 'min'), 1)

Cell.DataCL[u, 2:4] <- c(ELEV_Min, ELEV_Mean, ELEV_Max)
Cell.DataCL[u, 14] <- Iname
Cell.DataCLR[1:3, 1] <- c("Mean", "Max", "Min")
Cell.DataCLR[1:3, 2] <- c(ELEV_Mean, ELEV_Max, ELEV_Min)
Cell.DataCLR[1, 15] <- Iname
Cell.DataCLR[1, 16] <- UNIT_Ns[u]

Cell.DataCL
Cell.DataCLR

#BI_brksRH<-round(seq(RHLO, RHUP, length = 9),0)
colfuncEL <- colorRampPalette(brewer.pal(9, "PuBuGn"))(100)

# Load your raster (ELEV_CropI) and convert to a data frame
elev_df <- as.data.frame(ELEV_CropI, xy = TRUE, na.rm = TRUE)

# Convert the CoastM and HALE to sf objects if they aren't already
CoastM_sf <- st_as_sf(CoastM)
HALE_sf <- st_as_sf(HALE)


# Get map bounds (so the arrow scales with your map)
xrange <- range(elev_df$x, na.rm = TRUE)
yrange <- range(elev_df$y, na.rm = TRUE)

# Define arrow placement near top-right corner
x_arrow <- xrange[2] - 0.05 * diff(xrange)   # 5% from right edge
y_arrow_bottom <- yrange[1] + 0.01 * diff(yrange)  # base of arrow (15% up)
y_arrow_top <- yrange[1] + 0.12 * diff(yrange)     # arrow tip (30% up)
y_N_label <- y_arrow_top + 0.04 * diff(yrange)     # N label slightly above arrow

debug_print(paste0("ELMap: ", PATH_WITH_PROJECT_NAME, "ELMap.png"))
png(
  paste0(PATH_WITH_PROJECT_NAME, "ELMap.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

ggplot() +
  # Plot the elevation raster
  geom_raster(data = elev_df, aes(x = x, y = y, fill = ned_dem)) +
  scale_fill_gradientn(colors = colfuncEL, name = "Elev. (ft.)") +  # Use color ramp for elevation
  # Plot the CoastM polygon layer
  geom_sf(data = CoastM_sf, fill = NA, color = "black", size = 0.8) +
  geom_sf(data = HALE_sf, fill = NA, color = "darkred", linewidth = 1) +
  coord_sf() +  # Use coord_sf() to handle spatial data
  
  # Add the N label
  annotate("text", x = x_arrow, y = y_N_label, label = "N", size = 4, fontface = "bold") +
  
  # Add the arrow line
  annotate("segment",
           x = x_arrow, xend = x_arrow,
           y = y_arrow_bottom, yend = y_arrow_top,
           arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
           size = 1, color = "black") +
  
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(title = paste0("Elevation ", Iname, " (", ELUnit2, ")"))

dev.off()

########## Landcover Extraction and Map
debug_print("37, Landcover")

LC_MaskI <- mask(x = LC, mask = CoastM)
LC_CropI <- crop(x = LC_MaskI, y = extent(CoastM))

#Mask For Geography
LC_Mask <- mask(x = LC, mask = HALE)
LC_Crop <- crop(x = LC_Mask, y = extent(HALE))

#Convert raster to dataframe for stats
LC_Crop2 <- as.data.frame(LC_Crop, xy = T)

# remove LC "0" (no data)
LC_Crop3 <- LC_Crop2[which(LC_Crop2$LCMAP_HI_2020_V10_LCPRI_crs != 0), ]

#Count frequencies
LC_ct <- as.data.frame(table(LC_Crop3$LCMAP_HI_2020_V10_LCPRI_crs))
colnames(LC_ct) <- c("LC", "count")
LC_ct2 <- LC_ct[order(-LC_ct$count), ]

#convert into area (1 pixel = 900 square meters = 0.222395 acres)
LC_ct2$acres <- round(LC_ct2$count * 0.222395, 1)

#Calculate percentages
all <- sum(LC_ct2$count)

LC_ct2$pct <- round((LC_ct2$count / all) * 100, 0)

for (x in 1:nrow(LC_ct2)) {
  if(LC_ct2[x,]$pct == 0) {LC_ct2[x,]$pct<-"<1"}}

LC_ct2$pct <- paste0(LC_ct2$pct, "%")
LC_ct2

#define class names and colors
LC_ct2$class_name <- as.character(NA)
LC_ct2$color <- as.character(NA)

for (x in 1:nrow(LC_ct2)) {
  
  if(LC_ct2[x,]$LC == 1) {LC_ct2[x,]$class_name <- c("Developed")}
  if(LC_ct2[x,]$LC == 2) {LC_ct2[x,]$class_name <- c("Cropland")}
  if(LC_ct2[x,]$LC == 3) {LC_ct2[x,]$class_name <- c("Grass/Shrub")}
  if(LC_ct2[x,]$LC == 4) {LC_ct2[x,]$class_name <- c("Tree Cover")}
  if(LC_ct2[x,]$LC == 5) {LC_ct2[x,]$class_name <- c("Water")}
  if(LC_ct2[x,]$LC == 6) {LC_ct2[x,]$class_name <- c("Wetland")}
  if(LC_ct2[x,]$LC == 8) {LC_ct2[x,]$class_name <- c("Barren")}
  
}

LC_ct2

# write to table
write.csv(LC_ct2,
  paste0(PATH_WITH_PROJECT_NAME, "Landcover.csv"),
  row.names = F)

# set class names as leveled factors
LC_ct2$class_name<-factor(LC_ct2$class_name, 
  levels=c("Developed","Cropland","Grass/Shrub","Tree Cover",
    "Water","Wetland","Barren"))

# set first thru third landcover classes
LC1 <- paste0(LC_ct2[1, ]$class_name, "(", LC_ct2[1, ]$pct, ")")
LC2 <- paste0(LC_ct2[2, ]$class_name, "(", LC_ct2[2, ]$pct, ")")
LC3 <- paste0(LC_ct2[3, ]$class_name, "(", LC_ct2[3, ]$pct, ")")

### Plot Bar Graph
# set bargraph ylim
ylim <- max(LC_ct2$acres) + (max(LC_ct2$acres) * 0.15)

png(
  paste0(PATH_WITH_PROJECT_NAME, "LC_barchart.png"),
  width = 5 * dpi,
  height = 3 * dpi,
  res = dpi
)

print(ggplot(LC_ct2, aes(x=class_name, y=acres, fill=class_name)) +
    geom_col(width=0.9, position = position_dodge(0.7), color="black") +
    theme_bw() +
    ylim(0,ylim)+
    scale_fill_manual(values = c("Developed" = "black",
      "Cropland" = "darkgoldenrod1",
      "Grass/Shrub" = "darkolivegreen1",
      "Tree Cover" = "darkgreen",
      "Water" = "blue",
      "Wetland" = "aquamarine",
      "Barren" = "chocolate4")) +
    geom_text(aes(label = pct), vjust = -0.8, size=4.5) +
    theme(axis.text.x = element_text(angle=35, vjust = .65, size=14)) +
    labs(y="Area (acres)", x="", ) +
    theme(legend.position="none") +
    theme(axis.title.y=element_text(size=15), 
      axis.text.y=element_text(size=14)))

dev.off()

# export map

png(
  paste0(PATH_WITH_PROJECT_NAME, "LCMap.png"),
  width = 7 * dpi,
  height = 5 * dpi,
  res = dpi
)

# Define land cover values and their corresponding colors
landcover_values <- c(1, 2, 3, 4, 5, 6, 8)
landcover_colors <- c("black", "darkgoldenrod1", "darkolivegreen1", 
                      "darkgreen", "blue", "aquamarine", "chocolate4")

# Ensure correct color assignment by explicitly setting the breaks and labels
plot(LC_CropI, col=landcover_colors, breaks=c(0, landcover_values), 
     legend=FALSE, xaxt="n", yaxt="n", bty="n")
plot(st_geometry(CoastM), add=TRUE, border="black", col=NA, lwd=1)
plot(st_geometry(HALE), add=TRUE, border="red", col=NA, lwd=2)
title(main=paste0("Landcover", " ", Iname))

dev.off()


################################################################################
########## Water Sources
debug_print("38, Water Sources")

### Aquifers

# Ensure MOKU and AQU have the same CRS
HALE <- st_transform(HALE, st_crs(AQU))

# Perform the intersection using sf
AIt <- st_intersection(AQU, HALE[1])
AI2<-data.frame(AIt$objectid)
colnames(AI2)<-"ObjectID"

AIi <- st_intersection(AQU, CoastM)
AIi <- AIi[35]
AIi <- st_cast(AIi, "POLYGON")

# subset aquifers for study area
AI3<-AI2[AI2$ObjectID,]

# subset aquifers for study area
AQ3 <- AQU[AI2$ObjectID, ]

# clip polygons and calculate centroid coordinates
AQ3_sf<-st_as_sf(AIt)
AQ3_clipped<-AQ3_sf

AQ3_clipped <- st_make_valid(AQ3_clipped)

AQ3_centroids <- st_centroid(AQ3_clipped)

# Extract latitude and longitude from centroids
centroid_coords <- st_coordinates(AQ3_centroids)
AQ3_clipped <- AQ3_clipped %>%
  mutate(
    cent_long = centroid_coords[, 1],
    cent_lat = centroid_coords[, 2]
  )

# I added the following line because it was getting an error "Loop 0 is not valid: Edge 257 is degenerate (duplicate vertex)"
AQ3_validated <- st_make_valid(AQ3_clipped) # Fixes self-intersections and other issues
AQ3_centroids <- st_centroid(AQ3_validated)

# subset for hydrology
AQb<-AIt[which(AIt$typea1 == 1),]
AQh<-AIt[which(AIt$typea1 == 2),]  

# make data frame and only keep relevant columns
# AQd<-as.data.frame(AQ3_clipped)
AQd<-as.data.frame(AIt)
AQd<-as.data.frame(AQ3_clipped)
AQd<-subset(AQd, select = c("objectid","doh_aquife","typea1","typea3","stata1","strata3","cent_lat","cent_long"))

# translate type and status codes to text
AQd$Hydrology <- NA
AQd$Geology <- NA
AQd$Salinity <- NA
AQd$Use <- NA

for (i in 1:nrow(AQd)) {

  if(AQd[i,]$typea1 == 1) {AQd[i,]$Hydrology<-c("Basal")}
  if(AQd[i,]$typea1 == 2) {AQd[i,]$Hydrology<-c("High Level")}
  if(AQd[i,]$typea3 == 1) {AQd[i,]$Geology<-c("Flank")}
  if(AQd[i,]$typea3 == 2) {AQd[i,]$Geology<-c("Dike")}
  if(AQd[i,]$typea3 == 3) {AQd[i,]$Geology<-c("Flank/Dike")}
  if(AQd[i,]$typea3 == 4) {AQd[i,]$Geology<-c("Perched")}
  if(AQd[i,]$typea3 == 5) {AQd[i,]$Geology<-c("Dike/Perched")}
  if(AQd[i,]$typea3 == 6) {AQd[i,]$Geology<-c("Sedimentary")}
  if(AQd[i,]$stata1 == 1) {AQd[i,]$Use<-c("Currently used")}
  if(AQd[i,]$stata1 == 2) {AQd[i,]$Use<-c("Potential use")}
  if(AQd[i,]$stata1 == 3) {AQd[i,]$Use<-c("No potential use")}
  if(AQd[i,]$strata3 == 1) {AQd[i,]$Salinity<-c("Fresh")}
  if(AQd[i,]$strata3 == 2) {AQd[i,]$Salinity<-c("Low")}
  if(AQd[i,]$strata3 == 3) {AQd[i,]$Salinity<-c("Moderate")}
  if(AQd[i,]$strata3 == 4) {AQd[i,]$Salinity<-c("High")}
  if(AQd[i,]$strata3 == 5) {AQd[i,]$Salinity<-c("Seawater")}
  
}

# trim it down to the final table
AQe<-subset(AQd, select = c("doh_aquife", "Hydrology", "Geology", "Salinity", "Use"))
colnames(AQe)[1] = "DOH Aquifer"
write.csv(
  data.frame(AQe),
  paste0(PATH_WITH_PROJECT_NAME, "Aquifer.csv"),
  row.names = F
)

bbox <- st_bbox(HALE)  # Get bounding box
xm <- bbox["xmin"] * 1.0005
xma <- bbox["xmax"] *0.9995
ym <- bbox["ymin"] * 0.999
yma <- bbox["ymax"] * 1.001

# make map
png(
  paste0(PATH_WITH_PROJECT_NAME, "Aquifers.png"),
  width = 4 * dpi,
  height = 4 * dpi,
  res = dpi
)
par(mfrow = c(1, 1))

ggplot() +
    # Plot AIi as the bottom layer with conditional fill colors
    geom_sf(data = AIi, aes(fill = factor(typea1)), color = "black", linewidth = 1) +
    scale_fill_manual(values = c("1" = "lightblue", "2" = "grey")) +
    # Plot CoastM with no fill and black outline
    geom_sf(data = CoastM, fill = NA, color = "black") +
    # Plot AQb with lightblue fill and black outline
    geom_sf(data = AQb, fill = "lightblue", color = "black") +
    # Plot AQh with grey fill and black outline
    geom_sf(data = AQh, fill = "grey", color = "black") +
    # Plot HALE with no fill and red outline
    geom_sf(data = HALE, fill = NA, color = "red", linewidth = 1) +
    # Add text labels from AQd using cent_lat and cent_long
    geom_text(data = AQd, aes(x = cent_long, y = cent_lat, label = doh_aquife),
              size = 3, fontface = "bold") +
    # Set plot limits
    coord_sf(xlim = c(xm, xma), ylim = c(ym, yma)) +
    # Remove grid/ticks and legend
    theme_void() +
    theme(legend.position = "none")

dev.off()

################################################################################
### Streams, hydrologic features

debug_print("39, Streams")

bbox <- st_bbox(HALE)  # Get bounding box

xm <- as.numeric(bbox["xmin"]) * 1.00001
xma <- as.numeric(bbox["xmax"]) * 0.99999
ym <- as.numeric(bbox["ymin"]) * 0.9999
yma <- as.numeric(bbox["ymax"]) * 1.0001

# Create an sf bounding box object
bbox_poly <- st_as_sfc(st_bbox(c(xmin = xm, xmax = xma, ymin = ym, ymax = yma), crs = st_crs(HALE)))

bbox_extent <- extent(as.numeric(xm), as.numeric(xma), as.numeric(ym), as.numeric(yma))
HS_crop <- crop(HS, bbox_extent)

# subset for feature types
STRM <- st_transform(STRM, st_crs(HS))
STRMa<-STRM[which(STRM$fcode == 55800),]
STRMp<-STRM[which(STRM$fcode == 46006),]
STRMi<-STRM[which(STRM$fcode == 46003),]
STRMc<-STRM[which(STRM$Feature_Ty == "CANAL/DITCH"),]
STRMpi<-STRM[which(STRM$Feature_Ty == "PIPELINE"),]
STRMco<-STRM[which(STRM$Feature_Ty == "CONNECTOR"),]

# Clip each stream layer to the HS extent
STRMa_clipped <- st_intersection(STRMa, bbox_poly)
STRMp_clipped <- st_intersection(STRMp, bbox_poly)
STRMi_clipped <- st_intersection(STRMi, bbox_poly)
STRMc_clipped <- st_intersection(STRMc, bbox_poly)
STRMpi_clipped <- st_intersection(STRMpi, bbox_poly)

# make map
png(
  paste0(PATH_WITH_PROJECT_NAME, "Streams.png"),
  width = 4 * dpi,
  height = 4 * dpi,
  res = dpi
)

par(mar = c(0, 0, 0, 0))  # No margins

plot(HS_crop, col = gray.colors(255), legend = FALSE, 
     axes = FALSE, box = FALSE)

plot(st_geometry(STRMa_clipped), col = "green", lwd = 2, add = TRUE)
plot(st_geometry(STRMp_clipped), col = "blue", lwd = 2, add = TRUE)
plot(st_geometry(STRMi_clipped), col = "lightblue1", lwd = 2, add = TRUE)
plot(st_geometry(STRMc_clipped), col = "orange", lwd = 2, add = TRUE)
plot(st_geometry(STRMpi_clipped), col = "black", lwd = 2, add = TRUE)
plot(st_geometry(STRMc_clipped), col = "orange", lwd = 2, add = TRUE)

plot(st_geometry(HALE), add = TRUE, border = "red", lwd = 2)  # Add HALE boundary


dev.off()

# clip features to study area
# Perform the intersection using sf
ST <- st_intersection(STRM, HALE[1])

if(nrow(ST) == 0) {
  ht<-data.frame()
}

if(nrow(ST) != 0) {
  
  # get all hydrologic types present
  ht<-unique(ST$Feature_Ty)
  ht
  for (i in 1:length(ht)){
    if(ht[i] == "ARTIFICIAL PATH"){ht[i]<-"Managed Waterway"}
    if(ht[i] == "CANAL/DITCH"){ht[i]<-"Canal/Ditch"}
    if(ht[i] == "STREAM/RIVER"){ht[i]<-"Stream"}
    if(ht[i] == "PIPELINE"){ht[i]<-"Pipeline"}
    if(ht[i] == "CONNECTOR"){ht[i]<-"Pipeline"}
  }
}

# remove duplicate types before final output
ht <- unique(ht)
ht

write.csv(
  data.frame(ht),
  paste0(PATH_WITH_PROJECT_NAME, "Hydro_features.csv"),
  row.names = F
)

########## Rainfall Stations near AOI
debug_print("40, Rain Station Locations")

# load csv of station networks and website links
sn <- read.csv(paste0(INPUTS_FOLDER, "rain_stations_links.csv"))

# load stations point shapefile, tranform points and polygons to planar projection
rain_stations_file <- paste0(INPUTS_FOLDER, "rain_stations.shp")
debug_print(rain_stations_file)
sl <- st_read(rain_stations_file)

# Convert SpatialPointsDataFrame (sp object) to sf object
sl_sf <- st_as_sf(sl)

# Ensure both layers have the same CRS (preferably projected for accuracy)
projected_crs <- "+proj=utm +zone=4 +datum=NAD83 +units=m +no_defs"
sl_proj <- st_transform(sl_sf, crs = projected_crs)
HALE_proj <- st_transform(HALE, crs = projected_crs)

# Compute distances from each point in sl_proj to the HALE_proj polygon
distances <- st_distance(sl_proj, HALE_proj)

# Convert to a data frame
sd <- data.frame(SKN = sl_proj$SKN, distance = as.numeric(distances))
head(sd, 20)

# sort by distance
sd<-sd[order(sd$distance),]
head(sd)

# convert station points shp into dataframe
sld<-as.data.frame(sl_proj)
head(sld)

# make table of 3 closest stations and join network links
s1<-sld[which(sld$SKN == sd[1,]$SKN),]
s2<-sld[which(sld$SKN == sd[2,]$SKN),]
s3<-sld[which(sld$SKN == sd[3,]$SKN),]
st2<-rbind(s1,s2,s3)
st2
st<-st2[c("Station_Na","Network")]
st<-join(st,sn)
colnames(st)<-c("Station Name","Network","Website")

# export table
write.csv(st,
  paste0(PATH_WITH_PROJECT_NAME, "rain_stations.csv"),
  row.names = F)

# ### Make map of island with station and AOI
HS_masked <- mask(HS, CoastM)

bbox <- st_bbox(CoastM)  # Get bounding box

xm <- as.numeric(bbox["xmin"])
xma <- as.numeric(bbox["xmax"])
ym <- as.numeric(bbox["ymin"])
yma <- as.numeric(bbox["ymax"])

# make map
png(
  paste0(PATH_WITH_PROJECT_NAME, "rf_stations.png"),
  width = 5 * dpi,
  height = 3.5 * dpi,
  res = dpi
)

# Plot
par(mar = c(0, 2, 0, 0))  # bottom, left, top, right

plot(HS_masked, col = gray.colors(255), legend = FALSE, 
     xlim = c(xm, xma), ylim = c(ym, yma), 
     axes = FALSE, box = FALSE) 
points(sld$LON, sld$LAT, col = "black", pch = 16, cex = 0.7)
points(st2$LON, st2$LAT, col = "orange", pch = 16, cex = 0.7)
plot(HALE$geometry, add = TRUE, border = "blue", lwd = 2)

dev.off()

########## Mean Climate Extract
debug_print("41, Mean Climate")

#Create a climate matrix
Cell.matrix5 <- matrix(ncol = 14, nrow = 9)
Cell.CL_Year <- data.frame(Cell.matrix5)
colnames(Cell.CL_Year) <- c("Variable","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC","ANN")
Cell.CL_Year[1, 1] <-  paste0("RF [", RFUnit2, "]")
Cell.CL_Year[2, 1] <- paste0("Min TA [", TUnit, "]")
Cell.CL_Year[3, 1] <-  paste0("Mean TA [", TUnit, "]")
Cell.CL_Year[4, 1] <-  paste0("Max TA [", TUnit, "]")
Cell.CL_Year[5, 1] <- "RH [%]"
Cell.CL_Year[6, 1] <- "CF [%]"
Cell.CL_Year[7, 1] <- paste0("ET [", RFUnit2, "]")
Cell.CL_Year[8, 1] <- "SM [%]"
Cell.CL_Year[9, 1] <- "S [W m/2]"
Cell.CL_Year

##########   Solar Radiation
debug_print("42, Solar Radiation")

Jan <- raster(Mean_CLIM[58])
Jan_Mask <- mask(x = Jan, mask = HALE)
Jan_CropKD <- crop(x = Jan_Mask, y = extent(HALE))
JanMKD   <- round(cellStats(Jan_CropKD, 'mean'), 0)
JanMKDx   <- round(cellStats(Jan_CropKD, 'max'), 0)
JanMKDn   <- round(cellStats(Jan_CropKD, 'min'), 0)

Feb <- raster(Mean_CLIM[57])
Feb_Mask <- mask(x = Feb, mask = HALE)
Feb_CropKD <- crop(x = Feb_Mask, y = extent(HALE))
FebMKD   <- round(cellStats(Feb_CropKD, 'mean'), 0)
FebMKDx   <- round(cellStats(Feb_CropKD, 'max'), 0)
FebMKDn   <- round(cellStats(Feb_CropKD, 'min'), 0)

Mar <- raster(Mean_CLIM[61])
Mar_Mask <- mask(x = Mar, mask = HALE)
Mar_CropKD <- crop(x = Mar_Mask, y = extent(HALE))
MarMKD   <- round(cellStats(Mar_CropKD, 'mean'), 0)
MarMKDx   <- round(cellStats(Mar_CropKD, 'max'), 0)
MarMKDn   <- round(cellStats(Mar_CropKD, 'min'), 0)

Apr <- raster(Mean_CLIM[54])
Apr_Mask <- mask(x = Apr, mask = HALE)
Apr_CropKD <- crop(x = Apr_Mask, y = extent(HALE))
AprMKD   <- round(cellStats(Apr_CropKD, 'mean'), 0)
AprMKDx   <- round(cellStats(Apr_CropKD, 'max'), 0)
AprMKDn   <- round(cellStats(Apr_CropKD, 'min'), 0)

May <- raster(Mean_CLIM[62])
May_Mask <- mask(x = May, mask = HALE)
May_CropKD <- crop(x = May_Mask, y = extent(HALE))
MayMKD   <- round(cellStats(May_CropKD, 'mean'), 0)
MayMKDx   <- round(cellStats(May_CropKD, 'max'), 0)
MayMKDn   <- round(cellStats(May_CropKD, 'min'), 0)

Jun <- raster(Mean_CLIM[60])
Jun_Mask <- mask(x = Jun, mask = HALE)
Jun_CropKD <- crop(x = Jun_Mask, y = extent(HALE))
JunMKD   <- round(cellStats(Jun_CropKD, 'mean'), 0)
JunMKDx   <- round(cellStats(Jun_CropKD, 'max'), 0)
JunMKDn   <- round(cellStats(Jun_CropKD, 'min'), 0)

Jul <- raster(Mean_CLIM[59])
Jul_Mask <- mask(x = Jul, mask = HALE)
Jul_CropKD <- crop(x = Jul_Mask, y = extent(HALE))
JulMKD   <- round(cellStats(Jul_CropKD, 'mean'), 0)
JulMKDx   <- round(cellStats(Jul_CropKD, 'max'), 0)
JulMKDn   <- round(cellStats(Jul_CropKD, 'min'), 0)

Aug <- raster(Mean_CLIM[55])
Aug_Mask <- mask(x = Aug, mask = HALE)
Aug_CropKD <- crop(x = Aug_Mask, y = extent(HALE))
AugMKD   <- round(cellStats(Aug_CropKD, 'mean'), 0)
AugMKDx   <- round(cellStats(Aug_CropKD, 'max'), 0)
AugMKDn   <- round(cellStats(Aug_CropKD, 'min'), 0)

Sep <- raster(Mean_CLIM[65])
Sep_Mask <- mask(x = Sep, mask = HALE)
Sep_CropKD <- crop(x = Sep_Mask, y = extent(HALE))
SepMKD   <- round(cellStats(Sep_CropKD, 'mean'), 0)
SepMKDx   <- round(cellStats(Sep_CropKD, 'max'), 0)
SepMKDn   <- round(cellStats(Sep_CropKD, 'min'), 0)

Oct <- raster(Mean_CLIM[64])
Oct_Mask <- mask(x = Oct, mask = HALE)
Oct_CropKD <- crop(x = Oct_Mask, y = extent(HALE))
OctMKD   <- round(cellStats(Oct_CropKD, 'mean'), 0)
OctMKDx   <- round(cellStats(Oct_CropKD, 'max'), 0)
OctMKDn   <- round(cellStats(Oct_CropKD, 'min'), 0)

Nov <- raster(Mean_CLIM[63])
Nov_Mask <- mask(x = Nov , mask = HALE)
Nov_CropKD <- crop(x = Nov_Mask, y = extent(HALE))
NovMKD   <- round(cellStats(Nov_CropKD, 'mean'), 0)
NovMKDx   <- round(cellStats(Nov_CropKD, 'max'), 0)
NovMKDn   <- round(cellStats(Nov_CropKD, 'min'), 0)

Dec <- raster(Mean_CLIM[56])
Dec_Mask <- mask(x = Dec, mask = HALE)
Dec_CropKD <- crop(x = Dec_Mask, y = extent(HALE))
DecMKD   <- round(cellStats(Dec_CropKD, 'mean'), 0)
DecMKDx   <- round(cellStats(Dec_CropKD, 'max'), 0)
DecMKDn   <- round(cellStats(Dec_CropKD, 'min'), 0)

Ann <- raster(Mean_CLIM[53])
Ann_Mask <- mask(x = Ann, mask = HALE)
Ann_CropKD <- crop(x = Ann_Mask, y = extent(HALE))
AnnMKD   <- round(cellStats(Ann_CropKD, 'mean'), 0)
AnnMKDx   <- round(cellStats(Ann_CropKD, 'max'), 0)
AnnMKDn   <- round(cellStats(Ann_CropKD, 'min'), 0)

MEANKD2 <- c(JanMKD,FebMKD,MarMKD,AprMKD,MayMKD,JunMKD,JulMKD,AugMKD,SepMKD,OctMKD,NovMKD,DecMKD,AnnMKD)
MEANKD3 <- MEANKD2 * 100

#GKD thresholds for scale Bar
KDUP <- max(JanMKDx,FebMKDx,MarMKDx,AprMKDx,MayMKDx,JunMKDx,JulMKDx,AugMKDx,SepMKDx,OctMKDx,NovMKDx,DecMKDx,AnnMKDx)
KDLO <- min(JanMKDn,FebMKDn,MarMKDn,AprMKDn,MayMKDn,JunMKDn,JulMKDn,AugMKDn,SepMKDn,OctMKDn,NovMKDn,DecMKDn,AnnMKDn)

Cell.CL_Year[9, 2:14] <- MEANKD2
Cell.CL_Year

##########   Soil Moisture
debug_print("43, Soil Moisture")

Jan <- raster(Mean_CLIM[45])
Jan_Mask <- mask(x = Jan, mask = HALE)
Jan_CropSM <- crop(x = Jan_Mask, y = extent(HALE))
JanMSM   <- round(cellStats(Jan_CropSM, 'mean'), 2)
JanMSMx   <- round(cellStats(Jan_CropSM, 'max'), 2)
JanMSMn   <- round(cellStats(Jan_CropSM, 'min'), 2)

Feb <- raster(Mean_CLIM[44])
Feb_Mask <- mask(x = Feb, mask = HALE)
Feb_CropSM <- crop(x = Feb_Mask, y = extent(HALE))
FebMSM   <- round(cellStats(Feb_CropSM, 'mean'), 2)
FebMSMx   <- round(cellStats(Feb_CropSM, 'max'), 2)
FebMSMn   <- round(cellStats(Feb_CropSM, 'min'), 2)

Mar <- raster(Mean_CLIM[48])
Mar_Mask <- mask(x = Mar, mask = HALE)
Mar_CropSM <- crop(x = Mar_Mask, y = extent(HALE))
MarMSM   <- round(cellStats(Mar_CropSM, 'mean'), 2)
MarMSMx   <- round(cellStats(Mar_CropSM, 'max'), 2)
MarMSMn   <- round(cellStats(Mar_CropSM, 'min'), 2)

Apr <- raster(Mean_CLIM[41])
Apr_Mask <- mask(x = Apr, mask = HALE)
Apr_CropSM <- crop(x = Apr_Mask, y = extent(HALE))
AprMSM   <- round(cellStats(Apr_CropSM, 'mean'), 2)
AprMSMx   <- round(cellStats(Apr_CropSM, 'max'), 2)
AprMSMn   <- round(cellStats(Apr_CropSM, 'min'), 2)

May <- raster(Mean_CLIM[49])
May_Mask <- mask(x = May, mask = HALE)
May_CropSM <- crop(x = May_Mask, y = extent(HALE))
MayMSM   <- round(cellStats(May_CropSM, 'mean'), 2)
MayMSMx   <- round(cellStats(May_CropSM, 'max'), 2)
MayMSMn   <- round(cellStats(May_CropSM, 'min'), 2)

Jun <- raster(Mean_CLIM[47])
Jun_Mask <- mask(x = Jun, mask = HALE)
Jun_CropSM <- crop(x = Jun_Mask, y = extent(HALE))
JunMSM   <- round(cellStats(Jun_CropSM, 'mean'), 2)
JunMSMx   <- round(cellStats(Jun_CropSM, 'max'), 2)
JunMSMn   <- round(cellStats(Jun_CropSM, 'min'), 2)

Jul <- raster(Mean_CLIM[46])
Jul_Mask <- mask(x = Jul, mask = HALE)
Jul_CropSM <- crop(x = Jul_Mask, y = extent(HALE))
JulMSM   <- round(cellStats(Jul_CropSM, 'mean'), 2)
JulMSMx   <- round(cellStats(Jul_CropSM, 'max'), 2)
JulMSMn   <- round(cellStats(Jul_CropSM, 'min'), 2)

Aug <- raster(Mean_CLIM[42])
Aug_Mask <- mask(x = Aug, mask = HALE)
Aug_CropSM <- crop(x = Aug_Mask, y = extent(HALE))
AugMSM   <- round(cellStats(Aug_CropSM, 'mean'), 2)
AugMSMx   <- round(cellStats(Aug_CropSM, 'max'), 2)
AugMSMn   <- round(cellStats(Aug_CropSM, 'min'), 2)

Sep <- raster(Mean_CLIM[52])
Sep_Mask <- mask(x = Sep, mask = HALE)
Sep_CropSM <- crop(x = Sep_Mask, y = extent(HALE))
SepMSM   <- round(cellStats(Sep_CropSM, 'mean'), 2)
SepMSMx   <- round(cellStats(Sep_CropSM, 'max'), 2)
SepMSMn   <- round(cellStats(Sep_CropSM, 'min'), 2)

Oct <- raster(Mean_CLIM[51])
Oct_Mask <- mask(x = Oct, mask = HALE)
Oct_CropSM <- crop(x = Oct_Mask, y = extent(HALE))
OctMSM   <- round(cellStats(Oct_CropSM, 'mean'), 2)
OctMSMx   <- round(cellStats(Oct_CropSM, 'max'), 2)
OctMSMn   <- round(cellStats(Oct_CropSM, 'min'), 2)

Nov <- raster(Mean_CLIM[50])
Nov_Mask <- mask(x = Nov , mask = HALE)
Nov_CropSM <- crop(x = Nov_Mask, y = extent(HALE))
NovMSM   <- round(cellStats(Nov_CropSM, 'mean'), 2)
NovMSMx   <- round(cellStats(Nov_CropSM, 'max'), 2)
NovMSMn   <- round(cellStats(Nov_CropSM, 'min'), 2)

Dec <- raster(Mean_CLIM[43])
Dec_Mask <- mask(x = Dec, mask = HALE)
Dec_CropSM <- crop(x = Dec_Mask, y = extent(HALE))
DecMSM   <- round(cellStats(Dec_CropSM, 'mean'), 2)
DecMSMx   <- round(cellStats(Dec_CropSM, 'max'), 2)
DecMSMn   <- round(cellStats(Dec_CropSM, 'min'), 2)

Ann <- raster(Mean_CLIM[40])
Ann_Mask <- mask(x = Ann, mask = HALE)
Ann_CropSM <- crop(x = Ann_Mask, y = extent(HALE))
AnnMSM   <- round(cellStats(Ann_CropSM, 'mean'), 2)
AnnMSMx   <- round(cellStats(Ann_CropSM, 'max'), 2)
AnnMSMn   <- round(cellStats(Ann_CropSM, 'min'), 2)

MEANSM2 <- c(JanMSM,FebMSM,MarMSM,AprMSM,MayMSM,JunMSM,JulMSM,AugMSM,SepMSM,OctMSM,NovMSM,DecMSM,AnnMSM)
MEANSM3 <- MEANSM2 * 100

#GSM thresholds for scale Bar
SMUP <- max(JanMSMx,FebMSMx,MarMSMx,AprMSMx,MayMSMx,JunMSMx,JulMSMx,AugMSMx,SepMSMx,OctMSMx,NovMSMx,DecMSMx,AnnMSMx)
SMLo <- min(JanMSMn,FebMSMn,MarMSMn,AprMSMn,MayMSMn,JunMSMn,JulMSMn,AugMSMn,SepMSMn,OctMSMn,NovMSMn,DecMSMn,AnnMSMn)

Cell.CL_Year[8, 2:14] <- MEANSM3
Cell.CL_Year

##########   EVAPOTRANSPIRATION
debug_print("44, EVAPOTRANSPIRATION")

Jan <- raster(Mean_CLIM[19])
Jan_Mask <- mask(x = Jan, mask = HALE)
Jan_CropET <- crop(x = Jan_Mask, y = extent(HALE))
if (RFUnit == " in") { Jan_CropET <- Jan_CropET  * 0.0393701 }
JanMET   <- round(cellStats(Jan_CropET, 'mean'), 0)
JanMETx   <- round(cellStats(Jan_CropET, 'max'), 0)
JanMETn   <- round(cellStats(Jan_CropET, 'min'), 0)

Feb <- raster(Mean_CLIM[18])
Feb_Mask <- mask(x = Feb, mask = HALE)
Feb_CropET <- crop(x = Feb_Mask, y = extent(HALE))
if (RFUnit == " in") { Feb_CropET <- Feb_CropET  * 0.0393701 }
FebMET   <- round(cellStats(Feb_CropET, 'mean'), 0)
FebMETx   <- round(cellStats(Feb_CropET, 'max'), 0)
FebMETn   <- round(cellStats(Feb_CropET, 'min'), 0)

Mar <- raster(Mean_CLIM[22])
Mar_Mask <- mask(x = Mar, mask = HALE)
Mar_CropET <- crop(x = Mar_Mask, y = extent(HALE))
if (RFUnit == " in") { Mar_CropET <- Mar_CropET  * 0.0393701 }
MarMET   <- round(cellStats(Mar_CropET, 'mean'), 0)
MarMETx   <- round(cellStats(Mar_CropET, 'max'), 0)
MarMETn   <- round(cellStats(Mar_CropET, 'min'), 0)

Apr <- raster(Mean_CLIM[15])
Apr_Mask <- mask(x = Apr, mask = HALE)
Apr_CropET <- crop(x = Apr_Mask, y = extent(HALE))
if (RFUnit == " in") { Apr_CropET <- Apr_CropET  * 0.0393701 }
AprMET   <- round(cellStats(Apr_CropET, 'mean'), 0)
AprMETx   <- round(cellStats(Apr_CropET, 'max'), 0)
AprMETn   <- round(cellStats(Apr_CropET, 'min'), 0)

May <- raster(Mean_CLIM[23])
May_Mask <- mask(x = May, mask = HALE)
May_CropET <- crop(x = May_Mask, y = extent(HALE))
if (RFUnit == " in") { May_CropET <- May_CropET  * 0.0393701 }
MayMET   <- round(cellStats(May_CropET, 'mean'), 0)
MayMETx   <- round(cellStats(May_CropET, 'max'), 0)
MayMETn   <- round(cellStats(May_CropET, 'min'), 0)

Jun <- raster(Mean_CLIM[21])
Jun_Mask <- mask(x = Jun, mask = HALE)
Jun_CropET <- crop(x = Jun_Mask, y = extent(HALE))
if (RFUnit == " in") { Jun_CropET <- Jun_CropET  * 0.0393701 }
JunMET   <- round(cellStats(Jun_CropET, 'mean'), 0)
JunMETx   <- round(cellStats(Jun_CropET, 'max'), 0)
JunMETn   <- round(cellStats(Jun_CropET, 'min'), 0)

Jul <- raster(Mean_CLIM[20])
Jul_Mask <- mask(x = Jul, mask = HALE)
Jul_CropET <- crop(x = Jul_Mask, y = extent(HALE))
if (RFUnit == " in") { Jul_CropET <- Jul_CropET  * 0.0393701 }
JulMET   <- round(cellStats(Jul_CropET, 'mean'), 0)
JulMETx   <- round(cellStats(Jul_CropET, 'max'), 0)
JulMETn   <- round(cellStats(Jul_CropET, 'min'), 0)

Aug <- raster(Mean_CLIM[16])
Aug_Mask <- mask(x = Aug, mask = HALE)
Aug_CropET <- crop(x = Aug_Mask, y = extent(HALE))
if (RFUnit == " in") { Aug_CropET <- Aug_CropET  * 0.0393701 }
AugMET   <- round(cellStats(Aug_CropET, 'mean'), 0)
AugMETx   <- round(cellStats(Aug_CropET, 'max'), 0)
AugMETn   <- round(cellStats(Aug_CropET, 'min'), 0)

Sep <- raster(Mean_CLIM[26])
Sep_Mask <- mask(x = Sep, mask = HALE)
Sep_CropET <- crop(x = Sep_Mask, y = extent(HALE))
if (RFUnit == " in") { Sep_CropET <- Sep_CropET  * 0.0393701 }
SepMET   <- round(cellStats(Sep_CropET, 'mean'), 0)
SepMETx   <- round(cellStats(Sep_CropET, 'max'), 0)
SepMETn   <- round(cellStats(Sep_CropET, 'min'), 0)

Oct <- raster(Mean_CLIM[25])
Oct_Mask <- mask(x = Oct, mask = HALE)
Oct_CropET <- crop(x = Oct_Mask, y = extent(HALE))
if (RFUnit == " in") { Oct_CropET <- Oct_CropET  * 0.0393701 }
OctMET   <- round(cellStats(Oct_CropET, 'mean'), 0)
OctMETx   <- round(cellStats(Oct_CropET, 'max'), 0)
OctMETn   <- round(cellStats(Oct_CropET, 'min'), 0)

Nov <- raster(Mean_CLIM[24])
Nov_Mask <- mask(x = Nov , mask = HALE)
Nov_CropET <- crop(x = Nov_Mask, y = extent(HALE))
if (RFUnit == " in") { Nov_CropET <- Nov_CropET  * 0.0393701 }
NovMET   <- round(cellStats(Nov_CropET, 'mean'), 0)
NovMETx   <- round(cellStats(Nov_CropET, 'max'), 0)
NovMETn   <- round(cellStats(Nov_CropET, 'min'), 0)

Dec <- raster(Mean_CLIM[17])
Dec_Mask <- mask(x = Dec, mask = HALE)
Dec_CropET <- crop(x = Dec_Mask, y = extent(HALE))
if (RFUnit == " in") { Dec_CropET <- Dec_CropET  * 0.0393701 }
DecMET   <- round(cellStats(Dec_CropET, 'mean'), 0)
DecMETx   <- round(cellStats(Dec_CropET, 'max'), 0)
DecMETn   <- round(cellStats(Dec_CropET, 'min'), 0)

Ann <- raster(Mean_CLIM[14])
Ann_Mask <- mask(x = Ann, mask = HALE)
Ann_CropET <- crop(x = Ann_Mask, y = extent(HALE))
if (RFUnit == " in") { Ann_CropET <- Ann_CropET  * 0.0393701 }
AnnMET   <- round(cellStats(Ann_CropET, 'mean'), 0)
AnnMETx   <- round(cellStats(Ann_CropET, 'max'), 0)
AnnMETn   <- round(cellStats(Ann_CropET, 'min'), 0)

MEANET2 <- c(JanMET,FebMET,MarMET,AprMET,MayMET,JunMET,JulMET,AugMET,SepMET,OctMET,NovMET,DecMET,AnnMET)


#Get thresholds for scale Bar
ETUP <- max(JanMETx,FebMETx,MarMETx,AprMETx,MayMETx,JunMETx,JulMETx,AugMETx,SepMETx,OctMETx,NovMETx,DecMETx,AnnMETx)
ETLo <- min(JanMETn,FebMETn,MarMETn,AprMETn,MayMETn,JunMETn,JulMETn,AugMETn,SepMETn,OctMETn,NovMETn,DecMETn,AnnMETn)

Cell.CL_Year[7, 2:14] <- MEANET2
Cell.CL_Year

##########   Cloud Frequency
debug_print("45, Cloud Frequency")

Jan <- raster(Mean_CLIM[6])
Jan_Mask <- mask(x = Jan, mask = HALE)
Jan_CropCF <- crop(x = Jan_Mask, y = extent(HALE))
JanMCF   <- round(cellStats(Jan_CropCF, 'mean'), 2)
JanMCFx   <- round(cellStats(Jan_CropCF, 'max'), 2)
JanMCFn   <- round(cellStats(Jan_CropCF, 'min'), 2)

Feb <- raster(Mean_CLIM[5])
Feb_Mask <- mask(x = Feb, mask = HALE)
Feb_CropCF <- crop(x = Feb_Mask, y = extent(HALE))
FebMCF   <- round(cellStats(Feb_CropCF, 'mean'), 2)
FebMCFx   <- round(cellStats(Feb_CropCF, 'max'), 2)
FebMCFn   <- round(cellStats(Feb_CropCF, 'min'), 2)

Mar <- raster(Mean_CLIM[9])
Mar_Mask <- mask(x = Mar, mask = HALE)
Mar_CropCF <- crop(x = Mar_Mask, y = extent(HALE))
MarMCF   <- round(cellStats(Mar_CropCF, 'mean'), 2)
MarMCFx   <- round(cellStats(Mar_CropCF, 'max'), 2)
MarMCFn   <- round(cellStats(Mar_CropCF, 'min'), 2)

Apr <- raster(Mean_CLIM[2])
Apr_Mask <- mask(x = Apr, mask = HALE)
Apr_CropCF <- crop(x = Apr_Mask, y = extent(HALE))
AprMCF   <- round(cellStats(Apr_CropCF, 'mean'), 2)
AprMCFx   <- round(cellStats(Apr_CropCF, 'max'), 2)
AprMCFn   <- round(cellStats(Apr_CropCF, 'min'), 2)

May <- raster(Mean_CLIM[10])
May_Mask <- mask(x = May, mask = HALE)
May_CropCF <- crop(x = May_Mask, y = extent(HALE))
MayMCF   <- round(cellStats(May_CropCF, 'mean'), 2)
MayMCFx   <- round(cellStats(May_CropCF, 'max'), 2)
MayMCFn   <- round(cellStats(May_CropCF, 'min'), 2)

Jun <- raster(Mean_CLIM[8])
Jun_Mask <- mask(x = Jun, mask = HALE)
Jun_CropCF <- crop(x = Jun_Mask, y = extent(HALE))
JunMCF   <- round(cellStats(Jun_CropCF, 'mean'), 2)
JunMCFx   <- round(cellStats(Jun_CropCF, 'max'), 2)
JunMCFn   <- round(cellStats(Jun_CropCF, 'min'), 2)

Jul <- raster(Mean_CLIM[7])
Jul_Mask <- mask(x = Jul, mask = HALE)
Jul_CropCF <- crop(x = Jul_Mask, y = extent(HALE))
JulMCF   <- round(cellStats(Jul_CropCF, 'mean'), 2)
JulMCFx   <- round(cellStats(Jul_CropCF, 'max'), 2)
JulMCFn   <- round(cellStats(Jul_CropCF, 'min'), 2)

Aug <- raster(Mean_CLIM[3])
Aug_Mask <- mask(x = Aug, mask = HALE)
Aug_CropCF <- crop(x = Aug_Mask, y = extent(HALE))
AugMCF   <- round(cellStats(Aug_CropCF, 'mean'), 2)
AugMCFx   <- round(cellStats(Aug_CropCF, 'max'), 2)
AugMCFn   <- round(cellStats(Aug_CropCF, 'min'), 2)

Sep <- raster(Mean_CLIM[13])
Sep_Mask <- mask(x = Sep, mask = HALE)
Sep_CropCF <- crop(x = Sep_Mask, y = extent(HALE))
SepMCF   <- round(cellStats(Sep_CropCF, 'mean'), 2)
SepMCFx   <- round(cellStats(Sep_CropCF, 'max'), 2)
SepMCFn   <- round(cellStats(Sep_CropCF, 'min'), 2)

Oct <- raster(Mean_CLIM[12])
Oct_Mask <- mask(x = Oct, mask = HALE)
Oct_CropCF <- crop(x = Oct_Mask, y = extent(HALE))
OctMCF   <- round(cellStats(Oct_CropCF, 'mean'), 2)
OctMCFx   <- round(cellStats(Oct_CropCF, 'max'), 2)
OctMCFn   <- round(cellStats(Oct_CropCF, 'min'), 2)

Nov <- raster(Mean_CLIM[11])
Nov_Mask <- mask(x = Nov , mask = HALE)
Nov_CropCF <- crop(x = Nov_Mask, y = extent(HALE))
NovMCF   <- round(cellStats(Nov_CropCF, 'mean'), 2)
NovMCFx   <- round(cellStats(Nov_CropCF, 'max'), 2)
NovMCFn   <- round(cellStats(Nov_CropCF, 'min'), 2)

Dec <- raster(Mean_CLIM[4])
Dec_Mask <- mask(x = Dec, mask = HALE)
Dec_CropCF <- crop(x = Dec_Mask, y = extent(HALE))
DecMCF   <- round(cellStats(Dec_CropCF, 'mean'), 2)
DecMCFx   <- round(cellStats(Dec_CropCF, 'max'), 2)
DecMCFn   <- round(cellStats(Dec_CropCF, 'min'), 2)

Ann <- raster(Mean_CLIM[1])
Ann_Mask <- mask(x = Ann, mask = HALE)
Ann_CropCF <- crop(x = Ann_Mask, y = extent(HALE))
AnnMCF   <- round(cellStats(Ann_CropCF, 'mean'), 2)
AnnMCFx   <- round(cellStats(Ann_CropCF, 'max'), 2)
AnnMCFn   <- round(cellStats(Ann_CropCF, 'min'), 2)

MEANCF2 <- c(JanMCF,FebMCF,MarMCF,AprMCF,MayMCF,JunMCF,JulMCF,AugMCF,SepMCF,OctMCF,NovMCF,DecMCF,AnnMCF)
MEANCF3 <- MEANCF2 * 100

#Get thresholds for scale Bar
CFUP <- max(JanMCFx,FebMCFx,MarMCFx,AprMCFx,MayMCFx,JunMCFx,JulMCFx,AugMCFx,SepMCFx,OctMCFx,NovMCFx,DecMCFx,AnnMCFx)
CFLO <- min(JanMCFn,FebMCFn,MarMCFn,AprMCFn,MayMCFn,JunMCFn,JulMCFn,AugMCFn,SepMCFn,OctMCFn,NovMCFn,DecMCFn,AnnMCFn)

Cell.CL_Year[6, 2:14] <- MEANCF3
Cell.CL_Year

##########

##########   MEAN TEMPERATURE
debug_print("46, MEAN TEMPERATURE")

Jan <- raster(Mean_CLIM[71])
Jan_Mask <- mask(x = Jan, mask = HALE)
Jan_CropTA <- crop(x = Jan_Mask, y = extent(HALE))
if (TUnit == "°F") { Jan_CropTA <- (Jan_CropTA * 1.8) + 32 }
JanMTA   <- round(cellStats(Jan_CropTA, 'mean'), 1)
JanMTAx   <- round(cellStats(Jan_CropTA, 'max'), 1)
JanMTAn   <- round(cellStats(Jan_CropTA, 'min'), 1)

Feb <- raster(Mean_CLIM[70])
Feb_Mask <- mask(x = Feb, mask = HALE)
Feb_CropTA <- crop(x = Feb_Mask, y = extent(HALE))
if (TUnit == "°F") { Feb_CropTA <- (Feb_CropTA * 1.8) + 32 }
FebMTA   <- round(cellStats(Feb_CropTA, 'mean'), 1)
FebMTAx   <- round(cellStats(Feb_CropTA, 'max'), 1)
FebMTAn   <- round(cellStats(Feb_CropTA, 'min'), 1)

Mar <- raster(Mean_CLIM[74])
Mar_Mask <- mask(x = Mar, mask = HALE)
Mar_CropTA <- crop(x = Mar_Mask, y = extent(HALE))
if (TUnit == "°F") { Mar_CropTA <- (Mar_CropTA * 1.8) + 32 }
MarMTA   <- round(cellStats(Mar_CropTA, 'mean'), 1)
MarMTAx   <- round(cellStats(Mar_CropTA, 'max'), 1)
MarMTAn   <- round(cellStats(Mar_CropTA, 'min'), 1)

Apr <- raster(Mean_CLIM[67])
Apr_Mask <- mask(x = Apr, mask = HALE)
Apr_CropTA <- crop(x = Apr_Mask, y = extent(HALE))
if (TUnit == "°F") { Apr_CropTA <- (Apr_CropTA * 1.8) + 32 }
AprMTA   <- round(cellStats(Apr_CropTA, 'mean'), 1)
AprMTAx   <- round(cellStats(Apr_CropTA, 'max'), 1)
AprMTAn   <- round(cellStats(Apr_CropTA, 'min'), 1)

May <- raster(Mean_CLIM[75])
May_Mask <- mask(x = May, mask = HALE)
May_CropTA <- crop(x = May_Mask, y = extent(HALE))
if (TUnit == "°F") { May_CropTA <- (May_CropTA * 1.8) + 32 }
MayMTA   <- round(cellStats(May_CropTA, 'mean'), 1)
MayMTAx   <- round(cellStats(May_CropTA, 'max'), 1)
MayMTAn   <- round(cellStats(May_CropTA, 'min'), 1)

Jun <- raster(Mean_CLIM[73])
Jun_Mask <- mask(x = Jun, mask = HALE)
Jun_CropTA <- crop(x = Jun_Mask, y = extent(HALE))
if (TUnit == "°F") { Jun_CropTA <- (Jun_CropTA * 1.8) + 32 }
JunMTA   <- round(cellStats(Jun_CropTA, 'mean'), 1)
JunMTAx   <- round(cellStats(Jun_CropTA, 'max'), 1)
JunMTAn   <- round(cellStats(Jun_CropTA, 'min'), 1)

Jul <- raster(Mean_CLIM[72])
Jul_Mask <- mask(x = Jul, mask = HALE)
Jul_CropTA <- crop(x = Jul_Mask, y = extent(HALE))
if (TUnit == "°F") { Jul_CropTA <- (Jul_CropTA * 1.8) + 32 }
JulMTA   <- round(cellStats(Jul_CropTA, 'mean'), 1)
JulMTAx   <- round(cellStats(Jul_CropTA, 'max'), 1)
JulMTAn   <- round(cellStats(Jul_CropTA, 'min'), 1)

Aug <- raster(Mean_CLIM[68])
Aug_Mask <- mask(x = Aug, mask = HALE)
Aug_CropTA <- crop(x = Aug_Mask, y = extent(HALE))
if (TUnit == "°F") { Aug_CropTA <- (Aug_CropTA * 1.8) + 32 }
AugMTA   <- round(cellStats(Aug_CropTA, 'mean'), 1)
AugMTAx   <- round(cellStats(Aug_CropTA, 'max'), 1)
AugMTAn   <- round(cellStats(Aug_CropTA, 'min'), 1)

Sep <- raster(Mean_CLIM[78])
Sep_Mask <- mask(x = Sep, mask = HALE)
Sep_CropTA <- crop(x = Sep_Mask, y = extent(HALE))
if (TUnit == "°F") { Sep_CropTA <- (Sep_CropTA * 1.8) + 32 }
SepMTA   <- round(cellStats(Sep_CropTA, 'mean'), 1)
SepMTAx   <- round(cellStats(Sep_CropTA, 'max'), 1)
SepMTAn   <- round(cellStats(Sep_CropTA, 'min'), 1)

Oct <- raster(Mean_CLIM[77])
Oct_Mask <- mask(x = Oct, mask = HALE)
Oct_CropTA <- crop(x = Oct_Mask, y = extent(HALE))
if (TUnit == "°F") { Oct_CropTA <- (Oct_CropTA * 1.8) + 32 }
OctMTA   <- round(cellStats(Oct_CropTA, 'mean'), 1)
OctMTAx   <- round(cellStats(Oct_CropTA, 'max'), 1)
OctMTAn   <- round(cellStats(Oct_CropTA, 'min'), 1)

Nov <- raster(Mean_CLIM[76])
Nov_Mask <- mask(x = Nov , mask = HALE)
Nov_CropTA <- crop(x = Nov_Mask, y = extent(HALE))
if (TUnit == "°F") { Nov_CropTA <- (Nov_CropTA * 1.8) + 32 }
NovMTA   <- round(cellStats(Nov_CropTA, 'mean'), 1)
NovMTAx   <- round(cellStats(Nov_CropTA, 'max'), 1)
NovMTAn   <- round(cellStats(Nov_CropTA, 'min'), 1)

Dec <- raster(Mean_CLIM[69])
Dec_Mask <- mask(x = Dec, mask = HALE)
Dec_CropTA <- crop(x = Dec_Mask, y = extent(HALE))
if (TUnit == "°F") { Dec_CropTA <- (Dec_CropTA * 1.8) + 32 }
DecMTA   <- round(cellStats(Dec_CropTA, 'mean'), 1)
DecMTAx   <- round(cellStats(Dec_CropTA, 'max'), 1)
DecMTAn   <- round(cellStats(Dec_CropTA, 'min'), 1)

Ann <- raster(Mean_CLIM[66])
Ann_Mask <- mask(x = Ann, mask = HALE)
Ann_CropTA <- crop(x = Ann_Mask, y = extent(HALE))
if (TUnit == "°F") { Ann_CropTA <- (Ann_CropTA * 1.8) + 32 }
AnnMTA   <- round(cellStats(Ann_CropTA, 'mean'), 1)
AnnMTAx   <- round(cellStats(Ann_CropTA, 'max'), 1)
AnnMTAn   <- round(cellStats(Ann_CropTA, 'min'), 1)

MEANTA2 <- c(JanMTA,FebMTA,MarMTA,AprMTA,MayMTA,JunMTA,JulMTA,AugMTA,SepMTA,OctMTA,NovMTA,DecMTA,AnnMTA)

#Get thresholds for scale Bar
TaUP <- max(JanMTAx,FebMTAx,MarMTAx,AprMTAx,MayMTAx,JunMTAx,JulMTAx,AugMTAx,SepMTAx,OctMTAx,NovMTAx,DecMTAx,AnnMTAx)
TaLo <- min(JanMTAn,FebMTAn,MarMTAn,AprMTAn,MayMTAn,JunMTAn,JulMTAn,AugMTAn,SepMTAn,OctMTAn,NovMTAn,DecMTAn,AnnMTAn)

Cell.CL_Year[3, 2:14] <- MEANTA2
Cell.CL_Year

########## MAX TEMPERATURE
debug_print("47, MAX TEMPERATURE")

Jan <- raster(Mean_CLIM[84])
crs(Jan) <- crs(EXAMP)
Jan_Mask <- mask(x = Jan, mask = HALE)
Jan_CropTX <- crop(x = Jan_Mask, y = extent(HALE))
if (TUnit == "°F") { Jan_CropTX <- (Jan_CropTX * 1.8) + 32 }
JanMTX  <- round(cellStats(Jan_CropTX , 'mean'), 1)
JanMTXx   <- round(cellStats(Jan_CropTX, 'max'), 1)
JanMTXn   <- round(cellStats(Jan_CropTX, 'min'), 1)

Feb <- raster(Mean_CLIM[83])
crs(Feb) <- crs(EXAMP)
Feb_Mask <- mask(x = Feb, mask = HALE)
Feb_CropTX  <- crop(x = Feb_Mask, y = extent(HALE))
if (TUnit == "°F") { Feb_CropTX <- (Feb_CropTX * 1.8) + 32 }
FebMTX   <- round(cellStats(Feb_CropTX , 'mean'), 1)
FebMTXx   <- round(cellStats(Feb_CropTX, 'max'), 1)
FebMTXn   <- round(cellStats(Feb_CropTX, 'min'), 1)

Mar <- raster(Mean_CLIM[87])
crs(Mar) <- crs(EXAMP)
Mar_Mask <- mask(x = Mar, mask = HALE)
Mar_CropTX  <- crop(x = Mar_Mask, y = extent(HALE))
if (TUnit == "°F") { Mar_CropTX <- (Mar_CropTX * 1.8) + 32 }
MarMTX   <- round(cellStats(Mar_CropTX , 'mean'), 1)
MarMTXx   <- round(cellStats(Mar_CropTX, 'max'), 1)
MarMTXn   <- round(cellStats(Mar_CropTX, 'min'), 1)

Apr <- raster(Mean_CLIM[80])
crs(Apr) <- crs(EXAMP)
Apr_Mask <- mask(x = Apr, mask = HALE)
Apr_CropTX  <- crop(x = Apr_Mask, y = extent(HALE))
if (TUnit == "°F") { Apr_CropTX <- (Apr_CropTX * 1.8) + 32 }
AprMTX   <- round(cellStats(Apr_CropTX , 'mean'), 1)
AprMTXx   <- round(cellStats(Apr_CropTX, 'max'), 1)
AprMTXn   <- round(cellStats(Apr_CropTX, 'min'), 1)

May <- raster(Mean_CLIM[88])
crs(May) <- crs(EXAMP)
May_Mask <- mask(x = May, mask = HALE)
May_CropTX  <- crop(x = May_Mask, y = extent(HALE))
if (TUnit == "°F") { May_CropTX <- (May_CropTX * 1.8) + 32 }
MayMTX    <- round(cellStats(May_CropTX, 'mean'), 1)
MayMTXx   <- round(cellStats(May_CropTX, 'max'), 1)
MayMTXn   <- round(cellStats(May_CropTX, 'min'), 1)

Jun <- raster(Mean_CLIM[86])
crs(Jun) <- crs(EXAMP)
Jun_Mask <- mask(x = Jun, mask = HALE)
Jun_CropTX  <- crop(x = Jun_Mask, y = extent(HALE))
if (TUnit == "°F") { Jun_CropTX <- (Jun_CropTX * 1.8) + 32 }
JunMTX   <- round(cellStats(Jun_CropTX , 'mean'), 1)
JunMTXx   <- round(cellStats(Jun_CropTX, 'max'), 1)
JunMTXn   <- round(cellStats(Jun_CropTX, 'min'), 1)

Jul <- raster(Mean_CLIM[85])
crs(Jul) <- crs(EXAMP)
Jul_Mask <- mask(x = Jul, mask = HALE)
Jul_CropTX  <- crop(x = Jul_Mask, y = extent(HALE))
if (TUnit == "°F") { Jul_CropTX <- (Jul_CropTX * 1.8) + 32 }
JulMTX    <- round(cellStats(Jul_CropTX, 'mean'), 1)
JulMTXx   <- round(cellStats(Jul_CropTX, 'max'), 1)
JulMTXn   <- round(cellStats(Jul_CropTX, 'min'), 1)

Aug <- raster(Mean_CLIM[81])
crs(Aug) <- crs(EXAMP)
Aug_Mask <- mask(x = Aug, mask = HALE)
Aug_CropTX  <- crop(x = Aug_Mask, y = extent(HALE))
if (TUnit == "°F") { Aug_CropTX <- (Aug_CropTX * 1.8) + 32 }
AugMTX   <- round(cellStats(Aug_CropTX, 'mean'), 1)
AugMTXx   <- round(cellStats(Aug_CropTX, 'max'), 1)
AugMTXn   <- round(cellStats(Aug_CropTX, 'min'), 1)

Sep <- raster(Mean_CLIM[91])
crs(Sep) <- crs(EXAMP)
Sep_Mask <- mask(x = Sep, mask = HALE)
Sep_CropTX  <- crop(x = Sep_Mask, y = extent(HALE))
if (TUnit == "°F") { Sep_CropTX <- (Sep_CropTX * 1.8) + 32 }
SepMTX   <- round(cellStats(Sep_CropTX, 'mean'), 1)
SepMTXx   <- round(cellStats(Sep_CropTX, 'max'), 1)
SepMTXn   <- round(cellStats(Sep_CropTX, 'min'), 1)

Oct <- raster(Mean_CLIM[90])
crs(Oct) <- crs(EXAMP)
Oct_Mask <- mask(x = Oct, mask = HALE)
Oct_CropTX  <- crop(x = Oct_Mask, y = extent(HALE))
if (TUnit == "°F") { Oct_CropTX <- (Oct_CropTX * 1.8) + 32 }
OctMTX   <- round(cellStats(Oct_CropTX, 'mean'), 1)
OctMTXx   <- round(cellStats(Oct_CropTX, 'max'), 1)
OctMTXn   <- round(cellStats(Oct_CropTX, 'min'), 1)

Nov <- raster(Mean_CLIM[89])
crs(Nov) <- crs(EXAMP)
Nov_Mask <- mask(x = Nov , mask = HALE)
Nov_CropTX  <- crop(x = Nov_Mask, y = extent(HALE))
if (TUnit == "°F") { Nov_CropTX <- (Nov_CropTX * 1.8) + 32 }
NovMTX   <- round(cellStats(Nov_CropTX, 'mean'), 1)
NovMTXx   <- round(cellStats(Nov_CropTX, 'max'), 1)
NovMTXn   <- round(cellStats(Nov_CropTX, 'min'), 1)

Dec <- raster(Mean_CLIM[82])
crs(Dec) <- crs(EXAMP)
Dec_Mask <- mask(x = Dec, mask = HALE)
Dec_CropTX  <- crop(x = Dec_Mask, y = extent(HALE))
if (TUnit == "°F") { Dec_CropTX <- (Dec_CropTX * 1.8) + 32 }
DecMTX   <- round(cellStats(Dec_CropTX, 'mean'), 1)
DecMTXx   <- round(cellStats(Dec_CropTX, 'max'), 1)
DecMTXn   <- round(cellStats(Dec_CropTX, 'min'), 1)

Ann <- raster(Mean_CLIM[79])
crs(Ann) <- crs(EXAMP)
Ann_Mask <- mask(x = Ann, mask = HALE)
Ann_CropTX <- crop(x = Ann_Mask, y = extent(HALE))
if (TUnit == "°F") { Ann_CropTX <- (Ann_CropTX * 1.8) + 32 }
AnnMTX   <- round(cellStats(Ann_CropTX, 'mean'), 1)
AnnMTXx   <- round(cellStats(Ann_CropTX, 'max'), 1)
AnnMTXn   <- round(cellStats(Ann_CropTX, 'min'), 1)

MEANTX2 <- c(JanMTX,FebMTX,MarMTX,AprMTX,MayMTX,JunMTX,JulMTX,AugMTX,SepMTX,OctMTX,NovMTX,DecMTX,AnnMTX)

##########   Add to dataframe

Cell.CL_Year[4,2:14] <- MEANTX2

##########   Get thresholds for figures 

TxUP <- max(JanMTXx,FebMTXx,MarMTXx,AprMTXx,MayMTXx,JunMTXx,JulMTXx,AugMTXx,SepMTXx,OctMTXx,NovMTXx,DecMTXx,AnnMTXx)
TxLO <- min(JanMTXn,FebMTXn,MarMTXn,AprMTXn,MayMTXn,JunMTXn,JulMTXn,AugMTXn,SepMTXn,OctMTXn,NovMTXn,DecMTXn,AnnMTXn)

##########    MIN TEMPERATURE
debug_print("48, MIN TEMPERATURE")

Jan <- raster(Mean_CLIM[97])
crs(Jan) <- crs(EXAMP)
Jan_Mask <- mask(x = Jan, mask = HALE)
Jan_CropTN <- crop(x = Jan_Mask, y = extent(HALE))
if (TUnit == "°F") { Jan_CropTN <- (Jan_CropTN * 1.8) + 32 }
JanMTN   <- round(cellStats(Jan_CropTN, 'mean'), 1)
JanMTNx   <- round(cellStats(Jan_CropTN , 'max'), 1)
JanMTNn   <- round(cellStats(Jan_CropTN , 'min'), 1)
summary(Jan_CropTN)

plot(Jan_CropTN)
Feb <- raster(Mean_CLIM[96])
crs(Feb) <- crs(EXAMP)
Feb_Mask <- mask(x = Feb, mask = HALE)
Feb_CropTN <- crop(x = Feb_Mask, y = extent(HALE))
if (TUnit == "°F") { Feb_CropTN <- (Feb_CropTN * 1.8) + 32 }
FebMTN  <- round(cellStats(Feb_CropTN, 'mean'), 1)
FebMTNx   <- round(cellStats(Feb_CropTN , 'max'), 1)
FebMTNn   <- round(cellStats(Feb_CropTN , 'min'), 1)

Mar <- raster(Mean_CLIM[100])
crs(Mar) <- crs(EXAMP)
Mar_Mask <- mask(x = Mar, mask = HALE)
Mar_CropTN <- crop(x = Mar_Mask, y = extent(HALE))
if (TUnit == "°F") { Mar_CropTN <- (Mar_CropTN * 1.8) + 32 }
MarMTN   <- round(cellStats(Mar_CropTN, 'mean'), 1)
MarMTNx   <- round(cellStats(Mar_CropTN , 'max'), 1)
MarMTNn   <- round(cellStats(Mar_CropTN , 'min'), 1)

Apr <- raster(Mean_CLIM[93])
crs(Apr) <- crs(EXAMP)
Apr_Mask <- mask(x = Apr, mask = HALE)
Apr_CropTN <- crop(x = Apr_Mask, y = extent(HALE))
if (TUnit == "°F") { Apr_CropTN <- (Apr_CropTN * 1.8) + 32 }
AprMTN   <- round(cellStats(Apr_CropTN, 'mean'), 1)
AprMTNx   <- round(cellStats(Apr_CropTN , 'max'), 1)
AprMTNn   <- round(cellStats(Apr_CropTN , 'min'), 1)

May <- raster(Mean_CLIM[101])
crs(May) <- crs(EXAMP)
May_Mask <- mask(x = May, mask = HALE)
May_CropTN <- crop(x = May_Mask, y = extent(HALE))
if (TUnit == "°F") { May_CropTN <- (May_CropTN * 1.8) + 32 }
MayMTN   <- round(cellStats(May_CropTN, 'mean'), 1)
MayMTNx   <- round(cellStats(May_CropTN , 'max'), 1)
MayMTNn   <- round(cellStats(May_CropTN , 'min'), 1)

Jun <- raster(Mean_CLIM[99])
crs(Jun) <- crs(EXAMP)
Jun_Mask <- mask(x = Jun, mask = HALE)
Jun_CropTN <- crop(x = Jun_Mask, y = extent(HALE))
if (TUnit == "°F") { Jun_CropTN <- (Jun_CropTN * 1.8) + 32 }
JunMTN   <- round(cellStats(Jun_CropTN, 'mean'), 1)
JunMTNx   <- round(cellStats(Jun_CropTN , 'max'), 1)
JunMTNn   <- round(cellStats(Jun_CropTN , 'min'), 1)

Jul <- raster(Mean_CLIM[98])
crs(Jul) <- crs(EXAMP)
Jul_Mask <- mask(x = Jul, mask = HALE)
Jul_CropTN <- crop(x = Jul_Mask, y = extent(HALE))
if (TUnit == "°F") { Jul_CropTN <- (Jul_CropTN * 1.8) + 32 }
JulMTN   <- round(cellStats(Jul_CropTN, 'mean'), 1)
JulMTNx   <- round(cellStats(Jul_CropTN , 'max'), 1)
JulMTNn   <- round(cellStats(Jul_CropTN , 'min'), 1)

Aug <- raster(Mean_CLIM[94])
crs(Aug) <- crs(EXAMP)
Aug_Mask <- mask(x = Aug, mask = HALE)
Aug_CropTN <- crop(x = Aug_Mask, y = extent(HALE))
if (TUnit == "°F") { Aug_CropTN <- (Aug_CropTN * 1.8) + 32 }
AugMTN   <- round(cellStats(Aug_CropTN, 'mean'), 1)
AugMTNx   <- round(cellStats(Aug_CropTN , 'max'), 1)
AugMTNn   <- round(cellStats(Aug_CropTN , 'min'), 1)

Sep <- raster(Mean_CLIM[104])
crs(Sep) <- crs(EXAMP)
Sep_Mask <- mask(x = Sep, mask = HALE)
Sep_CropTN <- crop(x = Sep_Mask, y = extent(HALE))
if (TUnit == "°F") { Sep_CropTN <- (Sep_CropTN * 1.8) + 32 }
SepMTN   <- round(cellStats(Sep_CropTN, 'mean'), 1)
SepMTNx   <- round(cellStats(Sep_CropTN , 'max'), 1)
SepMTNn   <- round(cellStats(Sep_CropTN , 'min'), 1)

Oct <- raster(Mean_CLIM[103])
crs(Oct) <- crs(EXAMP)
Oct_Mask <- mask(x = Oct, mask = HALE)
Oct_CropTN <- crop(x = Oct_Mask, y = extent(HALE))
if (TUnit == "°F") { Oct_CropTN <- (Oct_CropTN * 1.8) + 32 }
OctMTN   <- round(cellStats(Oct_CropTN, 'mean'), 1)
OctMTNx   <- round(cellStats(Oct_CropTN , 'max'), 1)
OctMTNn   <- round(cellStats(Oct_CropTN , 'min'), 1)

Nov <- raster(Mean_CLIM[102])
crs(Nov) <- crs(EXAMP)
Nov_Mask <- mask(x = Nov , mask = HALE)
Nov_CropTN <- crop(x = Nov_Mask, y = extent(HALE))
if (TUnit == "°F") { Nov_CropTN <- (Nov_CropTN * 1.8) + 32 }
NovMTN   <- round(cellStats(Nov_CropTN, 'mean'), 1)
NovMTNx   <- round(cellStats(Nov_CropTN , 'max'), 1)
NovMTNn   <- round(cellStats(Nov_CropTN , 'min'), 1)

Dec <- raster(Mean_CLIM[95])
crs(Dec) <- crs(EXAMP)
Dec_Mask <- mask(x = Dec, mask = HALE)
Dec_CropTN <- crop(x = Dec_Mask, y = extent(HALE))
if (TUnit == "°F") { Dec_CropTN <- (Dec_CropTN * 1.8) + 32 }
DecMTN   <- round(cellStats(Dec_CropTN, 'mean'), 1)
DecMTNx   <- round(cellStats(Dec_CropTN , 'max'), 1)
DecMTNn   <- round(cellStats(Dec_CropTN , 'min'), 1)

Ann <- raster(Mean_CLIM[92])
crs(Ann) <- crs(EXAMP)
Ann_Mask <- mask(x = Ann, mask = HALE)
Ann_CropTN <- crop(x = Ann_Mask, y = extent(HALE))
if (TUnit == "°F") { Ann_CropTN <- (Ann_CropTN * 1.8) + 32 }
AnnMTN   <- round(cellStats(Ann_CropTN, 'mean'), 1)
AnnMTNx   <- round(cellStats(Ann_CropTN , 'max'), 1)
AnnMTNn   <- round(cellStats(Ann_CropTN , 'min'), 1)

########## Add to Table

MEANTN2 <- c(JanMTN,FebMTN,MarMTN,AprMTN,MayMTN,JunMTN,JulMTN,AugMTN,SepMTN,OctMTN,NovMTN,DecMTN,AnnMTN)
Cell.CL_Year[2,2:14] <- MEANTN2
Cell.CL_Year

########## Get thresholds for figures 

TnUP <- max(JanMTNx,FebMTNx,MarMTNx,AprMTNx,MayMTNx,JunMTNx,JulMTNx,AugMTNx,SepMTNx,OctMTNx,NovMTNx,DecMTNx,AnnMTNx)
TnLO <- min(JanMTNn,FebMTNn,MarMTNn,AprMTNn,MayMTNn,JunMTNn,JulMTNn,AugMTNn,SepMTNn,OctMTNn,NovMTNn,DecMTNn,AnnMTNn)

##########   RELATIVE HUMIDITY
debug_print("49, RELATIVE HUMIDITY")

Jan <- raster(Mean_CLIM[32])
Jan_Mask <- mask(x = Jan, mask = HALE)
Jan_CropRH <- crop(x = Jan_Mask, y = extent(HALE))
JanMRH    <- round(cellStats(Jan_CropRH , 'mean'), 0)
JanMRHx   <- round(cellStats(Jan_CropRH  , 'max'), 0)
JanMRHn   <- round(cellStats(Jan_CropRH  , 'min'), 0)

Feb <- raster(Mean_CLIM[31])
Feb_Mask <- mask(x = Feb, mask = HALE)
Feb_CropRH  <- crop(x = Feb_Mask, y = extent(HALE))
FebMRH    <- round(cellStats(Feb_CropRH , 'mean'), 0)
FebMRHx   <- round(cellStats(Feb_CropRH  , 'max'), 0)
FebMRHn   <- round(cellStats(Feb_CropRH  , 'min'), 0)

Mar <- raster(Mean_CLIM[35])
Mar_Mask <- mask(x = Mar, mask = HALE)
Mar_CropRH  <- crop(x = Mar_Mask, y = extent(HALE))
MarMRH    <- round(cellStats(Mar_CropRH , 'mean'), 0)
MarMRHx   <- round(cellStats(Mar_CropRH  , 'max'), 0)
MarMRHn   <- round(cellStats(Mar_CropRH  , 'min'), 0)

Apr <- raster(Mean_CLIM[28])
Apr_Mask <- mask(x = Apr, mask = HALE)
Apr_CropRH  <- crop(x = Apr_Mask, y = extent(HALE))
AprMRH    <- round(cellStats(Apr_CropRH , 'mean'), 0)
AprMRHx   <- round(cellStats(Apr_CropRH  , 'max'), 0)
AprMRHn   <- round(cellStats(Apr_CropRH  , 'min'), 0)

May <- raster(Mean_CLIM[36])
May_Mask <- mask(x = May, mask = HALE)
May_CropRH  <- crop(x = May_Mask, y = extent(HALE))
MayMRH   <- round(cellStats(May_CropRH, 'mean'), 0)
MayMRHx   <- round(cellStats(May_CropRH  , 'max'), 0)
MayMRHn   <- round(cellStats(May_CropRH  , 'min'), 0)

Jun <- raster(Mean_CLIM[34])
Jun_Mask <- mask(x = Jun, mask = HALE)
Jun_CropRH  <- crop(x = Jun_Mask, y = extent(HALE))
JunMRH    <- round(cellStats(Jun_CropRH, 'mean'), 0)
JunMRHx   <- round(cellStats(Jun_CropRH  , 'max'), 0)
JunMRHn   <- round(cellStats(Jun_CropRH  , 'min'), 0)

Jul <- raster(Mean_CLIM[33])
Jul_Mask <- mask(x = Jul, mask = HALE)
Jul_CropRH  <- crop(x = Jul_Mask, y = extent(HALE))
JulMRH   <- round(cellStats(Jul_CropRH, 'mean'), 0)
JulMRHx   <- round(cellStats(Jul_CropRH  , 'max'), 0)
JulMRHn   <- round(cellStats(Jul_CropRH  , 'min'), 0)

Aug <- raster(Mean_CLIM[29])
Aug_Mask <- mask(x = Aug, mask = HALE)
Aug_CropRH  <- crop(x = Aug_Mask, y = extent(HALE))
AugMRH   <- round(cellStats(Aug_CropRH, 'mean'), 0)
AugMRHx   <- round(cellStats(Aug_CropRH  , 'max'), 0)
AugMRHn   <- round(cellStats(Aug_CropRH  , 'min'), 0)

Sep <- raster(Mean_CLIM[39])
Sep_Mask <- mask(x = Sep, mask = HALE)
Sep_CropRH  <- crop(x = Sep_Mask, y = extent(HALE))
SepMRH <- round(cellStats(Sep_CropRH, 'mean'), 0)
SepMRHx   <- round(cellStats(Sep_CropRH  , 'max'), 0)
SepMRHn   <- round(cellStats(Sep_CropRH  , 'min'), 0)

Oct <- raster(Mean_CLIM[38])
Oct_Mask <- mask(x = Oct, mask = HALE)
Oct_CropRH  <- crop(x = Oct_Mask, y = extent(HALE))
OctMRH <- round(cellStats(Oct_CropRH, 'mean'), 0)
OctMRHx   <- round(cellStats(Oct_CropRH  , 'max'), 0)
OctMRHn   <- round(cellStats(Oct_CropRH  , 'min'), 0)

Nov <- raster(Mean_CLIM[37])
Nov_Mask <- mask(x = Nov , mask = HALE)
Nov_CropRH  <- crop(x = Nov_Mask, y = extent(HALE))
NovMRH <- round(cellStats(Nov_CropRH , 'mean'), 0)
NovMRHx   <- round(cellStats(Nov_CropRH  , 'max'), 0)
NovMRHn   <- round(cellStats(Nov_CropRH  , 'min'), 0)

Dec <- raster(Mean_CLIM[30])
Dec_Mask <- mask(x = Dec, mask = HALE)
Dec_CropRH  <- crop(x = Dec_Mask, y = extent(HALE))
DecMRH    <- round(cellStats(Dec_CropRH, 'mean'), 0)
DecMRHx   <- round(cellStats(Dec_CropRH  , 'max'), 0)
DecMRHn   <- round(cellStats(Dec_CropRH  , 'min'), 0)

ANN <- raster(Mean_CLIM[27])
ANN_Mask <- mask(x = ANN, mask = HALE)
ANN_CropRH  <- crop(x = ANN_Mask, y = extent(HALE))
ANNMRH    <- round(cellStats(ANN_CropRH, 'mean'), 0)
ANNMRHx   <- round(cellStats(ANN_CropRH  , 'max'), 0)
ANNMRHn   <- round(cellStats(ANN_CropRH  , 'min'), 0)

##########   Add to Table

MEANRH2 <- c(JanMRH,FebMRH ,MarMRH ,AprMRH ,MayMRH ,JunMRH ,JulMRH ,AugMRH ,SepMRH ,OctMRH ,NovMRH ,DecMRH, ANNMRH  )
Cell.CL_Year[5,2:14] <- MEANRH2 
Cell.CL_Year

##########   Get thresholds for figures 

RHUP <- max(JanMRHx,FebMRHx ,MarMRHx ,AprMRHx ,MayMRHx ,JunMRHx ,JulMRHx ,AugMRHx ,SepMRHx ,OctMRHx ,NovMRHx ,DecMRHx, ANNMRHx  )
RHLO <- min(JanMRHn,FebMRHn ,MarMRHn ,AprMRHn ,MayMRHn ,JunMRHn ,JulMRHn ,AugMRHn ,SepMRHn ,OctMRHn ,NovMRHn ,DecMRHn, ANNMRHn  )




#########   OTHER CLIMATE VARIABLES
debug_print("50, OTHER CLIMATE VARIABLES")

### Read in the Annual rasters for each variable only and extract min and max to show spatial range.
# Have extracted the min max from Annual layers already This next step can be streamlined

#See Varible "ANN_CropRH" for example
RH_P <- raster(Mean_CLIM[27])
RH_P_Mask <- mask(x = RH_P, mask = HALE)
RH_P_Crop <- crop(x = RH_P_Mask, y = extent(HALE))
RH_P_M   <- round(cellStats(RH_P_Crop, 'mean'), 0)
RHUPA   <- round(cellStats(RH_P_Crop  , 'max'), 1)
RHLOA   <- round(cellStats(RH_P_Crop  , 'min'), 1)

Tair_P <- raster(Mean_CLIM[66])
names(Tair_P) = "Mean Air Temp."
Tair_P_Mask <- mask(x = Tair_P, mask = HALE)
Tair_P_Crop <- crop(x = Tair_P_Mask, y = extent(HALE))
if (TUnit == "°F") { Tair_P_Crop <- (Tair_P_Crop * 1.8) + 32 }
Tair_P_M   <- round(cellStats(Tair_P_Crop, 'mean'), 1)
Tair_PName <- names(Tair_P)

SM_P <- raster(Mean_CLIM[40])
names(SM_P) = "Soil Moisture"
SM_P_Mask <- mask(x = SM_P, mask = HALE)
SM_P_Crop <- crop(x = SM_P_Mask, y = extent(HALE))
SM_P_M   <- round(cellStats(SM_P_Crop, 'mean'), 2)
SM_PName <- names(SM_P)
SMUP   <- round(cellStats(SM_P_Crop  , 'max'), 2)
SMLO   <- round(cellStats(SM_P_Crop  , 'min'), 2)

CF_P <- raster(Mean_CLIM[1])
names(CF_P) = "Cloud Freq"
CF_P_Mask <- mask(x = CF_P, mask = HALE)
CF_P_Crop <- crop(x = CF_P_Mask, y = extent(HALE))
CF_P_M   <- round(cellStats(CF_P_Crop, 'mean'), 2)
CF_PName <- names(CF_P)
CFUP   <- round(cellStats(CF_P_Crop  , 'max'), 2)
CFLO   <- round(cellStats(CF_P_Crop  , 'min'), 2)

Tmax_P <- raster(Mean_CLIM[79])
crs(Tmax_P)  <- crs(EXAMP)
names(Tmax_P) = "Max Air Temp."
Tmax_P_Mask <- mask(x = Tmax_P, mask = HALE)
Tmax_P_Crop <- crop(x = Tmax_P_Mask, y = extent(HALE))
if (TUnit == "°F") { Tmax_P_Crop <- (Tmax_P_Crop  * 1.8) + 32 }
Tmax_P_M   <- round(cellStats(Tmax_P_Crop, 'mean'), 1)
Tmax_PName <- names(Tmax_P)

Tmin_P <- raster(Mean_CLIM[92])
crs(Tmin_P)  <- crs(EXAMP)
names(Tmin_P) = "Mean Air Temp."
Tmin_P_Mask <- mask(x = Tmin_P, mask = HALE)
Tmin_P_Crop <- crop(x = Tmin_P_Mask, y = extent(HALE))
if (TUnit == "°F") { Tmin_P_Crop <- (Tmin_P_Crop  * 1.8) + 32 }
Tmin_P_M   <- round(cellStats(Tmin_P_Crop, 'mean'), 1)
Tmin_PName <- names(Tmin_P)

KD_P <- raster(Mean_CLIM[53])
names(KD_P) = "Solar Radiation."
KD_P_Mask <- mask(x = KD_P, mask = HALE)
KD_P_Crop <- crop(x = KD_P_Mask, y = extent(HALE))
KD_P_M   <- round(cellStats(KD_P_Crop, 'mean'), 0)
KD_PName <- names(KD_P)
KDUP   <- round(cellStats(KD_P_Crop  , 'max'), 0)
KDLO   <- round(cellStats(KD_P_Crop  , 'min'), 0)

ET_P <- raster(Mean_CLIM[14])
names(ET_P) = "Evaporation"
ET_P_Mask <- mask(x = ET_P, mask = HALE)
ET_P_Crop <- crop(x = ET_P_Mask, y = extent(HALE))
if (RFUnit == " in") { ET_P_Crop <- ET_P_Crop * 0.0393701 }
ET_P_M   <- round(cellStats(ET_P_Crop , 'mean'), 0)
ET_PName <- names(ET_P)
ETUP   <- round(cellStats(ET_P_Crop  , 'max'), 0)
ETLO   <- round(cellStats(ET_P_Crop  , 'min'), 0)

VPD_P <- raster(Mean_CLIM[105])
names(VPD_P) = "VPD"
VPD_P_Mask <- mask(x = VPD_P, mask = HALE)
VPD_P_Crop <- crop(x = VPD_P_Mask, y = extent(HALE))
if (RFUnit == " in") { VPD_P_Crop <- VPD_P_Crop * 0.0393701 }
VPD_P_M   <- round(cellStats(VPD_P_Crop, 'mean'), 0)
VPD_PName <- names(VPD_P)
VPDUP   <- round(cellStats(VPD_P_Crop  , 'max'), 0)
VPDLO   <- round(cellStats(VPD_P_Crop  , 'min'), 0)

WS_P <- raster(paste0(INPUTS_FOLDER, "Mean Climate/Wind/wind_sd_avg.tif"))
names(WS_P) = "Windspeed"
WS_P_Mask <- mask(x = WS_P, mask = HALE)
WS_P_Crop <- crop(x = WS_P_Mask, y = extent(HALE))
#convert from m/sec to mph
WS_P_Crop <- WS_P_Crop * 2.237
WS_P_M   <- round(cellStats(WS_P_Crop, 'mean'), 0)
WS_PName <- names(WS_P)
WSUP   <- round(cellStats(WS_P_Crop  , 'max'), 1)
WSLO   <- round(cellStats(WS_P_Crop  , 'min'), 1)

##########   MEAN RAINFALL Data
debug_print("51, MEAN RAINFALL Data")

Jan <- raster(MeanRF_ALL[1])
Jan_Mask <- mask(x = Jan, mask = HALE)
Jan_CropRF <- crop(x = Jan_Mask, y = extent(HALE))
if (RFUnit == " in") { Jan_CropRF  <- Jan_CropRF  * 0.0393701 }
JanMRF    <- round(cellStats(Jan_CropRF, 'mean'), 1)
JanMRFx   <- round(cellStats(Jan_CropRF  , 'max'), 1)
JanMRFn   <- round(cellStats(Jan_CropRF  , 'min'), 1)

Feb <- raster(MeanRF_ALL[2])
Feb_Mask <- mask(x = Feb, mask = HALE)
Feb_CropRF <- crop(x = Feb_Mask, y = extent(HALE))
if (RFUnit == " in") { Feb_CropRF  <- Feb_CropRF  * 0.0393701 }
FebMRF   <- round(cellStats(Feb_CropRF, 'mean'), 1)
FebMRFx   <- round(cellStats(Feb_CropRF  , 'max'), 1)
FebMRFn   <- round(cellStats(Feb_CropRF  , 'min'), 1)

Mar <- raster(MeanRF_ALL[3])
Mar_Mask <- mask(x = Mar, mask = HALE)
Mar_CropRF <- crop(x = Mar_Mask, y = extent(HALE))
if (RFUnit == " in") { Mar_CropRF  <- Mar_CropRF  * 0.0393701 }
MarMRF   <- round(cellStats(Mar_CropRF, 'mean'), 1)
MarMRFx   <- round(cellStats(Mar_CropRF  , 'max'), 1)
MarMRFn   <- round(cellStats(Mar_CropRF  , 'min'), 1)

Apr <- raster(MeanRF_ALL[4])
Apr_Mask <- mask(x = Apr, mask = HALE)
Apr_CropRF <- crop(x = Apr_Mask, y = extent(HALE))
if (RFUnit == " in") { Apr_CropRF  <- Apr_CropRF  * 0.0393701 }
AprMRF   <- round(cellStats(Apr_CropRF, 'mean'), 1)
AprMRFx   <- round(cellStats(Apr_CropRF  , 'max'), 1)
AprMRFn   <- round(cellStats(Apr_CropRF  , 'min'), 1)

May <- raster(MeanRF_ALL[5])
May_Mask <- mask(x = May, mask = HALE)
May_CropRF <- crop(x = May_Mask, y = extent(HALE))
if (RFUnit == " in") { May_CropRF  <- May_CropRF  * 0.0393701 }
MayMRF   <- round(cellStats(May_CropRF, 'mean'), 1)
MayMRFx   <- round(cellStats(May_CropRF  , 'max'), 1)
MayMRFn   <- round(cellStats(May_CropRF  , 'min'), 1)

Jun <- raster(MeanRF_ALL[6])
Jun_Mask <- mask(x = Jun, mask = HALE)
Jun_CropRF <- crop(x = Jun_Mask, y = extent(HALE))
if (RFUnit == " in") { Jun_CropRF  <- Jun_CropRF  * 0.0393701 }
JunMRF   <- round(cellStats(Jun_CropRF, 'mean'), 1)
JunMRFx   <- round(cellStats(Jun_CropRF  , 'max'), 1)
JunMRFn   <- round(cellStats(Jun_CropRF  , 'min'), 1)

Jul <- raster(MeanRF_ALL[7])
Jul_Mask <- mask(x = Jul, mask = HALE)
Jul_CropRF <- crop(x = Jul_Mask, y = extent(HALE))
if (RFUnit == " in") { Jul_CropRF  <- Jul_CropRF  * 0.0393701 }
JulMRF   <- round(cellStats(Jul_CropRF, 'mean'), 1)
JulMRFx   <- round(cellStats(Jul_CropRF  , 'max'), 1)
JulMRFn   <- round(cellStats(Jul_CropRF  , 'min'), 1)

Aug <- raster(MeanRF_ALL[8])
Aug_Mask <- mask(x = Aug, mask = HALE)
Aug_CropRF <- crop(x = Aug_Mask, y = extent(HALE))
if (RFUnit == " in") { Aug_CropRF  <- Aug_CropRF  * 0.0393701 }
AugMRF   <- round(cellStats(Aug_CropRF, 'mean'), 1)
AugMRFx   <- round(cellStats(Aug_CropRF  , 'max'), 1)
AugMRFn   <- round(cellStats(Aug_CropRF  , 'min'), 1)

Sep <- raster(MeanRF_ALL[9])
Sep_Mask <- mask(x = Sep, mask = HALE)
Sep_CropRF <- crop(x = Sep_Mask, y = extent(HALE))
if (RFUnit == " in") { Sep_CropRF  <- Sep_CropRF  * 0.0393701 }
SepMRF   <- round(cellStats(Sep_CropRF, 'mean'), 1)
SepMRFx   <- round(cellStats(Sep_CropRF  , 'max'), 1)
SepMRFn   <- round(cellStats(Sep_CropRF  , 'min'), 1)

Oct <- raster(MeanRF_ALL[10])
Oct_Mask <- mask(x = Oct, mask = HALE)
Oct_CropRF <- crop(x = Oct_Mask, y = extent(HALE))
if (RFUnit == " in") { Oct_CropRF  <- Oct_CropRF  * 0.0393701 }
OctMRF   <- round(cellStats(Oct_CropRF, 'mean'), 1)
OctMRFx   <- round(cellStats(Oct_CropRF  , 'max'), 1)
OctMRFn   <- round(cellStats(Oct_CropRF  , 'min'), 1)

Nov <- raster(MeanRF_ALL[11])
Nov_Mask <- mask(x = Nov , mask = HALE)
Nov_CropRF <- crop(x = Nov_Mask, y = extent(HALE))
if (RFUnit == " in") { Nov_CropRF  <- Nov_CropRF  * 0.0393701 }
NovMRF   <- round(cellStats(Nov_CropRF, 'mean'), 1)
NovMRFx   <- round(cellStats(Nov_CropRF  , 'max'), 1)
NovMRFn   <- round(cellStats(Nov_CropRF  , 'min'), 1)

Dec <- raster(MeanRF_ALL[12])
Dec_Mask <- mask(x = Dec, mask = HALE)
Dec_CropRF <- crop(x = Dec_Mask, y = extent(HALE))
if (RFUnit == " in") { Dec_CropRF  <- Dec_CropRF  * 0.0393701 }
DecMRF   <- round(cellStats(Dec_CropRF, 'mean'), 1)
DecMRFx   <- round(cellStats(Dec_CropRF  , 'max'), 1)
DecMRFn   <- round(cellStats(Dec_CropRF  , 'min'), 1)

ANN <- raster(MeanRF_ALL[13])
ANN_Mask <- mask(x = ANN, mask = HALE)
ANN_CropRF <- crop(x = ANN_Mask, y = extent(HALE))
if (RFUnit == " in") { ANN_CropRF  <- ANN_CropRF  * 0.0393701 }
ANNMRFM   <- round(cellStats(ANN_CropRF, 'mean'), 0)
ANNUP   <- round(cellStats(ANN_CropRF  , 'max'), 0)
ANNLO   <- round(cellStats(ANN_CropRF  , 'min'), 0)

##########   Add to Table

MEANRF2 <- c(JanMRF,FebMRF,MarMRF,AprMRF,MayMRF,JunMRF,JulMRF,AugMRF,SepMRF,OctMRF,NovMRF,DecMRF,ANNMRFM)
Cell.CL_Year[1,2:14] <- MEANRF2  

########## Get Threshold for figures

RFUP <- max(JanMRFx,FebMRFx,MarMRFx,AprMRFx,MayMRFx,JunMRFx,JulMRFx,AugMRFx,SepMRFx,OctMRFx,NovMRFx,DecMRFx)
RFLO <- min(JanMRFn,FebMRFn,MarMRFn,AprMRFn,MayMRFn,JunMRFn,JulMRFn,AugMRFn,SepMRFn,OctMRFn,NovMRFn,DecMRFn)

########## Write to Output Folder

write.csv(
  Cell.CL_Year,
  paste0(PATH_WITH_PROJECT_NAME, "Annual Climate.csv"),
  row.names = F
)
Cell.CL_Year

##########   CLIMO GRAPHS
debug_print("52, CLIMO GRAPHS")

#build data frame with temperature and precipitation data
df <- as.data.frame(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
colnames(df) <- c("month")
df$month <- factor(df$month, levels = month.abb)
PREC <- as.numeric(as.vector(Cell.CL_Year[1, 2:13]))
TA <- as.numeric(as.vector(Cell.CL_Year[3, 2:13]))

df$PREC <- PREC
df$TA <- TA

df

#Set X and Y Limits for Graph 20% higher than the max.
XPR <- max(PREC) + (max(PREC) * 0.2)
XTA <- max(TA) +   (max(TA) * 0.2)
NTA <- min(TA) - min(TA) * 0.2

#Make Climo Graph and Save to output folder

png(
  paste0(PATH_WITH_PROJECT_NAME, "Climograph2.png"),
  width = 5 * dpi,
  height = 3 * dpi,
  res = dpi
)

par(mar = c(4, 4, 4, 4))
my_bar <- barplot(df$PREC, border=F , names.arg=df$month ,
  las=2 ,
  col="darkblue" ,
  ylim=c(0,XPR) ,
  ylab = paste0("Rainfall [",RFUnit2,"]"))#,

title(paste0("Monthly Climate: ", UNIT_Ns[u]),
  line = 0.8,
  cex.main = 1.5)
# text(my_bar, df$PREC+13 ,df$PREC ,cex=1)
par(new = TRUE)
plot(df$month,df$TA ,pch=15 , lty = 0, axes = FALSE, xlab = "", ylab = "",col="red",ylim=c(NTA,XTA ))
lines(x = df$month, y = df$TA,lwd=2,col="red")
points(x = df$month, y = df$TA,col="red",pch=16,cex=2)
mtext(paste0("Temperature [",TUnit,"]"), side=4, line=2.5,col="red")
axis(4, ylim=c(NTA,XTA), col="red",col.axis="red",las=1)

# # Legend
# legend("topleft", legend = c("Rainfall", "Temperature"),
#        col = c("darkblue" , "red") ,
#        bty = "n", pch=c(15,20) , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.05, 0))

dev.off()

### Monthly Rainfall
debug_print("53, Monthly Rainfall")

png(
  paste0(PATH_WITH_PROJECT_NAME, "Climograph_RF.png"),
  width = 5 * dpi,
  height = 3 * dpi,
  res = dpi
)

barplot(df$PREC, border=F , names.arg=df$month , 
  las=2 , 
  col="darkblue" , 
  ylim=c(0,XPR) ,
  ylab = paste0("Rainfall (",RFUnit2,".)"))#,
title(paste0("Monthly Rainfall: ",UNIT_Ns[u]), line = 0.8,cex.main = 1.5)
par(new = TRUE)
plot(df$month,df$TA ,pch=15 , lty = 0, axes = FALSE, xlab = "", ylab = "",col="red",ylim=c(NTA,XTA ))

dev.off()

### Monthly Air Temp
col <- rgb(0.2, 0.2, 1, alpha = 0.15)

png(
  paste0(PATH_WITH_PROJECT_NAME, "Climograph_AT.png"),
  width = 5.2 * dpi,
  height = 3 * dpi,
  res = dpi
)

par(mar = c(4, 4, 4, 4))
barplot(df$PREC, border=F , names.arg=df$month , 
  las=2 , 
  col=col , 
  alpha = 0.5 ,
  ylim=c(0,XPR) ,
  ylab = paste0("Rainfall (",RFUnit2,".)"),
  col.lab = "grey50")#,
title(paste0("Monthly Air Temperature"), line = 0,cex.main = 1.3)
par(new = TRUE)
plot(df$month,df$TA ,pch=15 , lty = 0, axes = FALSE, xlab = "", ylab = "",col="red",ylim=c(NTA,XTA ))
lines(x = df$month, y = df$TA,lwd=2,col="red")
points(x = df$month, y = df$TA,col="red",pch=16,cex=2)
mtext(paste0("Temperature (",TUnit,")"), side=4, line=2.5,col="red")
axis(4, ylim=c(NTA,XTA), col="red",col.axis="red",las=1)

dev.off()

########## Figure for Other Climate Variables
# Decide on a Break for Rainfall

#Mean Rainfall
debug_print("53, Figure for Other Climate Variables ")

# Decide on a Break for Rainfall 
RNGERF <- RFUP-RFLO 
if (RNGERF > 8.99){
  BI_brksRF<-round(seq(RFLO - 0.1, RFUP + 0.1, length = 10),0);
  colfuncRF<-colorRampPalette(brewer.pal(10,"YlGnBu"))(50)
}
if (RNGERF < 9 && RNGERF > 7.99 ){
  BI_brksRF<-round(seq(RFLO - 0.1, RFUP + 0.1, length = 9),0);
  colfuncRF<-colorRampPalette(brewer.pal(9,"YlGnBu"))(50)
}
if (RNGERF < 8 && RNGERF > 6.99 ){
  BI_brksRF<-round(seq(RFLO - 0.1, RFUP + 0.1, length = 8),0);
  colfuncRF<-colorRampPalette(brewer.pal(8,"YlGnBu"))(50)
}
if (RNGERF < 7 && RNGERF > 5.99 ){
  BI_brksRF<-round(seq(RFLO - 0.1, RFUP + 0.1, length = 7),0);
  colfuncRF<-colorRampPalette(brewer.pal(7,"YlGnBu"))(50)
}
if (RNGERF < 6 && RNGERF > 4.99 ){
  BI_brksRF<-round(seq(RFLO - 0.1, RFUP + 0.1, length = 5),0);
  colfuncRF<-colorRampPalette(brewer.pal(5,"YlGnBu"))(50)
}
if (RNGERF < 5 && RNGERF > 3.99 ){
  BI_brksRF<-round(seq(RFLO - 0.1, RFUP + 0.1, length = 5),0);
  colfuncRF<-colorRampPalette(brewer.pal(5,"YlGnBu"))(50)
}
if (RNGERF < 4 && RNGERF > 2.99 ){
  BI_brksRF<-round(seq(RFLO - 0.1, RFUP + 0.1, length = 4),1);
  colfuncRF<-colorRampPalette(brewer.pal(4,"YlGnBu"))(50)
}
if (RNGERF < 3 && RNGERF > 1.99 ){
  BI_brksRF<-round(seq(RFLO - 0.1, RFUP + 0.1, length = 3),1);
  colfuncRF<-colorRampPalette(brewer.pal(3,"YlGnBu"))(50)
}
if (RNGERF < 2 && RNGERF > 0.99 ){
  BI_brksRF<-round(seq(RFLO - 0.1, RFUP + 0.1, length = 3),1);
  colfuncRF<-colorRampPalette(brewer.pal(3,"YlGnBu"))(50)
}
if (RNGERF < 2 && RNGERF > 0 ){
  BI_brksRF<-round(seq(RFLO - 0.1, RFUP + 0.1, length = 3),1);
  colfuncRF<-colorRampPalette(brewer.pal(3,"YlGnBu"))(50)
}

#Mean Rainfall
debug_print("53, Mean Rainfall")

# Decide on a Break for Rainfall 
RNGERFA <- ANNUP-ANNLO 

if (RNGERFA >= 10){
  BI_brksRFA<-round(seq(ANNLO - 0.1, ANNUP + 0.1, length = 11),0);
  colfuncRFA <- colorRampPalette((RColorBrewer::brewer.pal(9, "YlGnBu")))(length(BI_brksRFA) - 1)
}
if (RNGERFA < 10 && RNGERFA > 8.99 ){
  BI_brksRFA<-round(seq(ANNLO-0.1, ANNUP+0.1, length = 10),0);
  colfuncRFA<-colorRampPalette((RColorBrewer::brewer.pal(9, "YlGnBu")))(length(BI_brksRFA) - 1)
}
if (RNGERFA < 9 && RNGERFA > 7.99 ){
  BI_brksRFA<-round(seq(ANNLO-0.1, ANNUP+0.1, length = 9),0);
  colfuncRFA<-colorRampPalette((RColorBrewer::brewer.pal(9, "YlGnBu")))(length(BI_brksRFA) - 1)
}
if (RNGERFA < 8 && RNGERFA > 6.99 ){
  BI_brksRFA<-round(seq(ANNLO-0.1, ANNUP+0.1, length = 8),0);
  colfuncRFA<-colorRampPalette((RColorBrewer::brewer.pal(8, "YlGnBu")))(length(BI_brksRFA) - 1)
}
if (RNGERFA < 7 && RNGERFA > 5.99 ){
  BI_brksRFA<-round(seq(ANNLO-0.1, ANNUP+0.1, length = 7),0);
  colfuncRFA<-colorRampPalette((RColorBrewer::brewer.pal(7, "YlGnBu")))(length(BI_brksRFA) - 1)
}
if (RNGERFA < 6 && RNGERFA > 4.99 ){
  BI_brksRFA<-round(seq(ANNLO-0.1, ANNUP+0.1, length = 6),0);
  colfuncRFA<-colorRampPalette((RColorBrewer::brewer.pal(6, "YlGnBu")))(length(BI_brksRFA) - 1)
}
if (RNGERFA < 5 && RNGERFA > 3.99 ){
  BI_brksRFA<-round(seq(ANNLO-0.1, ANNUP+0.1, length = 5),0);
  colfuncRFA<-colorRampPalette((RColorBrewer::brewer.pal(5, "YlGnBu")))(length(BI_brksRFA) - 1)
}
if (RNGERFA < 4 && RNGERFA > 0.99 ){
  BI_brksRFA<-round(seq(ANNLO-0.1, ANNUP+0.1, length = 3),0);
  colfuncRFA<-colorRampPalette((RColorBrewer::brewer.pal(4, "YlGnBu")))(length(BI_brksRFA) - 1)
}
if (RNGERFA < 4 && RNGERFA < 0.99 ){
  BI_brksRFA<-round(seq(ANNLO-0.1, ANNUP+0.1, length = 3),0);
  colfuncRFA<-colorRampPalette((RColorBrewer::brewer.pal(3, "YlGnBu")))(length(BI_brksRFA) - 1)
}
if (RNGERFA < 3){
  BI_brksRFA<-round(seq(ANNLO-0.1, ANNUP+0.1, length = 3),0);
  colfuncRFA<-colorRampPalette((RColorBrewer::brewer.pal(3, "YlGnBu")))(length(BI_brksRFA) - 1)
}

#Use Same Scale for TA 
if(!is.infinite(TnLO)) {
  BI_brksTA<-round(seq(TnLO-0.1, TxUP+0.1, length = 9),0)
  colfuncTA <-colorRampPalette(brewer.pal(9,"YlOrRd"))(50)
}

if(is.infinite(TnLO)) {
  BI_brksTA<-round(seq(TaLo-0.1, TaUP+0.1, length = 9),0)
  colfuncTA <-colorRampPalette(brewer.pal(9,"YlOrRd"))(50) 
}

#For RH 
RNGERH <- RHUP+0.1-RHLO-0.1 
if (RNGERH > 8){
  BI_brksRH<-round(seq(RHLO - 0.1, RHUP + 0.1, length = 9),0)
  colfuncRH <-colorRampPalette(brewer.pal(9,"BuPu"))(50)
}
# RH
if (RNGERH < 9 && RNGERH > 7.99 ){
  BI_brksRH<-round(seq(RHLO - 0.1, RHUP + 0.1, length = 8),0)
  colfuncRH <-colorRampPalette(brewer.pal(8,"BuPu"))(50)
}
if (RNGERH < 8 && RNGERH > 6.99 ){
  BI_brksRH<-round(seq(RHLO - 0.1, RHUP + 0.1, length = 7),0)
  colfuncRH <-colorRampPalette(brewer.pal(7,"BuPu"))(50)
}
if (RNGERH < 7 && RNGERH > 5.99 ){
  BI_brksRH<-round(seq(RHLO - 0.1, RHUP + 0.1, length = 6),0)
  colfuncRH <-colorRampPalette(brewer.pal(6,"BuPu"))(50)
}
if (RNGERH < 6 && RNGERH > 4.99 ){
  BI_brksRH<-round(seq(RHLO - 0.1, RHUP + 0.1, length = 5),0)
  colfuncRH <-colorRampPalette(brewer.pal(5,"BuPu"))(50)
}
if (RNGERH < 5 && RNGERH > 3.99 ){
  BI_brksRH<-round(seq(RHLO - 0.1, RHUP + 0.1, length = 4),0)
  colfuncRH <-colorRampPalette(brewer.pal(4,"BuPu"))(50)
}
if (RNGERH < 4 && RNGERH > 2.99 ){
  BI_brksRH<-round(seq(RHLO - 0.1, RHUP + 0.1, length = 4),0)
  colfuncRH <-colorRampPalette(brewer.pal(4,"BuPu"))(50)
}
if (RNGERH < 2.99){
  l<-round((RHUP-RHLO),0)
  BI_brksRH<-round(seq(RHLO - 0.1, RHUP + 0.1, length = l),0)
  colfuncRH <-colorRampPalette(brewer.pal(l,"BuPu"))(50)
}

#For SM
BI_brksSM<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
colfuncSM <-colorRampPalette(brewer.pal(9,"YlGn"))(50)

#For KD (Solar)
if((KDUP - KDLO)>=9){
  BI_brksKD<-round(seq(KDLO - 0.1, KDUP + 0.1, length = 9),0)
  colfuncKD <-colorRampPalette(brewer.pal(9,"OrRd"))(50)
}
if((KDUP - KDLO)<9){
  BI_brksKD<-round(seq(KDLO - 0.1, KDUP + 0.1, length = (KDUP +1 - KDLO -1)),0)
  colfuncKD <-colorRampPalette(brewer.pal(3,"OrRd"))(50)
}
if((KDUP - KDLO)<=1){
  KDUP<-KDLO+1
  KDLO<-KDLO-1
  BI_brksKD<-round(seq(KDLO, KDUP, length = 2),0)
  colfuncKD <-colorRampPalette(brewer.pal(3,"OrRd"))(50)
}

#For CF
BI_brksCF<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
colfuncCF <-colorRampPalette(brewer.pal(9,"PuRd"))(50)
#For ET
BI_brksET<-round(seq(ETLO-0.1, ETUP+0.1, length = 9),0)
colfuncET <-colorRampPalette(brewer.pal(9,"PuBu"))(50)
#For WS
BI_brksWS<-round(seq(WSLO, WSUP, length = 9),2)
colfuncWS <-colorRampPalette(brewer.pal(9,"Purples"))(50)

ANN_CropRF_df <- as.data.frame(ANN_CropRF, xy = TRUE, na.rm = TRUE)
brks<- round(((max(BI_brksRFA) - min(BI_brksRFA))/4),1)
colnames(ANN_CropRF_df)[3] <- "staterf_mmann"

png(
  paste0(PATH_WITH_PROJECT_NAME, "Climate_less_RF.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
ggplot() +
  # Plot the classified raster
  geom_raster(data = ANN_CropRF_df, aes(x = x, y = y, fill = staterf_mmann)) +
  scale_fill_gradientn(colors = colfuncRFA, limits = range(BI_brksRFA),
                       na.value = "transparent", breaks = seq(min(BI_brksRFA), max(BI_brksRFA), by = brks),
                       guide = guide_colorbar(barheight = unit(0.6, "npc"))) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  # Customize the plot
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 20),  # Increase legend text size
    panel.grid = element_blank(),  # Remove gridlines
    axis.text = element_blank(),   # Remove axis text (tick labels)
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(hjust = 0.5, size = 24)
  ) +
  labs(
    title = paste0("Rainfall (", ANNMRFM, RFUnit, ")"),
    fill = NULL  # Remove legend title
  )
dev.off()

Tair_P_Crop_df<-as.data.frame(Tair_P_Crop, xy = TRUE, na.rm = TRUE)
brks<- round(((max(BI_brksTA) - min(BI_brksTA))/4),1)

png(
  paste0(PATH_WITH_PROJECT_NAME, "Climate_less_TA.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
ggplot() +
  # Plot the classified raster
  geom_raster(data = Tair_P_Crop_df, aes(x = x, y = y, fill = Mean.Air.Temp.)) +
  scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                       na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                       guide = guide_colorbar(barheight = unit(0.6, "npc"))) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  # Customize the plot
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 20),  # Increase legend text size
    panel.grid = element_blank(),  # Remove gridlines
    axis.text = element_blank(),   # Remove axis text (tick labels)
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(hjust = 0.5, size = 24)
  ) +
  labs(
    title = paste0("Air Temperature (",Tair_P_M ,TUnit2,")"),
    fill = NULL  # Remove legend title
  )

dev.off()

RH_P_Crop_df<-as.data.frame(RH_P_Crop, xy = TRUE, na.rm = TRUE)
brks<- round(((max(BI_brksRH) - min(BI_brksRH))/4),1)

png(
  paste0(PATH_WITH_PROJECT_NAME, "Climate_less_RH.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
ggplot() +
  # Plot the classified raster
  geom_raster(data = RH_P_Crop_df, aes(x = x, y = y, fill = rh_ann)) +
  scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                       na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                       guide = guide_colorbar(barheight = unit(0.6, "npc"))) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  # Customize the plot
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 20),  # Increase legend text size
    panel.grid = element_blank(),  # Remove gridlines
    axis.text = element_blank(),   # Remove axis text (tick labels)
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(hjust = 0.5, size = 24)
  ) +
  labs(
    title = paste0("Relative Humidity (",RH_P_M ," %)"),
    fill = NULL  # Remove legend title
  )
dev.off()

KD_P_Crop_df<-as.data.frame(KD_P_Crop, xy = TRUE, na.rm = TRUE)
brks<- round(((max(BI_brksKD) - min(BI_brksKD))/4),1)

png(
  paste0(PATH_WITH_PROJECT_NAME, "Climate_less_SR.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
ggplot() +
  # Plot the classified raster
  geom_raster(data = KD_P_Crop_df, aes(x = x, y = y, fill = Solar.Radiation.)) +
  scale_fill_gradientn(colors = colfuncKD, limits = range(BI_brksKD),
                       na.value = "transparent", breaks = seq(min(BI_brksKD), max(BI_brksKD), by = brks),
                       guide = guide_colorbar(barheight = unit(0.6, "npc"))) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  # Customize the plot
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 20),  # Increase legend text size
    panel.grid = element_blank(),  # Remove gridlines
    axis.text = element_blank(),   # Remove axis text (tick labels)
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(hjust = 0.5, size = 22)
  ) +
  labs(
    title = paste0("Solar Radiation (",KD_P_M ," W/m2)"),
    fill = NULL  # Remove legend title
  )
dev.off()

SM_P_Crop_df<-as.data.frame(SM_P_Crop, xy = TRUE, na.rm = TRUE)
brks<- round(((max(BI_brksSM) - min(BI_brksSM))/4),1)
png(
  paste0(PATH_WITH_PROJECT_NAME, "Climate_less_SM.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
ggplot() +
  # Plot the classified raster
  geom_raster(data = SM_P_Crop_df, aes(x = x, y = y, fill = Soil.Moisture)) +
  scale_fill_gradientn(colors = colfuncSM, limits = range(BI_brksSM),
                       na.value = "transparent", breaks = seq(min(BI_brksSM), max(BI_brksSM), by = brks),
                       guide = guide_colorbar(barheight = unit(0.6, "npc"))) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  # Customize the plot
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 20),  # Increase legend text size
    panel.grid = element_blank(),  # Remove gridlines
    axis.text = element_blank(),   # Remove axis text (tick labels)
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(hjust = 0.5, size = 24)
  ) +
  labs(
    title = paste0("Soil Moisture (",SM_P_M ,")"),
    fill = NULL  # Remove legend title
  )
dev.off()

ET_P_Crop_df<-as.data.frame(ET_P_Crop, xy = TRUE, na.rm = TRUE)
brks<- round(((max(BI_brksET) - min(BI_brksET))/4),1)
png(
  paste0(PATH_WITH_PROJECT_NAME, "Climate_less_ET.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
ggplot() +
  # Plot the classified raster
  geom_raster(data = ET_P_Crop_df, aes(x = x, y = y, fill = Evaporation)) +
  scale_fill_gradientn(colors = colfuncET, limits = range(BI_brksET),
                       na.value = "transparent", breaks = seq(min(BI_brksET), max(BI_brksET), by = brks),
                       guide = guide_colorbar(barheight = unit(0.6, "npc"))) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  # Customize the plot
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 20),  # Increase legend text size
    panel.grid = element_blank(),  # Remove gridlines
    axis.text = element_blank(),   # Remove axis text (tick labels)
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(hjust = 0.5, size = 18)
  ) +
  labs(
    title = paste0("Evapotranspiration (",ET_P_M ,RFUnit,")"),
    fill = NULL  # Remove legend title
  )
dev.off()

WS_P_Crop_df<-as.data.frame(WS_P_Crop, xy = TRUE, na.rm = TRUE)
brks<- round(((max(BI_brksWS) - min(BI_brksWS))/4),1)
if(brks == 0){brks <- 1}

png(
  paste0(PATH_WITH_PROJECT_NAME, "Climate_less_WS.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
ggplot() +
  # Plot the classified raster
  geom_raster(data = WS_P_Crop_df, aes(x = x, y = y, fill = Windspeed)) +
  scale_fill_gradientn(colors = colfuncWS, limits = range(BI_brksWS),
                       na.value = "transparent", breaks = seq(min(BI_brksWS), max(BI_brksWS), by = brks),
                       guide = guide_colorbar(barheight = unit(0.6, "npc"))) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  # Customize the plot
  coord_sf() +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 20),  # Increase legend text size
    panel.grid = element_blank(),  # Remove gridlines
    axis.text = element_blank(),   # Remove axis text (tick labels)
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.title = element_blank(),  # Remove axis titles
    plot.title = element_text(hjust = 0.5, size = 24)
  ) +
  labs(
    title = paste0("Windspeed (",WS_P_M," mph)"),
    fill = NULL  # Remove legend title
  )
dev.off()

##########  Temperature Maps Figure
debug_print("54, Temperature Maps Figure ")

Jan_CropTA_df<-as.data.frame(Jan_CropTA, xy = TRUE, na.rm = TRUE)
Feb_CropTA_df<-as.data.frame(Feb_CropTA, xy = TRUE, na.rm = TRUE)
Mar_CropTA_df<-as.data.frame(Mar_CropTA, xy = TRUE, na.rm = TRUE)
Apr_CropTA_df<-as.data.frame(Apr_CropTA, xy = TRUE, na.rm = TRUE)
May_CropTA_df<-as.data.frame(May_CropTA, xy = TRUE, na.rm = TRUE)
Jun_CropTA_df<-as.data.frame(Jun_CropTA, xy = TRUE, na.rm = TRUE)
Jul_CropTA_df<-as.data.frame(Jul_CropTA, xy = TRUE, na.rm = TRUE)
Aug_CropTA_df<-as.data.frame(Aug_CropTA, xy = TRUE, na.rm = TRUE)
Sep_CropTA_df<-as.data.frame(Sep_CropTA, xy = TRUE, na.rm = TRUE)
Oct_CropTA_df<-as.data.frame(Oct_CropTA, xy = TRUE, na.rm = TRUE)
Nov_CropTA_df<-as.data.frame(Nov_CropTA, xy = TRUE, na.rm = TRUE)
Dec_CropTA_df<-as.data.frame(Dec_CropTA, xy = TRUE, na.rm = TRUE)

brks<- round(((max(BI_brksTA) - min(BI_brksTA))/2),1)

  png(
    paste0(PATH_WITH_PROJECT_NAME, "TA12.png"),
    width = 5 * dpi,
    height = 5 * dpi,
    res = dpi
  )
  
  title1=textGrob(paste("Monthly Temperature:", UNIT_Ns[u]),gp=gpar(col="darkred",fontface="bold",fontsize=15))  
  grid.arrange(top = title1,
   #JAN
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Jan_CropTA_df, aes(x = x, y = y, fill = tair_jan)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("JAN (",JanMTA,TUnit2,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #FEB
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Feb_CropTA_df, aes(x = x, y = y, fill = tair_feb)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("FEB (",FebMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #MAR
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Mar_CropTA_df, aes(x = x, y = y, fill = tair_mar)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("MAR (",MarMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #APR
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Apr_CropTA_df, aes(x = x, y = y, fill = tair_apr)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("APR (",AprMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #MAY
   ggplot() +
     # Plot the classified raster
     geom_raster(data = May_CropTA_df, aes(x = x, y = y, fill = tair_may)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("MAY (",MayMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #JUN
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Jun_CropTA_df, aes(x = x, y = y, fill = tair_jun)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("JUN (",JunMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #JUL
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Jul_CropTA_df, aes(x = x, y = y, fill = tair_jul)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("JUL (",JulMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #AUG
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Aug_CropTA_df, aes(x = x, y = y, fill = tair_aug)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("AUG (",AugMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #SEP
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Sep_CropTA_df, aes(x = x, y = y, fill = tair_sep)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("SEP (",SepMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #OCT
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Oct_CropTA_df, aes(x = x, y = y, fill = tair_oct)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("OCT (",OctMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #NOV
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Nov_CropTA_df, aes(x = x, y = y, fill = tair_nov)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("NOV (",NovMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #DEC
   ggplot() +
     # Plot the classified raster
     geom_raster(data = Dec_CropTA_df, aes(x = x, y = y, fill = tair_dec)) +
     scale_fill_gradientn(colors = colfuncTA, limits = range(BI_brksTA),
                          na.value = "transparent", breaks = seq(min(BI_brksTA), max(BI_brksTA), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(
       title = paste0("DEC (",DecMTA,TUnit2,")"),
       fill = NULL  # Remove legend title
     ) +
     # Add bounding box around the plot
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1))
  
  dev.off()

##### Hottest and Coldest month maps
debug_print("55, find hottest and coldest months")

# find hottest and coldest months
Cell.CL_Year
t <- Cell.CL_Year[3, 2:13]
min.col <- function(m, ...) max.col(-m, ...)
hm<-names(t)[which.max(t[1,])]
cm<-names(t)[which.min(t[1,])]

# Select correct month map based on dry and wet months
if(hm == "JAN") {hmap<-Jan_CropTA ; mnh<-"January"}
if(hm == "FEB") {hmap<-Feb_CropTA ; mnh<-"February"}
if(hm == "MAR") {hmap<-Mar_CropTA ; mnh<-"March"}
if(hm == "APR") {hmap<-Apr_CropTA ; mnh<-"April"}
if(hm == "MAY") {hmap<-May_CropTA ; mnh<-"May"}
if(hm == "JUN") {hmap<-Jun_CropTA ; mnh<-"June"}
if(hm == "JUL") {hmap<-Jul_CropTA ; mnh<-"July"}
if(hm == "AUG") {hmap<-Aug_CropTA ; mnh<-"August"}
if(hm == "SEP") {hmap<-Sep_CropTA ; mnh<-"September"}
if(hm == "OCT") {hmap<-Oct_CropTA ; mnh<-"October"}
if(hm == "NOV") {hmap<-Nov_CropTA ; mnh<-"November"}
if(hm == "DEC") {hmap<-Dec_CropTA ; mnh<-"December"}

if(cm == "JAN") {cmap<-Jan_CropTA ; mnc<-"January"}
if(cm == "FEB") {cmap<-Feb_CropTA ; mnc<-"February"}
if(cm == "MAR") {cmap<-Mar_CropTA ; mnc<-"March"}
if(cm == "APR") {cmap<-Apr_CropTA ; mnc<-"April"}
if(cm == "MAY") {cmap<-May_CropTA ; mnc<-"May"}
if(cm == "JUN") {cmap<-Jun_CropTA ; mnc<-"June"}
if(cm == "JUL") {cmap<-Jul_CropTA ; mnc<-"July"}
if(cm == "AUG") {cmap<-Aug_CropTA ; mnc<-"August"}
if(cm == "SEP") {cmap<-Sep_CropTA ; mnc<-"September"}
if(cm == "OCT") {cmap<-Oct_CropTA ; mnc<-"October"}
if(cm == "NOV") {cmap<-Nov_CropTA ; mnc<-"November"}
if(cm == "DEC") {cmap<-Dec_CropTA ; mnc<-"December"}

### make plots
hmap_df<-as.data.frame(hmap, xy = TRUE, na.rm = TRUE)
colnames(hmap_df)[3] <- "val"

brks<- round(((max(BI_brksTA) - min(BI_brksTA))/4),1)

png(
  paste0(PATH_WITH_PROJECT_NAME, "TA_hm.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
ggplot() +
  # Plot the classified raster
  geom_raster(data = hmap_df, aes(x = x, y = y, fill = val)) +
  scale_fill_gradientn(
    colors = colfuncTA,
    limits = range(BI_brksTA),
    na.value = "transparent",
    breaks = pretty(range(BI_brksTA), n = 5),   # <-- exactly 5 breaks
    guide = guide_colorbar(barheight = unit(0.6, "npc"))
  ) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  coord_sf() + theme_minimal() +
  theme(legend.position = "right",legend.text = element_text(size = 18), panel.grid = element_blank(),axis.text = element_blank(),  
        axis.ticks = element_blank(),  axis.title = element_blank(),  
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)) +
  labs(
    title = paste0(mnh," Temp (°F)"),
    fill = NULL  # Remove legend titl
  ) +
  # Add bounding box around the plot
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           color = "black", fill = NA, size = 1)
dev.off()

# COLDEST MONTH MAP
cmap_df<-as.data.frame(cmap, xy=TRUE, na.rm=TRUE)
colnames(cmap_df)[3] <- "val"

brks<- round(((max(BI_brksTA) - min(BI_brksTA))/4),1)

png(
  paste0(PATH_WITH_PROJECT_NAME, "TA_cm.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
ggplot() +
  # Plot the classified raster
  geom_raster(data = cmap_df, aes(x = x, y = y, fill = val)) +
  scale_fill_gradientn(
    colors = colfuncTA,
    limits = range(BI_brksTA),
    na.value = "transparent",
    breaks = pretty(range(BI_brksTA), n = 5),   # <-- exactly 5 breaks
    guide = guide_colorbar(barheight = unit(0.6, "npc"))
  ) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  coord_sf() + theme_minimal() +
  theme(legend.position = "right",legend.text = element_text(size = 18), panel.grid = element_blank(),axis.text = element_blank(),  
        axis.ticks = element_blank(),  axis.title = element_blank(),  
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)) +
  labs(
    title = paste0(mnc," Temp (°F)"),
    fill = NULL  # Remove legend titl
  ) +
  # Add bounding box around the plot
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           color = "black", fill = NA, size = 1)
dev.off()

##########   MEAN ANNUAL Relative Humidity 12-maps
debug_print("56, MEAN ANNUAL Relative Humidity 12-maps")

# Vector of month names
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Use lapply to process each month dynamically
CropRH_dfs <- lapply(months, function(m) {
  as.data.frame(get(paste0(m, "_CropRH")), xy = TRUE, na.rm = TRUE)
})

# Assign names to the list elements
names(CropRH_dfs) <- paste0(months, "_CropRH_df")

# Rename third column to "val"
CropRH_dfs <- lapply(CropRH_dfs, function(df) {
  colnames(df)[3] <- "val"  # Rename the third column to "val"
  return(df)
})

# If you still want individual dataframes in the environment:
list2env(CropRH_dfs, envir = .GlobalEnv)

brks<- round(((max(BI_brksRH) - min(BI_brksRH))/4),1)

  png(
    paste0(PATH_WITH_PROJECT_NAME, "RH12.png"),
    width = 5 * dpi,
    height = 5 * dpi,
    res = dpi
  )
  
  title1=textGrob(paste("Monthly Relative Humidity:", UNIT_Ns[u]),gp=gpar(col="darkred",fontface="bold",fontsize=15))   
  grid.arrange(top = title1,
   #JAN             
   ggplot() +
     geom_raster(data = Jan_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("JAN (",JanMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #FEB
   ggplot() +
     geom_raster(data = Feb_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("Feb (",FebMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #MAR            
   ggplot() +
     geom_raster(data = Mar_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("MAR (",MarMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #APR        
   ggplot() +
     geom_raster(data = Apr_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("APR (",AprMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #MAY
   ggplot() +
     geom_raster(data = May_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("MAY (",MayMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #JUN             
   ggplot() +
     geom_raster(data = Jun_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("JUN (",JunMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #JUL         
   ggplot() +
     geom_raster(data = Jul_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("JUL (",JulMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #AUG          
   ggplot() +
     geom_raster(data = Aug_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("AUG (",AugMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #SEP            
   ggplot() +
     geom_raster(data = Sep_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("SEP (",SepMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #OCT            
   ggplot() +
     geom_raster(data = Oct_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("OCT (",OctMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #NOV          
   ggplot() +
     geom_raster(data = Nov_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("NOV (",NovMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #DEC          
   ggplot() +
     geom_raster(data = Dec_CropRH_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRH, limits = range(BI_brksRH),
                          na.value = "transparent", breaks = seq(min(BI_brksRH), max(BI_brksRH), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("DEC (",DecMRH ,"%)"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1))
  dev.off()

dpi = 300
##########   MEAN ANNUAL Rainfall 12-maps
debug_print("57, MEAN ANNUAL Rainfall 12-maps ")

# Vector of month names
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Use lapply to process each month dynamically
CropRF_dfs <- lapply(months, function(m) {
  as.data.frame(get(paste0(m, "_CropRF")), xy = TRUE, na.rm = TRUE)
})

# Assign names to the list elements
names(CropRF_dfs) <- paste0(months, "_CropRF_df")

# Rename third column to "val"
CropRF_dfs <- lapply(CropRF_dfs, function(df) {
  colnames(df)[3] <- "val"  # Rename the third column to "val"
  return(df)
})

# If you still want individual dataframes in the environment:
list2env(CropRF_dfs, envir = .GlobalEnv)

brks<- round(((max(BI_brksRF) - min(BI_brksRF))/4),1)

debug_print(paste0("length(values(Jan_CropRF)): ", length(values(Jan_CropRF))))

  png(
    paste0(PATH_WITH_PROJECT_NAME, "RF12.png"),
    width = 5 * dpi,
    height = 5 * dpi,
    res = dpi
  )
  
  title1=textGrob(paste("Monthly Rainfall:", UNIT_Ns[u]),gp=gpar(col="darkred",fontface="bold",fontsize=15))   
  grid.arrange(top = title1,
   #JAN
   ggplot() +
     geom_raster(data = Jan_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("JAN (",JanMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #FEB
   ggplot() +
     geom_raster(data = Feb_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("FEB (",FebMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #MAR             
   ggplot() +
     geom_raster(data = Mar_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("MAR (",MarMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #APR            
   ggplot() +
     geom_raster(data = Apr_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("APR (",AprMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #MAY             
   ggplot() +
     geom_raster(data = May_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("MAY (",MayMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #JUN            
   ggplot() +
     geom_raster(data = Jun_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("JUN (",JunMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #JUL            
   ggplot() +
     geom_raster(data = Jul_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("JUL (",JulMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #AUG            
   ggplot() +
     geom_raster(data = Aug_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("AUG (",AugMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #SEP             
   ggplot() +
     geom_raster(data = Sep_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("SEP (",SepMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #OCT            
   ggplot() +
     geom_raster(data = Oct_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("OCT (",OctMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #NOV           
   ggplot() +
     geom_raster(data = Nov_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("NOV (",NovMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1),
   #DEC            
   ggplot() +
     geom_raster(data = Dec_CropRF_df, aes(x = x, y = y, fill = val)) +
     scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                          na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                          guide = guide_colorbar(barheight = unit(0.12, "npc"))) +
     geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
     coord_sf() + theme_minimal() +
     theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
           axis.ticks = element_blank(),  axis.title = element_blank(),  
           plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
           plot.margin = margin(10, 10, 10, 10)) +
     labs(title = paste0("DEC (",DecMRF ,RFUnit,")"),fill = NULL) +
     annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
              color = "black", fill = NA, size = 1))
  dev.off()

##### Driest and Wettest month maps
debug_print("58, Driest and Wettest month maps")

# find driest and wettest months
Cell.CL_Year
d <- Cell.CL_Year[1, 2:13]
dm<-names(d)[which.min(d[1,])]
wm<-names(d)[which.max(d[1,])]

# Select correct month map based on dry and wet months
if(dm == "JAN") {dmap<-Jan_CropRF_df ; mn<-"January"}
if(dm == "FEB") {dmap<-Feb_CropRF_df ; mn<-"February"}
if(dm == "MAR") {dmap<-Mar_CropRF_df ; mn<-"March"}
if(dm == "APR") {dmap<-Apr_CropRF_df ; mn<-"April"}
if(dm == "MAY") {dmap<-May_CropRF_df ; mn<-"May"}
if(dm == "JUN") {dmap<-Jun_CropRF_df ; mn<-"June"}
if(dm == "JUL") {dmap<-Jul_CropRF_df ; mn<-"July"}
if(dm == "AUG") {dmap<-Aug_CropRF_df ; mn<-"August"}
if(dm == "SEP") {dmap<-Sep_CropRF_df ; mn<-"September"}
if(dm == "OCT") {dmap<-Oct_CropRF_df ; mn<-"October"}
if(dm == "NOV") {dmap<-Nov_CropRF_df ; mn<-"November"}
if(dm == "DEC") {dmap<-Dec_CropRF_df ; mn<-"December"}

if(wm == "JAN") {wmap<-Jan_CropRF_df ; mnw<-"January"}
if(wm == "FEB") {wmap<-Feb_CropRF_df ; mnw<-"February"}
if(wm == "MAR") {wmap<-Mar_CropRF_df ; mnw<-"March"}
if(wm == "APR") {wmap<-Apr_CropRF_df ; mnw<-"April"}
if(wm == "MAY") {wmap<-May_CropRF_df ; mnw<-"May"}
if(wm == "JUN") {wmap<-Jun_CropRF_df ; mnw<-"June"}
if(wm == "JUL") {wmap<-Jul_CropRF_df ; mnw<-"July"}
if(wm == "AUG") {wmap<-Aug_CropRF_df ; mnw<-"August"}
if(wm == "SEP") {wmap<-Sep_CropRF_df ; mnw<-"September"}
if(wm == "OCT") {wmap<-Oct_CropRF_df ; mnw<-"October"}
if(wm == "NOV") {wmap<-Nov_CropRF_df ; mnw<-"November"}
if(wm == "DEC") {wmap<-Dec_CropRF_df ; mnw<-"December"}

### make plots

brks<- round(((max(BI_brksRF) - min(BI_brksRF))/4),1)

png(
  paste0(PATH_WITH_PROJECT_NAME, "RF_dm.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

ggplot() +
  geom_raster(data = dmap, aes(x = x, y = y, fill = val)) +
  scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                       na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                       guide = guide_colorbar(barheight = unit(0.6, "npc"))) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  coord_sf() + theme_minimal() +
  theme(legend.position = "right",legend.text = element_text(size = 18), panel.grid = element_blank(),axis.text = element_blank(),  
        axis.ticks = element_blank(),  axis.title = element_blank(),  
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)) +
  labs(title = paste0(mn," Rainfall (in.)"),fill = NULL) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           color = "black", fill = NA, size = 1)

dev.off()

brks<- round(((max(BI_brksRF) - min(BI_brksRF))/4),1)

png(
  paste0(PATH_WITH_PROJECT_NAME, "RF_wm.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
ggplot() +
  geom_raster(data = wmap, aes(x = x, y = y, fill = val)) +
  scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksRF),
                       na.value = "transparent", breaks = seq(min(BI_brksRF), max(BI_brksRF), by = brks),
                       guide = guide_colorbar(barheight = unit(0.6, "npc"))) +
  geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
  coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 18), panel.grid = element_blank(),axis.text = element_blank(),  
        axis.ticks = element_blank(),  axis.title = element_blank(),  
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.margin = margin(10, 10, 10, 10)) +
  labs(title = paste0(mnw," Rainfall (in.)"),fill = NULL) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
           color = "black", fill = NA, size = 1)

dev.off()

#########   Seasonal Rainfall Maps
debug_print("59, Seasonal Rainfall Maps ")

DrySeasonRF <- (May_CropRF + Jun_CropRF + Jul_CropRF + Aug_CropRF + Sep_CropRF + Oct_CropRF)
WetSeasonRF <- (Nov_CropRF + Dec_CropRF + Jan_CropRF + Feb_CropRF + Mar_CropRF + Apr_CropRF)

### Make wet and dry season average monthly rainfall maps and values
WetSM <- (WetSeasonRF / 6)
WetSMV <- round(cellStats(WetSM, 'mean'), 1)

DrySM <- (DrySeasonRF / 6)
DrySMV <- round(cellStats(DrySM, 'mean'), 1)

# get min and max values for figure scale
DryUPm   <- ceiling(cellStats(DrySM, 'max')* 10) / 10
DryLOm   <- floor(cellStats(DrySM, 'min')* 10) / 10
WetUPm   <- ceiling(cellStats(WetSM, 'max')* 10) / 10
WetLOm  <- floor(cellStats(WetSM, 'min')* 10) / 10

DryUPMOm   <- round(DryUPm / 6, 1)
DryLOMOm   <- round(DryLOm / 6, 1)
WetUPMOm   <- round(WetUPm / 6, 1)
WetLOMOm   <- round(WetLOm / 6, 1)

SEAUPm <- max(DryUPm, WetUPm)
SEALOm <- min(DryLOm, WetLOm)
BI_brksSEAm <- round(seq(SEALOm, SEAUPm, length = 10), 0)
if ((max(BI_brksSEAm) - min(BI_brksSEAm)) < 10) {
  BI_brksSEAm <- round(seq(SEALOm, SEAUPm, length = 5), 1)
}

BI_brksSEAm
# export plot below this next section

#####
DryUP   <- round(cellStats(DrySeasonRF, 'max'), 1)
DryLO   <- round(cellStats(DrySeasonRF, 'min'), 1)
WetUP   <- round(cellStats(WetSeasonRF, 'max'), 1)
WetLO   <- round(cellStats(WetSeasonRF, 'min'), 1)

DryUPMO   <- round(DryUP / 6, 1)
DryLOMO   <- round(DryLO / 6, 1)
WetUPMO   <- round(WetUP / 6, 1)
WetLOMO   <- round(WetLO / 6, 1)

SEAUP <- max(DryUP, WetUP)
SEALO <- min(DryLO, WetLO)

BI_brksSEA <- round(seq(SEALO - 0.1, SEAUP + 0.1, length = 10),0)

RNGESEA <- SEAUP-SEALO 
if (RNGESEA > 8.99){
  BI_brksSEA<-round(seq(SEALO - 0.1, SEAUP + 0.1, length = 11),1);
  colfuncRF<-colorRampPalette(brewer.pal(9,"YlGnBu"))(50)
}
if (RNGESEA < 9 && RNGESEA > 7.99 ){
  BI_brksSEA<-round(seq(SEALO - 0.1, SEAUP + 0.1, length = 9),0);
  colfuncRF<-colorRampPalette(brewer.pal(9,"YlGnBu"))(50)
}
if (RNGESEA < 8 && RNGESEA > 6.99 ){
  BI_brksSEA<-round(seq(SEALO - 0.1, SEAUP + 0.1, length = 8),0);
  colfuncRF<-colorRampPalette(brewer.pal(8,"YlGnBu"))(50)
}
if (RNGESEA < 7 && RNGESEA > 5.99 ){
  BI_brksSEA<-round(seq(SEALO - 0.1, SEAUP + 0.1, length = 7),0);
  colfuncRF<-colorRampPalette(brewer.pal(7,"YlGnBu"))(50)
}
if (RNGESEA < 6 && RNGESEA > 4.99 ){
  BI_brksSEA<-round(seq(SEALO - 0.1, SEAUP + 0.1, length = 6),0);
  colfuncRF<-colorRampPalette(brewer.pal(6,"YlGnBu"))(50)
}
if (RNGESEA < 5 && RNGESEA > 3.99 ){
  BI_brksSEA<-round(seq(SEALO - 0.1, SEAUP + 0.1, length = 5),0);
  colfuncRF<-colorRampPalette(brewer.pal(5,"YlGnBu"))(50)
}
if (RNGESEA < 4){
  BI_brksSEA<-round(seq(SEALO - 0.1, SEAUP + 0.1, length = 4),0);
  colfuncRF<-colorRampPalette(brewer.pal(4,"YlGnBu"))(50)
}

DSeaMRF   <- round(cellStats(DrySeasonRF, 'mean'), 1)
WSeaMRF   <- round(cellStats(WetSeasonRF, 'mean'), 1)
DSeaTRF   <- round(cellStats(DrySeasonRF, 'sum'), 1)

DryMEMO   <- as.character(round(DSeaMRF / 6, 1))
DryUPMO   <- as.character(round(DryUP / 6, 1))
DryLOMO   <- as.character(round(DryLO / 6, 1))
WetMEMO   <- as.character(round(WSeaMRF / 6, 1))
WetUPMO   <- as.character(round(WetUP / 6, 1))
WetLOMO   <- as.character(round(WetLO / 6, 1))

Cell.DataCLR[1, 13] <- DSeaMRF
Cell.DataCLR[2, 13] <- DryUP
Cell.DataCLR[3, 13] <- DryLO
Cell.DataCLR[1, 14] <- WSeaMRF
Cell.DataCLR[2, 14] <- WetUP
Cell.DataCLR[3, 14] <- WetLO
Cell.DataCLR

WetSeasonRF_df<-as.data.frame(WetSeasonRF, xy=TRUE, na.rm=TRUE)
DrySeasonRF_df<-as.data.frame(DrySeasonRF, xy=TRUE, na.rm=TRUE)

png(
  paste0(PATH_WITH_PROJECT_NAME, "SeaRF.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

title1=textGrob(paste("Seasonal Rainfall:", UNIT_Ns[u]),gp=gpar(col="darkred",fontface="bold",fontsize=15)) 
title1=textGrob(paste("Seasonal Rainfall:", UNIT_Ns[u]),gp=gpar(col="darkred",fontface="bold",fontsize=15)) 
grid.arrange(top = title1,
    ggplot() +
        geom_raster(data = WetSeasonRF_df, aes(x = x, y = y, fill = layer)) +
        scale_fill_gradientn(
          colors = colfuncRF,
          limits = range(BI_brksSEA),
          na.value = "transparent",
          breaks = pretty(range(BI_brksSEA), n = 5),   # <-- exactly 5 breaks
          guide = guide_colorbar(barheight = unit(4, "cm"))
          ) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 18), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Wet Season (NOV-APR) ", WSeaMRF  ,RFUnit),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1),
 
 ggplot() +
   geom_raster(data = DrySeasonRF_df, aes(x = x, y = y, fill = layer)) +
   scale_fill_gradientn(
     colors = colfuncRF,
     limits = range(BI_brksSEA),
     na.value = "transparent",
     breaks = pretty(range(BI_brksSEA), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(4, "cm"))
   ) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 18), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Dry Season (MAY-OCT) ", DSeaMRF  ,RFUnit),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1))
dev.off()

# export seasonal monthly rainfall figure from above
WetSM_df<-as.data.frame(WetSM, xy=TRUE, na.rm=TRUE)
DrySM_df<-as.data.frame(DrySM, xy=TRUE, na.rm=TRUE)

brks<- round(((max(BI_brksSEAm) - min(BI_brksSEAm))/2),1)

png(
  paste0(PATH_WITH_PROJECT_NAME, "SeaMRF.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

title1=textGrob(paste("Monthly Rainfall:", UNIT_Ns[u]),gp=gpar(col="darkred",fontface="bold",fontsize=15)) 
grid.arrange(top = title1,
 ggplot() +
   geom_raster(data = WetSM_df, aes(x = x, y = y, fill = layer)) +
   scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksSEAm),
                        na.value = "transparent", breaks = seq(min(BI_brksSEAm), max(BI_brksSEAm), by = brks),
                        guide = guide_colorbar(barheight = unit(0.2, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 18), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Wet Season (NOV-APR) ", WetSMV  ,RFUnit),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1),
 ggplot() +
   geom_raster(data = DrySM_df, aes(x = x, y = y, fill = layer)) +
   scale_fill_gradientn(colors = colfuncRF, limits = range(BI_brksSEAm),
                        na.value = "transparent", breaks = seq(min(BI_brksSEAm), max(BI_brksSEAm), by = brks),
                        guide = guide_colorbar(barheight = unit(0.2, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 18), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Dry Season (MAY-OCT) ", DrySMV  ,RFUnit),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1))
dev.off()

### calculate the seasonal avg. monthly rainfall percentile vs. whole state
stateRF
head(stateRFMd)
summary(stateRFMd)

per <- ecdf(stateRFMd$layer) 
per
summary(per)

WetSMP <- round(per(WetSMV) * 100)
DrySMP <- round(per(DrySMV) * 100)

P <- data.frame(WetSMP, DrySMP)
P

# write seasonal monthly rainfall percentiles to csv
write.csv(P,
  paste0(PATH_WITH_PROJECT_NAME, "RF percentiles.csv"),
  row.names = F)

########## Add Data To tABLE

Cell.DataCL[u,5:13] <- c(ANNMRFM ,AnnMTA, AnnMTX,AnnMTN,ANNMRH, SM_P_M,KD_P_M,ET_P_M,CF_P_M)

Cell.DataCLR[1,3:12] <- c(ANNMRFM ,AnnMTA, AnnMTX,AnnMTN,ANNMRH, SM_P_M,KD_P_M,ET_P_M,CF_P_M,WS_P_M)
Cell.DataCLR[2,3:12] <- c(ANNUP,AnnMTAx,AnnMTXx,AnnMTNx,ANNMRHx, SMUP, KDUP,ETUP,CFUP,WSUP)
Cell.DataCLR[3,3:12] <- c(ANNLO,AnnMTAn,AnnMTXn,AnnMTNn,ANNMRHn, SMLO, KDLO,ETLO,CFLO,WSLO)

Cell.RF_Year[u,3:15] <-  MEANRF2

write.csv(
  Cell.DataCLR,
  paste0(PATH_WITH_PROJECT_NAME, "Mean Climate.csv"),
  row.names = F
)
Cell.DataCLR


#######################################################################################################################





##########  Downscaling

##########  Dynamical and Statistical Downscaling
debug_print("60, Downscaling")

PWW_T <- HALE #Name Variable UNIT

##########  Create a matrix for each cell

Cell.matrix <- matrix(nrow = 36, ncol = 2)
Cell.Data_DS <- data.frame(Cell.matrix)
colnames(Cell.Data_DS) <- c("Metric", "Value")

########## Percent Change DyDs RCP 4.5 & 8.5 2100
DyDs_P_4.5_ANN <- raster(DyDs_files[1])
Dy4.5A_M <- mask(x = DyDs_P_4.5_ANN, mask = PWW_T)
Dy4.5A_C <- crop(x = Dy4.5A_M, y = extent(PWW_T))

DyDs_P_4.5_Dry <- raster(DyDs_files[2])
Dy4.5D_M <- mask(x = DyDs_P_4.5_Dry, mask = PWW_T)
Dy4.5D_C <- crop(x = Dy4.5D_M, y = extent(PWW_T))

DyDs_P_4.5_Wet <- raster(DyDs_files[3])
Dy4.5W_M <- mask(x = DyDs_P_4.5_Wet, mask = PWW_T)
Dy4.5W_C <- crop(x = Dy4.5W_M, y = extent(PWW_T))

DyDs_P_8.5_ANN <- raster(DyDs_files[4])
Dy8.5A_M <- mask(x = DyDs_P_8.5_ANN, mask = PWW_T)
Dy8.5A_C <- crop(x = Dy8.5A_M, y = extent(PWW_T))

DyDs_P_8.5_Dry <- raster(DyDs_files[5])
Dy8.5D_M <- mask(x = DyDs_P_8.5_Dry, mask = PWW_T)
Dy8.5D_C <- crop(x = Dy8.5D_M, y = extent(PWW_T))

DyDs_P_8.5_Wet <- raster(DyDs_files[6])
Dy8.5W_M <- mask(x = DyDs_P_8.5_Wet, mask = PWW_T)
Dy8.5W_C <- crop(x = Dy8.5W_M, y = extent(PWW_T))

Dy4.5A_MM <- round(cellStats(Dy4.5A_C , 'mean'), 0)
Cell.Data_DS[1, 1:2] <- c("Dy4.5A_MM", Dy4.5A_MM)

Dy4.5D_MM <- round(cellStats(Dy4.5D_C  , 'mean'), 0)
Cell.Data_DS[2, 1:2] <- c("Dy4.5D_MM", Dy4.5D_MM)

Dy4.5W_MM <- round(cellStats(Dy4.5W_C , 'mean'), 0)
Cell.Data_DS[3, 1:2] <- c("Dy4.5W_MM", Dy4.5W_MM)

Dy8.5A_MM <- round(cellStats(Dy8.5A_C , 'mean'), 0)
Cell.Data_DS[4, 1:2] <- c("Dy8.5A_MM", Dy8.5A_MM)

Dy8.5D_MM <- round(cellStats(Dy8.5D_C , 'mean'), 0)
Cell.Data_DS[5, 1:2] <- c("Dy8.5D_MM", Dy8.5D_MM)

Dy8.5W_MM <- round(cellStats(Dy8.5W_C , 'mean'), 0)
Cell.Data_DS[6, 1:2] <- c("Dy8.5W_MM", Dy8.5W_MM)

Dy4.5A_x <- round(cellStats(Dy4.5A_C , 'max'), 0)
Dy4.5D_x <- round(cellStats(Dy4.5D_C , 'max'), 0)
Dy4.5W_x <- round(cellStats(Dy4.5W_C , 'max'), 0)
Dy8.5A_x <- round(cellStats(Dy8.5A_C , 'max'), 0)
Dy8.5D_x <- round(cellStats(Dy8.5D_C , 'max'), 0)
Dy8.5W_x <- round(cellStats(Dy8.5W_C , 'max'), 0)

Dy4.5A_n <- round(cellStats(Dy4.5A_C , 'min'), 0)
Dy4.5D_n <- round(cellStats(Dy4.5D_C , 'min'), 0)
Dy4.5W_n <- round(cellStats(Dy4.5W_C , 'min'), 0)
Dy8.5A_n <- round(cellStats(Dy8.5A_C , 'min'), 0)
Dy8.5D_n <- round(cellStats(Dy8.5D_C , 'min'), 0)
Dy8.5W_n <- round(cellStats(Dy8.5W_C , 'min'), 0)

DyRF100UP <- max(c(Dy4.5A_x, Dy4.5D_x, Dy4.5W_x, Dy8.5A_x, Dy8.5D_x, Dy8.5W_x))
DyRF100LO <- min(c(Dy4.5A_n, Dy4.5D_n, Dy4.5W_n, Dy8.5A_n, Dy8.5D_n, Dy8.5W_n))

########## DyDs Temperature 2100
#Percent Change DyDs RCP 4.5 & 8.5
debug_print("61, Percent Change DyDs RCP 4.5 & 8.5")

DyDs_files
DyDs_P_4.5_T_ChgF  <- raster(DyDs_files[37])
if (TUnit == "\u00B0F") {
  DyDs_P_4.5_T_ChgF  <- raster(DyDs_files[39])
}
Dy4.5A_MT <- mask(x = DyDs_P_4.5_T_ChgF, mask = PWW_T)
Dy4.5A_CT <- crop(x = Dy4.5A_MT, y = extent(PWW_T))

DyDs_P_8.5_T_ChgF  <- raster(DyDs_files[38])
if (TUnit == "\u00B0F") {
  DyDs_P_4.5_T_ChgF  <- raster(DyDs_files[40])
}
Dy8.5A_MT <- mask(x = DyDs_P_8.5_T_ChgF, mask = PWW_T)
Dy8.5A_CT <- crop(x = Dy8.5A_MT, y = extent(PWW_T))

###########   Calculate Mean
debug_print("62, Calculate Mean ")

Dy4.5A_MMT <- round(cellStats(Dy4.5A_CT , 'mean'), 1)
Cell.Data_DS[7, 1:2] <- c("Dy4.5A_MMT", Dy4.5A_MMT)
Dy8.5A_MMT <- round(cellStats(Dy8.5A_CT , 'mean'), 1)
Cell.Data_DS[8, 1:2] <- c("Dy8.5A_MMT", Dy8.5A_MMT)

Dy4.5A_x <- round(cellStats(Dy4.5A_CT , 'max'), 1)
Dy8.5A_x <- round(cellStats(Dy8.5A_CT , 'max'), 1)
Dy4.5A_n <- round(cellStats(Dy4.5A_CT , 'min'), 1)
Dy8.5A_n <- round(cellStats(Dy8.5A_CT , 'min'), 1)

##########   StDs Rainfall 2100
debug_print("63, StDs Rainfall 2100")

##########   Percent Change StDs RCP 4.5 & 8.5

StDs_P_4.5_ANN <- raster(DyDs_files[50])
St4.5A_M <- mask(x = StDs_P_4.5_ANN, mask = PWW_T)
St4.5A_C <- crop(x = St4.5A_M, y = extent(PWW_T))

StDs_P_4.5_Dry <- raster(DyDs_files[52])
St4.5D_M <- mask(x = StDs_P_4.5_Dry, mask = PWW_T)
St4.5D_C <- crop(x = St4.5D_M, y = extent(PWW_T))

StDs_P_4.5_Wet <- raster(DyDs_files[54])
St4.5W_M <- mask(x = StDs_P_4.5_Wet, mask = PWW_T)
St4.5W_C <- crop(x = St4.5W_M, y = extent(PWW_T))

StDs_P_8.5_ANN <- raster(DyDs_files[56])
St8.5A_M <- mask(x = StDs_P_8.5_ANN, mask = PWW_T)
St8.5A_C <- crop(x = St8.5A_M, y = extent(PWW_T))

StDs_P_8.5_Dry <- raster(DyDs_files[58])
St8.5D_M <- mask(x = StDs_P_8.5_Dry, mask = PWW_T)
St8.5D_C <- crop(x = St8.5D_M, y = extent(PWW_T))

StDs_P_8.5_Wet <- raster(DyDs_files[60])
St8.5W_M <- mask(x = StDs_P_8.5_Wet, mask = PWW_T)
St8.5W_C <- crop(x = St8.5W_M, y = extent(PWW_T))

St4.5A_MM <- round(cellStats(St4.5A_C , 'mean'), 0)
Cell.Data_DS[9, 1:2] <- c("St4.5A_MM", St4.5A_MM)

St4.5D_MM <- round(cellStats(St4.5D_C , 'mean'), 0)
Cell.Data_DS[10, 1:2] <- c("St4.5D_MM", St4.5D_MM)

St4.5W_MM <- round(cellStats(St4.5W_C , 'mean'), 0)
Cell.Data_DS[11, 1:2] <- c("St4.5W_MM", St4.5W_MM)

St8.5A_MM <- round(cellStats(St8.5A_C , 'mean'), 0)
Cell.Data_DS[12, 1:2] <- c("St8.5A_MM", St8.5A_MM)

St8.5D_MM <- round(cellStats(St8.5D_C , 'mean'), 0)
Cell.Data_DS[13, 1:2] <- c("St8.5D_MM", St8.5D_MM)

St8.5W_MM <- round(cellStats(St8.5W_C , 'mean'), 0)
Cell.Data_DS[14, 1:2] <- c("St8.5W_MM", St8.5W_MM)

St4.5A_x <- round(cellStats(St4.5A_C , 'max'), 0)
St4.5D_x <- round(cellStats(St4.5D_C , 'max'), 0)
St4.5W_x <- round(cellStats(St4.5W_C , 'max'), 0)
St8.5A_x <- round(cellStats(St8.5A_C , 'max'), 0)
St8.5D_x <- round(cellStats(St8.5D_C , 'max'), 0)
St8.5W_x <- round(cellStats(St8.5W_C , 'max'), 0)

St4.5A_n <- round(cellStats(St4.5A_C , 'min'), 0)
St4.5D_n <- round(cellStats(St4.5D_C , 'min'), 0)
St4.5W_n <- round(cellStats(St4.5W_C , 'min'), 0)
St8.5A_n <- round(cellStats(St8.5A_C , 'min'), 0)
St8.5D_n <- round(cellStats(St8.5D_C , 'min'), 0)
St8.5W_n <- round(cellStats(St8.5W_C , 'min'), 0)

StRF100UP <- max(c(St4.5A_x, St4.5D_x, St4.5W_x, St8.5A_x, St8.5D_x, St8.5W_x))
StRF100LO <- min(c(St4.5A_n, St4.5D_n, St4.5W_n, St8.5A_n, St8.5D_n, St8.5W_n))


##########   StDs Rainfall 2040-2070
debug_print("64, StDs Rainfall 2040-2070 ")

##########   Percent Change StDs RCP 4.5 & 8.5

StDs_P_4.5_ANN40 <- raster(DyDs_files[49])
St4.5A_M40 <- mask(x = StDs_P_4.5_ANN40 , mask = PWW_T)
St4.5A_C40 <- crop(x = St4.5A_M40, y = extent(PWW_T))

StDs_P_4.5_Dry40 <- raster(DyDs_files[51])
St4.5D_M40 <- mask(x = StDs_P_4.5_Dry40 , mask = PWW_T)
St4.5D_C40 <- crop(x = St4.5D_M40, y = extent(PWW_T))

StDs_P_4.5_Wet40 <- raster(DyDs_files[53])
St4.5W_M40 <- mask(x = StDs_P_4.5_Wet40 , mask = PWW_T)
St4.5W_C40 <- crop(x = St4.5W_M40, y = extent(PWW_T))

StDs_P_8.5_ANN40 <- raster(DyDs_files[55])
St8.5A_M40 <- mask(x = StDs_P_8.5_ANN40 , mask = PWW_T)
St8.5A_C40 <- crop(x = St8.5A_M40, y = extent(PWW_T))

StDs_P_8.5_Dry40 <- raster(DyDs_files[57])
St8.5D_M40 <- mask(x = StDs_P_8.5_Dry40 , mask = PWW_T)
St8.5D_C40 <- crop(x = St8.5D_M40, y = extent(PWW_T))

StDs_P_8.5_Wet40 <- raster(DyDs_files[59])
St8.5W_M40 <- mask(x = StDs_P_8.5_Wet40 , mask = PWW_T)
St8.5W_C40 <- crop(x = St8.5W_M40, y = extent(PWW_T))

St4.5A_MM40 <- round(cellStats(St4.5A_C40 , 'mean'), 0)
Cell.Data_DS[15, 1:2] <- c("St4.5A_MM40", St4.5A_MM40)
St4.5D_MM40 <- round(cellStats(St4.5D_C40 , 'mean'), 0)
Cell.Data_DS[16, 1:2] <- c("St4.5D_MM40", St4.5D_MM40)
St4.5W_MM40 <- round(cellStats(St4.5W_C40 , 'mean'), 0)
Cell.Data_DS[17, 1:2] <- c("St4.5W_MM40", St4.5W_MM40)
St8.5A_MM40 <- round(cellStats(St8.5A_C40 , 'mean'), 0)
Cell.Data_DS[18, 1:2] <- c("St8.5A_MM40", St8.5A_MM40)
St8.5D_MM40 <- round(cellStats(St8.5D_C40 , 'mean'), 0)
Cell.Data_DS[19, 1:2] <- c("St8.5D_MM40", St8.5D_MM40)
St8.5W_MM40 <- round(cellStats(St8.5W_C40 , 'mean'), 0)
Cell.Data_DS[20, 1:2] <- c("St8.5W_MM40", St8.5W_MM40)

##########   StDs Temperature 2100 ################################
debug_print("65, StDs Temperature 2100")

##########   Percent Change StDs RCP 4.5 & 8.5

StDs_P_4.5_T_ChgF  <- raster(DyDs_files[116])
if (TUnit == "\u00B0F") {
  StDs_P_4.5_T_ChgF  <- raster(DyDs_files[120])
}
St4.5A_MT <- mask(x = StDs_P_4.5_T_ChgF, mask = PWW_T)
St4.5A_CT <- crop(x = St4.5A_MT, y = extent(PWW_T))

StDs_P_8.5_T_ChgF  <- raster(DyDs_files[118])
if (TUnit == "\u00B0F") {
  StDs_P_8.5_T_ChgF  <- raster(DyDs_files[122])
}
St8.5A_MT <- mask(x = StDs_P_8.5_T_ChgF , mask = PWW_T)
St8.5A_CT <- crop(x = St8.5A_MT, y = extent(PWW_T))

#########   Calculate Mean
St4.5A_MMT <- round(cellStats(St4.5A_CT , 'mean'), 1)
Cell.Data_DS[21, 1:2] <- c("St4.5A_MMT", St4.5A_MMT)
St8.5A_MMT <- round(cellStats(St8.5A_CT , 'mean'), 1)
Cell.Data_DS[22, 1:2] <- c("St8.5A_MMT", St8.5A_MMT)

St4.5A_x <- round(cellStats(St4.5A_CT , 'max'), 1)
St8.5A_x <- round(cellStats(St8.5A_CT , 'max'), 1)
St4.5A_n <- round(cellStats(St4.5A_CT , 'min'), 1)
St8.5A_n <- round(cellStats(St8.5A_CT , 'min'), 1)

DsTAUP <- max(Dy4.5A_x, Dy8.5A_x, Dy4.5A_n, Dy8.5A_n, St4.5A_x, St8.5A_x, St4.5A_n, St8.5A_n)
DsTALO <- min(Dy4.5A_x, Dy8.5A_x, Dy4.5A_n, Dy8.5A_n, St4.5A_x, St8.5A_x, St4.5A_n, St8.5A_n)
DsTaup <- max(abs(DsTAUP ),abs(DsTALO ))
DsTAlo <-  DsTaup*-1

#Percent Change StDs RCP 4.5 & 8.5
StDs_P_4.5_T_ChgF40  <- raster(DyDs_files[115])
if (TUnit == "\u00B0F") {
  StDs_P_4.5_T_ChgF40  <- raster(DyDs_files[119])
}
St4.5A_MT40 <- mask(x = StDs_P_4.5_T_ChgF40 , mask = PWW_T)
St4.5A_CT40 <- crop(x = St4.5A_MT40, y = extent(PWW_T))

StDs_P_8.5_T_ChgF40  <- raster(DyDs_files[117])
if (TUnit == "\u00B0F") {
  StDs_P_8.5_T_ChgF40  <- raster(DyDs_files[121])
}
St8.5A_MT40 <- mask(x = StDs_P_8.5_T_ChgF40 , mask = PWW_T)
St8.5A_CT40 <- crop(x = St8.5A_MT40, y = extent(PWW_T))

#Calculate Mean
St4.5A_MMT40 <- round(cellStats(St4.5A_CT40 , 'mean'), 1)
Cell.Data_DS[23, 1:2] <- c("St4.5A_MMT40", St4.5A_MMT40)
St8.5A_MMT40 <- round(cellStats(St8.5A_CT40 , 'mean'), 1)
Cell.Data_DS[24, 1:2] <- c("St8.5A_MMT40", St8.5A_MMT40)


############################################
###########################################
###########################################


StRF100UP <- max(c(St4.5A_x, St4.5D_x, St4.5W_x, St8.5A_x, St8.5D_x, St8.5W_x))
StRF100LO <- min(c(St4.5A_n, St4.5D_n, St4.5W_n, St8.5A_n, St8.5D_n, St8.5W_n))
RFdsup <- max(c(abs(DyRF100UP),abs(DyRF100LO),abs(StRF100UP),abs(StRF100LO)))
RFdslo <- RFdsup * -1

##########   Write all DS results

write.csv(
  Cell.Data_DS,
  paste0(PATH_WITH_PROJECT_NAME, "Downscaling.csv"),
  row.names = F
)

##########   FIGURES
debug_print("68, FIGURES")

##########   Dynamical and Statistical Downscaling


############################################################################

#These graphs make use of North arrow and scale bar from RF Maps
colfunc2 <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
BI_brks2 <- round(seq(RFdslo, RFdsup , length = 9), 0)
brks <- 5

Dy4.5A_C_df<-as.data.frame(Dy4.5A_C, xy=TRUE, na.rm=TRUE)
colnames(Dy4.5A_C_df)[3] <- "val"
head(Dy4.5A_C_df)
St4.5A_C_df<-as.data.frame(St4.5A_C, xy=TRUE, na.rm=TRUE)
colnames(St4.5A_C_df)[3] <- "val"
Dy8.5A_C_df<-as.data.frame(Dy8.5A_C, xy=TRUE, na.rm=TRUE)
colnames(Dy8.5A_C_df)[3] <- "val"
St8.5A_C_df<-as.data.frame(St8.5A_C, xy=TRUE, na.rm=TRUE)
colnames(St8.5A_C_df)[3] <- "val"

png(
  paste0(PATH_WITH_PROJECT_NAME, "DS_RF_8.5_v2.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

title1=textGrob("Changes in Annual Rainfall by 2100", gp=gpar(col="darkred",fontface="bold",fontsize=15))  

title1=textGrob("Changes in Annual Rainfall by 2100", gp=gpar(col="darkred",fontface="bold",fontsize=15))  

grid.arrange(top = title1,
 ggplot() +
   geom_raster(data = Dy4.5A_C_df, aes(x = x, y = y, fill = val)) +
   # scale_fill_gradientn(colors = colfunc2, limits = range(BI_brks2),
   #                      na.value = "transparent", breaks = seq(min(BI_brks2), max(BI_brks2), by = brks),
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 7),   # <-- exactly 5 breaks
                        guide = guide_colorbar(barheight = unit(0.3, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Dynamical 4.5 (",Dy4.5A_MM ,"% change)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = St4.5A_C_df, aes(x = x, y = y, fill = val)) +
  scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 7),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.3, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Statistical 4.5 (",St4.5A_MM ,"% change)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = Dy8.5A_C_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 7),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.3, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Dynamical 8.5 (",Dy8.5A_MM ,"% change)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = St8.5A_C_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 7),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.3, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Statistical 8.5 (",St8.5A_MM ,"% change)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1))

dev.off()

########## Dynamical and Statistical Downscaling
debug_print("69, Dynamical and Statistical Downscaling ")

########## Downscaling Compare RF RCP 4.5

Dy4.5A_C_df<-as.data.frame(Dy4.5A_C, xy=TRUE, na.rm=TRUE)
colnames(Dy4.5A_C_df)[3]<-"val"
St4.5A_C_df<-as.data.frame(St4.5A_C, xy=TRUE, na.rm=TRUE)
colnames(St4.5A_C_df)[3]<-"val"
Dy4.5D_C_df<-as.data.frame(Dy4.5D_C, xy=TRUE, na.rm=TRUE)
colnames(Dy4.5D_C_df)[3]<-"val"
St4.5D_C_df<-as.data.frame(St4.5D_C, xy=TRUE, na.rm=TRUE)
colnames(St4.5D_C_df)[3]<-"val"
Dy4.5W_C_df<-as.data.frame(Dy4.5W_C, xy=TRUE, na.rm=TRUE)
colnames(Dy4.5W_C_df)[3]<-"val"
St4.5W_C_df<-as.data.frame(St4.5W_C, xy=TRUE, na.rm=TRUE)
colnames(St4.5W_C_df)[3]<-"val"

png(
  paste0(PATH_WITH_PROJECT_NAME, "DS_RF_2100_4.5.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

title1=textGrob("Dynamical and Statistical DS RCP 4.5, 2100", gp=gpar(col="darkred",fontface="bold",fontsize=15))  
grid.arrange(top = title1,
 ggplot() +
   geom_raster(data = Dy4.5A_C_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.2, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("DyDs ANN (",Dy4.5A_MM ,"%)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = St4.5A_C_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.2, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("StDs ANN (",St4.5A_MM ,"%)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = Dy4.5D_C_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.2, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("DyDs DRY (",Dy4.5D_MM ,"%)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = St4.5D_C_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.2, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("StDs Dry (",St4.5D_MM ,"%)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 , 
 ggplot() +
   geom_raster(data = Dy4.5W_C_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.2, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("DyDs WET (",Dy4.5W_MM ,"%)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = St4.5W_C_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.2, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 10), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("StDs WET (",St4.5W_MM ,"%)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
)              

dev.off()        

########### Statistical 2040-2070
debug_print("70, Statistical 2040-2070")

St4.5A_C40_df<-as.data.frame(St4.5A_C40, xy=TRUE, na.rm=TRUE)
colnames(St4.5A_C40_df)[3]<-"val"
St8.5A_C40_df<-as.data.frame(St8.5A_C40, xy=TRUE, na.rm=TRUE)
colnames(St8.5A_C40_df)[3]<-"val"

brks<- round(((max(BI_brks2) - min(BI_brks2))/4),1)

png(
  paste0(PATH_WITH_PROJECT_NAME, "StDsRF2040.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

title1=textGrob("Changes in Annual Rainfall by Mid-Century", gp=gpar(col="darkred",fontface="bold",fontsize=15))  

grid.arrange(top = title1,
 ggplot() +
   geom_raster(data = St4.5A_C40_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 7),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.5, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 12), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("RCP 4.5 (",St4.5A_MM40 ,"% change)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = St8.5A_C40_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc2,
     limits = range(BI_brks2),
     na.value = "transparent",
     breaks = pretty(range(BI_brks2), n = 7),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.5, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 12), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("RCP 8.5 (",St8.5A_MM40 ,"% change)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1),
 ncol=2)

dev.off()

##########   Downscaling Compare TEMP RCP 8.5 & 4.5 2100
debug_print("71, Downscaling Compare TEMP RCP 8.5 & 4.5 2100")

colfunc3 <- colorRampPalette(brewer.pal(9, "Reds"))(100)

if(TUnit == "°C") {BI_brksF<-round(seq(0, 7 , length = 7),0)}
if(TUnit == "°F") {BI_brksF<-round(seq(0, 9 , length = 9),0)}

brkst<-round(((max(BI_brksF) - min(BI_brksF))/4),1)

Dy4.5A_CT_df<-as.data.frame(Dy4.5A_CT, xy=TRUE, na.rm=TRUE)
colnames(Dy4.5A_CT_df)[3]<-"val"
St4.5A_CT_df<-as.data.frame(St4.5A_CT, xy=TRUE, na.rm=TRUE)
colnames(St4.5A_CT_df)[3]<-"val"
Dy8.5A_CT_df<-as.data.frame(Dy8.5A_CT, xy=TRUE, na.rm=TRUE)
colnames(Dy8.5A_CT_df)[3]<-"val"
St8.5A_CT_df<-as.data.frame(St8.5A_CT, xy=TRUE, na.rm=TRUE)
colnames(St8.5A_CT_df)[3]<-"val"

png(
  paste0(PATH_WITH_PROJECT_NAME, "DS_Temp2100.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

title4=textGrob("Changes in Air Temperature by 2100", gp=gpar(col="darkred",fontface="bold",fontsize=15))
grid.arrange(top = title4,
 ggplot() +
   geom_raster(data = Dy4.5A_CT_df, aes(x = x, y = y, fill = val)) +
   # scale_fill_gradientn(colors = colfunc3, limits = range(BI_brksF),
   #                      na.value = "transparent", breaks = seq(min(BI_brksF), max(BI_brksF), by = brkst),
   #                      guide = guide_colorbar(barheight = unit(0.2, "npc"))) +
   scale_fill_gradientn(
     colors = colfunc3,
     limits = range(BI_brksF),
     na.value = "transparent",
     breaks = pretty(range(BI_brksF), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.3, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 12), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Dynamical 4.5 (",Dy4.5A_MMT ,TUnit2," average)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = St4.5A_CT_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc3,
     limits = range(BI_brksF),
     na.value = "transparent",
     breaks = pretty(range(BI_brksF), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.3, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 12), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Statistical 4.5 (",St4.5A_MMT ,TUnit2," average)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = Dy8.5A_CT_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc3,
     limits = range(BI_brksF),
     na.value = "transparent",
     breaks = pretty(range(BI_brksF), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.3, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 12), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Dynamical 8.5 (",Dy8.5A_MMT ,TUnit2," average)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 ,
 ggplot() +
   geom_raster(data = St8.5A_CT_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc3,
     limits = range(BI_brksF),
     na.value = "transparent",
     breaks = pretty(range(BI_brksF), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.3, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 12), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("Statistical 8.5 (",St8.5A_MMT ,TUnit2," average)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
)

dev.off()


######################## AIR TEMP 2040-2070
debug_print("72, AIR TEMP 2040-2070 ")

St4.5A_CT40_df<-as.data.frame(St4.5A_CT40, xy=TRUE, na.rm=TRUE)
colnames(St4.5A_CT40_df)[3]<-"val"
St8.5A_CT40_df<-as.data.frame(St8.5A_CT40, xy=TRUE, na.rm=TRUE)
colnames(St8.5A_CT40_df)[3]<-"val"

png(
  paste0(PATH_WITH_PROJECT_NAME, "StDs_Temp2040.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

title4=textGrob("Changes in Air Temperature by Mid-Century", gp=gpar(col="darkred",fontface="bold",fontsize=15))

grid.arrange(top = title4,
 ggplot() +
   geom_raster(data = St4.5A_CT40_df, aes(x = x, y = y, fill = val)) +
   # scale_fill_gradientn(colors = colfunc3, limits = range(BI_brksF),
   #                      na.value = "transparent", breaks = seq(min(BI_brksF), max(BI_brksF), by = brkst),
   #                      guide = guide_colorbar(barheight = unit(0.5, "npc"))) +
   scale_fill_gradientn(
     colors = colfunc3,
     limits = range(BI_brksF),
     na.value = "transparent",
     breaks = pretty(range(BI_brksF), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.5, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 12), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("RCP 4.5 (",St4.5A_MMT40 ,TUnit2," change)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 , 
 ggplot() +
   geom_raster(data = St8.5A_CT40_df, aes(x = x, y = y, fill = val)) +
   scale_fill_gradientn(
     colors = colfunc3,
     limits = range(BI_brksF),
     na.value = "transparent",
     breaks = pretty(range(BI_brksF), n = 5),   # <-- exactly 5 breaks
     guide = guide_colorbar(barheight = unit(0.5, "npc"))) +
   geom_sf(data = HALE_sf, fill = NA, color = "black", linewidth = 1) +
   coord_sf() + theme_minimal() +
   theme(legend.position = "right",legend.text = element_text(size = 12), panel.grid = element_blank(),axis.text = element_blank(),  
         axis.ticks = element_blank(),  axis.title = element_blank(),  
         plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
         plot.margin = margin(10, 10, 10, 10)) +
   labs(title = paste0("RCP 8.5 (",St8.5A_MMT40 ,TUnit2," change)"),fill = NULL) +
   annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, 
            color = "black", fill = NA, size = 1)
 , 
 ncol=2)

dev.off()

##########################################################################################################################

##########  Create A Monthly Rainfall Time Series
debug_print("73, Rainfall Extract ")

debug_print("Rainfall Extract")
debug_print("Frazier et al 2016")

#########   Extract RF Data From Monthly Maps (2 DATASETS)
#########   Call in the Shapefile
PWW_T <- HALE

#########   Load RF MAPS
RF_Map_Path_A
RF_Tif_files = dir(RF_Map_Path_A, pattern="*.tif", recursive=T, full.names=T)  #Monthly RF
nfiles <- length(RF_Tif_files)

#Create a matrix for each cell
Cell.AF_Maps <- data.frame(matrix(nrow = 1116, ncol = 4))
colnames(Cell.AF_Maps) <- c("Date", "Year", "Month", "RF")

for (i in 1:nfiles) {
  if(i == 1){print("AF")}
  if(i == 250){print(i)} #Just to see where Im at in the process 
  if(i == 500){print(i)}
  if(i == 750){print(i)}
  if(i == 1000){print(i)}
  
  #Open Map as raster and change projection
  RF_Map <- RF_Tif_files[i]
  RF_Map2 <- raster(RF_Map)

  #Get the Name (date) from each map (derek's version)
  name <- basename(RF_Map)
  nameSplit <- unlist(strsplit(name, "_"))
  nameSplit2 <- unlist(strsplit(nameSplit[8] , ".t"))
  Year <- nameSplit[7]
  Month <- nameSplit2[1]
  MY <- paste0(Year, "/", Month)
  
  
  projection(RF_Map2) <- projection(EXAMP)
  RFM_Mask <- mask(x = RF_Map2, mask = PWW_T)
  RFM_Crop <- crop(x = RFM_Mask, y = extent(PWW_T))
  RFM_M   <- round(cellStats(RFM_Crop, 'mean'), 1)
  
  Cell.AF_Maps[i, 1:4] <- c(MY, Year, Month, RFM_M)
  
}

head(Cell.AF_Maps)
tail(Cell.AF_Maps)

##########   Matty's Maps
debug_print("74, Lucas")

#Load Daily RF MAPS
RF_Map_Path
RF_Tif_files = dir(RF_Map_Path, pattern="*.tif", recursive=T, full.names=T)  #Monthly RF
nfiles <- length(RF_Tif_files)

#Create a matrix for each cell
Cell.ML_Maps <- data.frame(matrix(nrow = 396, ncol = 4))
colnames(Cell.ML_Maps) <- c("Date", "Year", "Month", "RF")

for (i in 1:nfiles) {
  if(i == 1){print("ML")}
  if(i == 100){print(i)} #Just to see where Im at in the process 
  if(i == 200){print(i)}
  if(i == 300){print(i)}
  
  ###########   Open Map as raster and change projection
  
  RF_Map <- RF_Tif_files[i]
  RF_Map2 <- raster(RF_Map)
  plot(RF_Map2)
  
  ##########   #Get the Name (date) from each map
  # Derek's version
  name <- basename(RF_Map)
  nameSplit <- unlist(strsplit(name, "_"))
  nameSplit2 <- unlist(strsplit(nameSplit[8] , ".t"))
  D1 <- nameSplit2[1]
  Year <- nameSplit[7]
  Month <- nameSplit2[1]
  MY <- paste0(Year, "/", Month)
  
  projection(RF_Map2) <- projection(EXAMP)
  RFM_Mask <- mask(x = RF_Map2, mask = PWW_T)
  RFM_Crop <- crop(x = RFM_Mask, y = extent(PWW_T))
  RFM_M   <- round(cellStats(RFM_Crop, 'mean'), 1)
  
  ##########  ADD Data To Matrix
  
  Cell.ML_Maps[i, 1:4] <- c(MY, Year, Month, RFM_M)
  
} #End Month RF Extract

head(Cell.ML_Maps)
tail(Cell.ML_Maps)

###

### Abby's data
# subset for 1990 - 2012
MRF_A2 =  Cell.AF_Maps[c(841:1116), ]
head(MRF_A2)

# just keep RF values
MRF_A3 =  as.numeric(MRF_A2[, 4])

# convert from mm to inches
if (RFUnit == " in") { MRF_A3 <-  MRF_A3 * 0.0393701 }
summary(MRF_A3)

### Matty's data
# subet for 1990 - 2012
MRF_N2 =  Cell.ML_Maps[c(1:276), ]
head(MRF_N2)

# just keep RF values
MRF_N3 =  as.numeric(MRF_N2[, 4])

# convert mm to inches
if (RFUnit == " in") { MRF_N3  <-  MRF_N3  * 0.0393701 }
summary(MRF_N3)



MBE <- round(mean(MRF_A3 - MRF_N3), 1)
MAE <- round(mean(abs(MRF_A3 - MRF_N3)), 1)

D_Comp <- cbind(MRF_A3, MRF_N3)
Mx <- max(D_Comp, na.rm = T)

FNAME <- paste0("RF_Compare_23_", UNIT_Ns[u], ".csv")

##########   Comparison Figures For the two datasets
debug_print("75, Comparison Figures For the two datasets")

LM1 <- lm(MRF_A3 ~ MRF_N3)
LM1P <- round(coefficients(summary(LM1))[2, 4], 4)
LM1R <- round(summary(LM1)$r.squared, 2)
TITLE = paste(UNIT_Ns[u], " 23-yr RF Compare (", RFUnit2, ")")

dpi = 300

png(
  paste0(PATH_WITH_PROJECT_NAME, "23yr_RF_Compare.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

plot(MRF_A3~MRF_N3,ylim=c(0,Mx),xlim=c(0,Mx),main=TITLE,
  ylab = "Frazier et al. (2016)",xlab="Lucas et al. (2022)")
abline(0,1)
legend("topleft",c(paste("R2 = ",LM1R),paste("MBE = ",MBE),paste("MAE = ",MAE)))

dev.off()

##########   Merge Datasets to create full time period (1920 - current)
debug_print("76, Merge Datasets to create full time period (1920 - current)")

nrows <- nrow(Cell.ML_Maps)

MRF_ND3 =  Cell.ML_Maps[c(1:nrows), ]
head(MRF_ND3)
tail(MRF_ND3)

MRF_AD3 =  Cell.AF_Maps[c(1:840), ]
head(MRF_AD3)
tail(MRF_AD3)


colnames(MRF_AD3) <- c("Date", "Year", "Month", "RF")
colnames(MRF_ND3) <- c("Date", "Year", "Month", "RF")

MRF100 <- rbind(MRF_AD3, MRF_ND3)
if (RFUnit == " in") { MRF100$RF <-  as.numeric(MRF100$RF) * 0.0393701 }

write.csv(
  MRF100,
  paste0(PATH_WITH_PROJECT_NAME, "Monthly Rainfall_", RFUnit2, ".csv"),
  row.names = F
)

# ## Can read in the csv from above to start script here
# MRF100<-read.csv(paste0(OUTPUTS_FOLDER,UNIT_N[u],"/",UNIT_N[u],"_Monthly Rainfall_in.csv"))

head(MRF100)
# fix month values
MRF100$Month <- sub(".*/", "", MRF100$Date)

# get end date
ey <- as.numeric(MRF100[nrow(MRF100), ]$Year)
em <- as.numeric(MRF100[nrow(MRF100), ]$Month)
ed <- as.Date(MRF100[nrow(MRF100), ]$Ndate)

MRF100[nrow(MRF100), ]

head(MRF100)
tail(MRF100)

ggplot(data = MRF100, aes(x = Year, y = RF)) +
  geom_line()

RF_IN <- as.numeric(MRF100[, 4])
summary(RF_IN)

MRF100$Day <- 1
MRF100$Ndate <- as.Date(with(MRF100, paste(Year, Month, Day, sep = "-")), "%Y-%m-%d")

MRF_Max <- round(max(RF_IN, na.rm = T), 1)
MRF_Min <- round(min(RF_IN, na.rm = TRUE), 1)
MRF_MED <- round(median(RF_IN, na.rm = TRUE), 1)
MRF_MEAN <- round(mean(RF_IN, na.rm = TRUE), 1)

MoRF.ts <- ts(RF_IN, c(1920, 1), end = c(ey, em), frequency = 12)

myts1 <- as.vector(window(MoRF.ts, start = c(1920, 1), end = c(ey, em)))
myts2 <- as.vector(window(MoRF.ts, start = c(1940, 1), end = c(ey, em)))
myts3 <- as.vector(window(MoRF.ts, start = c(1960, 1), end = c(ey, em)))
myts4 <- as.vector(window(MoRF.ts, start = c(1980, 1), end = c(ey, em)))
myts5 <- as.vector(window(MoRF.ts, start = c(2000, 1), end = c(ey, em)))
myts6 <- as.vector(window(MoRF.ts, start = c(2010, 1), end = c(ey, em)))

ed <- as.Date(MRF100[nrow(MRF100), ]$Ndate)

DateT1 <- seq(as.Date("1920-01-01"), as.Date(ed), by = "months")
DateT2 <- seq(as.Date("1940-01-01"), as.Date(ed), by = "months")
DateT3 <- seq(as.Date("1960-01-01"), as.Date(ed), by = "months")
DateT4 <- seq(as.Date("1980-01-01"), as.Date(ed), by = "months")
DateT5 <- seq(as.Date("2000-01-01"), as.Date(ed), by = "months")
DateT6 <- seq(as.Date("2010-01-01"), as.Date(ed), by = "months")

LM1 <- lm(myts1 ~ DateT1)
LM2 <- lm(myts2 ~ DateT2)
LM3 <- lm(myts3 ~ DateT3)
LM4 <- lm(myts4 ~ DateT4)
LM5 <- lm(myts5 ~ DateT5)
LM6 <- lm(myts6 ~ DateT6)

### Get direction of trendlines for PPT
if(coef(LM1)[[2]] > 0) {LM1s <- c("Increase")}
if(coef(LM4)[[2]] > 0) {LM4s <- c("Increase")}
if(coef(LM6)[[2]] > 0) {LM6s <- c("Increase")}
if(coef(LM1)[[2]] < 0) {LM1s <- c("Decrease")}
if(coef(LM4)[[2]] < 0) {LM4s <- c("Decrease")}
if(coef(LM6)[[2]] < 0) {LM6s <- c("Decrease")}

l <- paste0("1920 - ", ey)
m <- paste0("1980 - ", ey)
s <- paste0("2010 - ", ey)

RFT <- data.frame(Period = c(l, m, s))
RFT
RFT$Trend <- c(LM1s, LM4s, LM6s)
RFT
write.csv(
  RFT,
  paste0(PATH_WITH_PROJECT_NAME, "RF_trend_directions.csv"),
  row.names = F
)
###

T1 <- round(coefficients(summary(LM1))[2, 1], 2)
T2 <- round(coefficients(summary(LM2))[2, 1], 2)
T3 <- round(coefficients(summary(LM3))[2, 1], 2)
T4 <- round(coefficients(summary(LM4))[2, 1], 2)
T5 <- round(coefficients(summary(LM5))[2, 1], 2)
T6 <- round(coefficients(summary(LM6))[2, 1], 2)

LM1P <- coefficients(summary(LM1))[2, 4]
LM2P <- coefficients(summary(LM2))[2, 4]
LM3P <- coefficients(summary(LM3))[2, 4]
LM4P <- coefficients(summary(LM4))[2, 4]
LM5P <- coefficients(summary(LM5))[2, 4]
LM6P <- coefficients(summary(LM6))[2, 4]

LM1R <- summary(LM1)$r.squared
LM2R <- summary(LM2)$r.squared
LM3R <- summary(LM3)$r.squared
LM4R <- summary(LM4)$r.squared
LM5R <- summary(LM5)$r.squared
LM6R <- summary(LM6)$r.squared

# ########## 2003-2021
#
# MA02 <- round(mean(myts6),1)
# ME02 <- round(median(myts6),1)
# MX02 <- round(max(myts6),1)
# MI02 <- round(min(myts6),1)

##########   Aggregate Month to Year

##########   Aggregate From Monthly to annual Average
debug_print("77, Aggregate From Monthly to annual Average")

rm(mean)

# this does daily to monthly. Not needed but keeping here
short.date_M = strftime(MRF100$Ndate, "%Y/%m")
Mean_M_RF = aggregate(as.numeric(MRF100$RF) ~ short.date_M, FUN = mean)
colnames(Mean_M_RF) <- c("Date", "RF")
Mean_M_RF

# monthly to annual mean rainfall
short.date_Y = strftime(as.Date(MRF100$Ndate), "%Y")
Mean_Y_RF = aggregate(as.numeric(MRF100$RF) ~ short.date_Y, FUN = mean)
colnames(Mean_Y_RF) <- c("Date", "RF")
head(Mean_Y_RF, 20)

write.csv(
  Mean_Y_RF,
  paste0(PATH_WITH_PROJECT_NAME, "Annual_RF_in.csv"),
  row.names = F
)


##########   Plot Annual RF
debug_print("78, Plot Annual RF")

YrRF.ts <- ts(Mean_Y_RF$RF, c(1920), end = c(ey), frequency = 1)
myts1Y <- as.vector(window(YrRF.ts, start = c(1920), end = c(ey)))
myts2Y <- as.vector(window(YrRF.ts, start = c(1940), end = c(ey)))
myts3Y <- as.vector(window(YrRF.ts, start = c(1960), end = c(ey)))
myts4Y <- as.vector(window(YrRF.ts, start = c(1980), end = c(ey)))
myts5Y <- as.vector(window(YrRF.ts, start = c(2000), end = c(ey)))
myts6Y <- as.vector(window(YrRF.ts, start = c(2010), end = c(ey)))

YDateT1 <- seq(as.Date("1920-01-01"), as.Date(ed), by = "years")
YDateT2 <- seq(as.Date("1940-01-01"), as.Date(ed), by = "years")
YDateT3 <- seq(as.Date("1960-01-01"), as.Date(ed), by = "years")
YDateT4 <- seq(as.Date("1980-01-01"), as.Date(ed), by = "years")
YDateT5 <- seq(as.Date("2000-01-01"), as.Date(ed), by = "years")
YDateT6 <- seq(as.Date("2010-01-01"), as.Date(ed), by = "years")

LM1Y <- lm(myts1Y ~ YDateT1)
LM2Y <- lm(myts2Y ~ YDateT2)
LM3Y <- lm(myts3Y ~ YDateT3)
LM4Y <- lm(myts4Y ~ YDateT4)
LM5Y <- lm(myts5Y ~ YDateT5)
LM6Y <- lm(myts6Y ~ YDateT6)

T1Y <- round(coefficients(summary(LM1Y))[2, 1], 3)
T1Y <- round(coefficients(summary(LM1Y))[2, 1], 3)
T2Y <- round(coefficients(summary(LM2Y))[2, 1], 3)
T3Y <- round(coefficients(summary(LM3Y))[2, 1], 3)
T4Y <- round(coefficients(summary(LM4Y))[2, 1], 3)
T5Y <- round(coefficients(summary(LM5Y))[2, 1], 3)
T6Y <- round(coefficients(summary(LM6Y))[2, 1], 3)

LM1PY <- round(coefficients(summary(LM1Y))[2, 4], 2)
LM2PY <- round(coefficients(summary(LM2Y))[2, 4], 2)
LM3PY <- round(coefficients(summary(LM3Y))[2, 4], 2)
LM4PY <- round(coefficients(summary(LM4Y))[2, 4], 2)
LM5PY <- round(coefficients(summary(LM5Y))[2, 4], 2)
LM6PY <- round(coefficients(summary(LM6Y))[2, 4], 2)

LM1RY <- round(summary(LM1Y)$r.squared, 2)
LM2RY <- round(summary(LM2Y)$r.squared, 2)
LM3RY <- round(summary(LM3Y)$r.squared, 2)
LM4RY <- round(summary(LM4Y)$r.squared, 2)
LM5RY <- round(summary(LM5Y)$r.squared, 2)
LM6RY <- round(summary(LM6Y)$r.squared, 2)

##########   Seasonal RF
debug_print("79, Seasonal RF")

# Wet Season
MRF100a <- MRF100
MRF100a$Month <- as.numeric(MRF100a$Month)
head(MRF100a)

WET_RF <- subset(MRF100a, Month == c('1', '2', '3', '4', '11', '12'))
head(WET_RF)

# Add in two rows at beginning for beginning of 1920 wet season (Nov and Dec. 1919)
WET_RF2 <- rbind(c(NA, NA, "12", NA, NA, NA), WET_RF)
WET_RF3 <- rbind(c(NA, NA, "11", NA, NA, NA), WET_RF2)
head(WET_RF3, 12)
WET_RF4 <- as.numeric(WET_RF3$RF)

#Get seasonal average
WET_RF5 <- as.vector(tapply(WET_RF4, gl(length(WET_RF4) / 6, 6), mean, na.rm = T))

# Dry Season
DRY_RF <- MRF100[MRF100a$Month != 1 & MRF100a$Month  != 02 &
    MRF100a$Month  != 3 & MRF100a$Month  != 04 &
    MRF100a$Month  != 11 & MRF100a$Month  != 12, ]
head(DRY_RF, 12)
DRY_RF2 <- as.numeric(DRY_RF$RF)

#Get seasonal average
DRY_RF3 <- as.vector(tapply(DRY_RF2, gl(length(DRY_RF2) / 6, 6), mean, na.rm = T))

########## Season Stats
##########   Wet
W_MRF_Max <- round(max(WET_RF5, na.rm = T), 0)
W_MRF_Min <- round(min(WET_RF5, na.rm = TRUE), 0)
W_MRF_MED <- round(median(WET_RF5, na.rm = TRUE), 0)
W_MRF_MEAN <- round(mean(WET_RF5, na.rm = TRUE), 0)
##########   Dry
D_MRF_Max <- round(max(DRY_RF3, na.rm = T), 0)
D_MRF_Min <- round(min(DRY_RF3, na.rm = TRUE), 0)
D_MRF_MED <- round(median(DRY_RF3, na.rm = TRUE), 0)
D_MRF_MEAN <- round(mean(DRY_RF3, na.rm = TRUE), 0)




short.date_Y = strftime(DRY_RF$Ndate, "%Y")

# Dry_RF = aggregate(as.numeric(DRY_RF$RF) ~ short.date_Y, FUN = mean)
# colnames(Dry_RF) <- c("Date","RF")

####### Seasonal Trends ########################
debug_print("80, Seasonal Trends")

#WET Season
YrRF.tsW <- ts(WET_RF5, c(1920), end = c(ey), frequency = 1)
myts1YW <- as.vector(window(YrRF.tsW, start = c(1920), end = c(ey)))
myts2YW <- as.vector(window(YrRF.tsW, start = c(1940), end = c(ey)))
myts3YW <- as.vector(window(YrRF.tsW, start = c(1960), end = c(ey)))
myts4YW <- as.vector(window(YrRF.tsW, start = c(1980), end = c(ey)))
myts5YW <- as.vector(window(YrRF.tsW, start = c(2000), end = c(ey)))
myts6YW <- as.vector(window(YrRF.tsW, start = c(2010), end = c(ey)))

#DRy Season
YrRF.tsD <- ts(DRY_RF3, c(1920), end = c(ey), frequency = 1)
myts1YD <- as.vector(window(YrRF.tsD, start = c(1920), end = c(ey)))
myts2YD <- as.vector(window(YrRF.tsD, start = c(1940), end = c(ey)))
myts3YD <- as.vector(window(YrRF.tsD, start = c(1960), end = c(ey)))
myts4YD <- as.vector(window(YrRF.tsD, start = c(1980), end = c(ey)))
myts5YD <- as.vector(window(YrRF.tsD, start = c(2000), end = c(ey)))
myts6YD <- as.vector(window(YrRF.tsD, start = c(2010), end = c(ey)))

# WET and DRy # Annual
YDateT1 <- seq(as.Date("1920-01-01"), as.Date(ed), by = "years")
YDateT2 <- seq(as.Date("1940-01-01"), as.Date(ed), by = "years")
YDateT3 <- seq(as.Date("1960-01-01"), as.Date(ed), by = "years")
YDateT4 <- seq(as.Date("1980-01-01"), as.Date(ed), by = "years")
YDateT5 <- seq(as.Date("2000-01-01"), as.Date(ed), by = "years")
YDateT6 <- seq(as.Date("2010-01-01"), as.Date(ed), by = "years")

#WET Regression
LM1YW <- lm(myts1YW ~ YDateT1)
LM2YW <- lm(myts2YW ~ YDateT2)
LM3YW <- lm(myts3YW ~ YDateT3)
LM4YW <- lm(myts4YW ~ YDateT4)
LM5YW <- lm(myts5YW ~ YDateT5)
LM6YW <- lm(myts6YW ~ YDateT6)

# Dry Regression
LM1YD <- lm(myts1YD ~ YDateT1)
LM2YD <- lm(myts2YD ~ YDateT2)
LM3YD <- lm(myts3YD ~ YDateT3)
LM4YD <- lm(myts4YD ~ YDateT4)
LM5YD <- lm(myts5YD ~ YDateT5)
LM6YD <- lm(myts6YD ~ YDateT6)

#WET SLOPE
T1YW <- round(coefficients(summary(LM1YW))[2, 1], 3)
T2YW <- round(coefficients(summary(LM2YW))[2, 1], 3)
T3YW <- round(coefficients(summary(LM3YW))[2, 1], 3)
T4YW <- round(coefficients(summary(LM4YW))[2, 1], 3)
T5YW <- round(coefficients(summary(LM5YW))[2, 1], 3)
T6YW <- round(coefficients(summary(LM6YW))[2, 1], 3)

#DRY SLOPE
T1YD <- round(coefficients(summary(LM1YD))[2, 1], 3)
T2YD <- round(coefficients(summary(LM2YD))[2, 1], 3)
T3YD <- round(coefficients(summary(LM3YD))[2, 1], 3)
T4YD <- round(coefficients(summary(LM4YD))[2, 1], 3)
T5YD <- round(coefficients(summary(LM5YD))[2, 1], 3)
T6YD <- round(coefficients(summary(LM6YD))[2, 1], 3)

#Wet Pvale
LM1PYW <- round(coefficients(summary(LM1YW))[2, 4], 2)
LM2PYW <- round(coefficients(summary(LM2YW))[2, 4], 2)
LM3PYW <- round(coefficients(summary(LM3YW))[2, 4], 2)
LM4PYW <- round(coefficients(summary(LM4YW))[2, 4], 2)
LM5PYW <- round(coefficients(summary(LM5YW))[2, 4], 2)
LM6PYW <- round(coefficients(summary(LM6YW))[2, 4], 2)

#Dry Pvale
LM1PYD <- round(coefficients(summary(LM1YD))[2, 4], 2)
LM2PYD <- round(coefficients(summary(LM2YD))[2, 4], 2)
LM3PYD <- round(coefficients(summary(LM3YD))[2, 4], 2)
LM4PYD <- round(coefficients(summary(LM4YD))[2, 4], 2)
LM5PYD <- round(coefficients(summary(LM5YD))[2, 4], 2)
LM6PYD <- round(coefficients(summary(LM6YD))[2, 4], 2)

#Wet Regression
LM1RYW <- round(summary(LM1YW)$r.squared, 2)
LM2RYW <- round(summary(LM2YW)$r.squared, 2)
LM3RYW <- round(summary(LM3YW)$r.squared, 2)
LM4RYW <- round(summary(LM4YW)$r.squared, 2)
LM5RYW <- round(summary(LM5YW)$r.squared, 2)
LM6RYW <- round(summary(LM6YW)$r.squared, 2)

#Dry Regression
LM1RYD <- round(summary(LM1YD)$r.squared, 2)
LM2RYD <- round(summary(LM2YD)$r.squared, 2)
LM3RYD <- round(summary(LM3YD)$r.squared, 2)
LM4RYD <- round(summary(LM4YD)$r.squared, 2)
LM5RYD <- round(summary(LM5YD)$r.squared, 2)
LM6RYD <- round(summary(LM6YD)$r.squared, 2)


dpi <- 300
png(
  paste0(PATH_WITH_PROJECT_NAME, "RF_Trend.png"),
  width = 6.3 * dpi,
  height = 7 * dpi,
  res = dpi
)

##########   Annual and Seasonal Plot
debug_print("81, Annual and Seasonal Plot")

par(mfrow = c(3, 1))
par(mar = c(4, 4, 4, 2))

MLIM1 <- max(c(myts1Y, myts1YW, myts1YD), na.rm = T)
YLIM <-  min(myts1Y)
MLIM <-  (MLIM1 + (MLIM1 * 0.45))

par(mai=c(0.3,0.6,0.2,0.2))
plot(myts1Y~YDateT1,ylab = paste0("Average Rainfall (",RFUnit2,"/month)"),type="l",col="blue",xlab="",xaxt="n",ylim=c(YLIM,MLIM),cex.axis =1.3,las=1)
title(paste("Rainfall Trend 1920-",ey,":",UNIT_Ns[u]), line = 0.5,cex.main = 1.5)

legend("topright", c(paste0("1920-",ey," R2 = ",LM1RY, " p = ",LM1PY),
  #paste0("1940-",ey," Trend =", T2Y,", R2 =",LM2RY, " p = ",LM2PY),
  #paste0("1960-",ey," R2 = ",LM3RY, " p = ",LM3PY),
  paste0("1980-",ey," R2 = ",LM4RY, " p = ",LM4PY),
  #paste0("2000-",ey," R2 = ",LM5RY, " p = ",LM5PY)),
  paste0("2010-",ey," R2 = ",LM6RY, " p = ",LM6PY)),
  #lty = 1, col = c("darkred","darkorange","darkgreen","darkblue","purple","darkcyan"), lwd = 3)
  lty = 1, col = c("grey70","grey30","grey1"), lwd = 3)

legend("topleft",c("Annual"),cex=1.5,text.font=2, bty = "n")

ablineclip(lm(myts1Y~YDateT1),x1=-19000,x2=20000,col="grey70",lwd=3)
ablineclip(lm(myts4Y~YDateT4),x1=3500,x2=20000,col="grey30",lwd=3)
ablineclip(lm(myts6Y~YDateT6),x1=14000,x2=20000,col="grey1",lwd=3)

####### Wet Season ###########
debug_print("82, Wet Season")

par(mai = c(0.3, 0.6, 0.2, 0.2))
YLIM <-  min(myts1YW, na.rm = T)
plot(myts1YW~YDateT1,ylab = paste0("Average Rainfall (",RFUnit2,"/month)"),type="l",col="blue",xlab="",xaxt="n",ylim=c(YLIM,MLIM),cex.axis =1.3,las=1)

legend("topright", c(paste0("1920-",ey," R2 = ",LM1RYW, " p = ",LM1PYW),
  #paste0("1940-",ey," Trend =", T2YW,", R2 =",LM2RYW, " p = ",LM2PYW),
  #paste0("1960-",ey," R2 = ",LM3RYW, " p = ",LM3PYW),
  paste0("1980-",ey," R2 = ",LM4RYW, " p = ",LM4PYW),
  #paste0("2000-",ey," R2 = ",LM5RYW, " p = ",LM5PYW)),
  paste0("2010-",ey," R2 = ",LM6RYW, " p = ",LM6PYW)),
  #lty = 1, col = c("darkred","darkorange","darkgreen","darkblue","purple","darkcyan"), lwd = 3,cex=1)
  lty = 1, col = c("grey70","grey30","grey1"), lwd = 3)

legend("topleft",c("Wet Season"),cex=1.5,text.font=2, bty = "n")

ablineclip(lm(myts1YW~YDateT1),x1=-19000,x2=20000,col="grey70",lwd=3)
ablineclip(lm(myts4YW~YDateT4),x1=3500,x2=20000,col="grey30",lwd=3)
ablineclip(lm(myts6YW~YDateT6),x1=14000,x2=20000,col="grey1",lwd=3)

########## Dry Season
debug_print("83, Dry Season")

par(mai = c(0.3, 0.6, 0.2, 0.2))
YLIM <-  min(myts1YD, na.rm = T)
plot(myts1YD~YDateT1,ylab = paste0("Average Rainfall (",RFUnit2,"/month)"),type="l",col="blue",xlab="",ylim=c(YLIM,MLIM),cex.axis =1.3,las=1)
# axis(1, labels = T)
# title(main = "Average Wet Season Rainfall Pu'u Wa'awa'a (1920-2021)", line = 1)
legend("topright", c(paste0("1920-",ey," R2 = ",LM1RYD, " p = ",LM1PYD),
  #paste0("1940-",ey," Trend =", T2YD,", R2 =",LM2RYD, " p = ",LM2PYD),
  #paste0("1960-",ey," R2 = ",LM3RYD, " p = ",LM3PYD),
  paste0("1980-",ey," R2 = ",LM4RYD, " p = ",LM4PYD),
  #paste0("2000-",ey," R2 = ",LM5RYD, " p = ",LM5PYD)),
  paste0("2010-",ey," R2 = ",LM6RYD, " p = ",LM6PYD)),
  #lty = 1, col = c("darkred","darkorange","darkgreen","darkblue","purple","darkcyan"), lwd = 3)
  lty = 1, col = c("grey70","grey30","grey1"), lwd = 3)

legend("topleft",c("Dry Season"),cex=1.5,text.font=2, bty = "n")

ablineclip(lm(myts1YD~YDateT1),x1=-20000,x2=20000,col="grey70",lwd=3)
ablineclip(lm(myts4YD~YDateT4),x1=3500,x2=20000,col="grey30",lwd=3)
ablineclip(lm(myts6YD~YDateT6),x1=14000,x2=20000,col="grey1",lwd=3)

dev.off()


### SPI/DROUGHT ###

#############################################################################
# Count All Drought Events Long-term 12 and Short term 12 and 3.
debug_print("84, SPI/DROUGHT")

Cell.SPICNT <- data.frame(matrix(ncol = 4, nrow = 4))
colnames(Cell.SPICNT) <- c("Total", "SPI_12_Long", "SPI_3_Short", "SPI_12_Short")
Cell.SPICNT[1:4, 1] <- c("Drought Events", "Moderate", "Severe", "Extreme")

### Ryan's code
debug_print("SPI FIG 12")
RF <- as.numeric(MRF100$RF)
RFDATA <- cbind(MRF100[, c(2, 3)])

##################################
#SPI-12 Calculation
SPI12 <- spi(RF, scale = 12, distribution = 'Gamma')
spi12 <- spi(ts(RF, freq = 12, start = c(1920, 1)), scale = 12)
#SPI-3 Calculation
SPI3 <- spi(RF, scale = 3, distribution = 'Gamma')

# extract the fitted values
spi_vals <- as.numeric(spi12$fitted)
time_vals <- time(spi12$fitted)

#PlOT ALL SPI
png(
  paste0(PATH_WITH_PROJECT_NAME, "SPI.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

# custom barplot style
plot(time_vals, spi_vals, type = "h",
     col = ifelse(spi_vals >= 0, "blue", "red"),
     main = paste0("SPI-12 1920-", ey, ": ", UNIT_Ns[u]),
     xlab = "Year", ylab = "SPI")
abline(h = 0, lwd = 2)

dev.off()

##########   Table Metrics SPI 12
debug_print("85, Table Metrics SPI 12")

Cell.DataSPI <- data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI) <- c("Start","End", "Duration","A Intensity", "P Intensity", "Magnitude")
Cell.DataSPI

##### Derek's Drought Code from Guam
debug_print("86, Derek's Drought Code from Guam")

# load rainfall dataset created around line 2828 above
wd <- paste0(PATH, "/")
setwd(wd)
csv <- paste0(UNIT_N[u], "_Monthly Rainfall_", RFUnit2, ".csv")

RF_Data <- read.csv(csv)
head(RF_Data)
tail(RF_Data)

# calculate SPI and SPEI
spi <- spi(ts(RF_Data$RF, freq = 12, start = c(1920, 1)), scale = 12)
plot(spi)

spi3 <- spi(ts(RF_Data$RF, freq = 12, start = c(1920, 1)), scale = 3)
plot(spi3)

# write spi to dataframe
spi_val <- spi$fitted
spi_val <- data.frame(spi = as.matrix(spi_val), date = time(spi_val))
spi_val$m.scale <- 12
spi_val$date <- sub("\\..*", "", spi_val$date)

head(spi_val, 20)
tail(spi_val)

spi3_val <- spi3$fitted
spi3_val <- data.frame(spi = as.matrix(spi3_val), date = time(spi3_val))
spi3_val$m.scale <- 3
spi3_val$date <- sub("\\..*", "", spi3_val$date)

head(spi3_val, 20)
tail(spi3_val)

# combine dataframes
SPI_ALL <- rbind(spi_val, spi3_val)
colnames(SPI_ALL) <- c("spi", "date", "m.scale")
head(SPI_ALL, 20)
tail(SPI_ALL, 20)

# determine if last year is complete, or how many months are missing
# complete will be 24 months in 1 year because two time scales
ly <- nrow(SPI_ALL[which(SPI_ALL$date == max(SPI_ALL$date)), ])
ly
m <- 24 - ly
m

# make extra rows so the last year is complete
if (m != 0) {
  extra <- data.frame(matrix(ncol = 3, nrow = m))
  colnames(extra) <- c("spi", "date", "m.scale")
  max(SPI_ALL$date)
  extra$date <- max(SPI_ALL$date)
  extra[1:(m / 2), ]$m.scale <- 3
  extra[((m / 2) + 1):nrow(extra), ]$m.scale <- 12
  # add extra rows to dataset
  SPI_ALL <- rbind(SPI_ALL, extra)
}

head(SPI_ALL)

# sort by year
SPI_ALL <- SPI_ALL[order(SPI_ALL$date, SPI_ALL$m.scale), ]

# fix date column
SPI_ALL$date2 <- NA
SPI_ALL$date2 <- rep(c(1:12))
SPI_ALL$date <- as.Date(paste0(SPI_ALL$date, "/", SPI_ALL$date2, "/01"), format = "%Y/%m/%d")
head(SPI_ALL, 20)
tail(SPI_ALL, 30)

SPI_ALL <- SPI_ALL[1:3]
colnames(SPI_ALL) <- c("SPI", "date", "m.scale")

# sort by date
SPI_ALL <- SPI_ALL[order(as.Date(SPI_ALL$date)), ]

# Convert drought events (negative values) to positive
SPI_ALL$spi_negs <- ifelse(SPI_ALL$SPI > 0, 0, SPI_ALL$SPI)
SPI_ALL$spi_negs <- abs(SPI_ALL$spi_negs)

summary(SPI_ALL$spi_negs)
head(SPI_ALL, 50)
tail(SPI_ALL, 50)

# Save SPI_ALL drought intensity (inverted SPI) dataset
write.csv(
  SPI_ALL,
  paste0(PATH_WITH_PROJECT_NAME, "SPI_NEGS_ALL.csv"),
  row.names = F
)

### Drought Event Count for SPI-12 ###

spi2 <- SPI_ALL[which(SPI_ALL$m.scale == 12), ]
head(spi2)

# create binary column of drought (SPI>0) yes or no
spi2$drought <- ifelse(spi2$spi_negs > 0, 1, 0)
head(spi2, 50)
tail(spi2, 50)

spi2$DRGHT_ct<-with(spi2, (drought == 1)*
    ave(drought, rleid(drought == 1),
      FUN=seq_along))

summary(spi2$DRGHT_ct)

# subset data for drought event months only
spi3 <- spi2[which(spi2$DRGHT_ct >= 1), ]
head(spi3, 20)
tail(spi3, 50)

# create empty event count column
spi3$event_ct <- 0

### fill in event_ct column
for (r in 1:nrow(spi3)) {
  
  SPI_M<-spi3[r,]
  SPI_M
  
  if(r == 1) {spi3[r,]$event_ct<-r}
  if(SPI_M$DRGHT_ct > 1) {spi3[r,]$event_ct<-spi3[r-1,]$event_ct}
  if(SPI_M$DRGHT_ct == 1 && r > 1) {spi3[r,]$event_ct<-spi3[r-1,]$event_ct + 1}
}

head(spi3, 100)
tail(spi3, 50)
summary(spi3$event_ct)

### keep only events with peak spi_negs >=1
event <- unique(spi3$event_ct)

for (x in event) {
  
  sub<-spi3[which(spi3$event_ct == x),]
  
  if(max(sub$spi_negs)<1) {spi3[which(spi3$event_ct == x),]$event_ct<-0}
}

head(spi3, 100)
tail(spi3, 50)

### make final event_ct column
# subset for events only
spi4 <- spi3[which(spi3$event_ct > 0), ]
spi4$event_ct2 <- NA
spi4

event <- data.frame(unique(spi4$event_ct))
event

# fill event_Ct2 column with event number
for (x in 1:nrow(event)) {
  
  # get event number
  e <- event[x, ]
  
  # for each event, add event_ct2 value
  spi4[which(spi4$event_ct == e), ]$event_ct2 <- x
}

tail(spi4, 10)

# merge event_ct back into full SPI dataset
head(spi2)
head(spi4)

spi5 <- merge(spi2, spi4, all = TRUE)
spi5 <- spi5[order(as.Date(spi5$date)), ]
head(spi5, 10)
tail(spi5, 50)
summary(spi5$event_ct2)

# get rid of wrong event column
spi5 <- subset(spi5, select = -c(event_ct))

colnames(spi5)[which(names(spi5) == "event_ct2")] <- "event_ct"
head(spi5, 50)

### create column which labels each drought event by intensity
# create new copy of dataframe
spi6 <- spi5

# get total events count as list
events <- as.list(1:max(spi6$event_ct, na.rm = T))

# loop through each event and label by maximum SPI value
for (x in events) {
  
  SPI_I <- spi5[which(spi5$event_ct == x), ]
  SPI_I$months <- NA
  SPI_I$intensity <- NA
  SPI_I$peak <- NA
  SPI_I$mean <- NA
  SPI_I$mag <- NA
  SPI_I
  
  SPI_I$months <- nrow(SPI_I)
  
  if(max(SPI_I$spi_negs) <= 1.5) {SPI_I$intensity<-1}
  if(max(SPI_I$spi_negs) > 1.5 && max(SPI_I$spi_negs) <= 2) {SPI_I$intensity<-2}
  if(max(SPI_I$spi_negs) > 2) {SPI_I$intensity<-3}
  
  SPI_I$peak <- max(SPI_I$spi_negs)
  SPI_I$mean <- mean(SPI_I$spi_negs)
  SPI_I$mag <- sum(SPI_I$spi_negs)
  
  # keep only these columns
  SPI_I <- SPI_I[c("date", "months", "intensity", "peak", "mean", "mag")]
  SPI_I
  
  spi6 <- merge(spi6, SPI_I, by = "date", all.x = TRUE)
  
  # merge intensity columns together
  if (x>1) {spi6$months<-ifelse(is.na(spi6$months.x), spi6$months.y, spi6$months.x)}
  if (x>1) {spi6$intensity<-ifelse(is.na(spi6$intensity.x), spi6$intensity.y, spi6$intensity.x)}
  if (x>1) {spi6$peak<-ifelse(is.na(spi6$peak.x), spi6$peak.y, spi6$peak.x)}
  if (x>1) {spi6$mean<-ifelse(is.na(spi6$mean.x), spi6$mean.y, spi6$mean.x)}
  if (x>1) {spi6$mag<-ifelse(is.na(spi6$mag.x), spi6$mag.y, spi6$mag.x)}
  
  if (x>1) {spi6<-subset(spi6, select=-c(months.x, months.y, intensity.x, 
    intensity.y, peak.x, peak.y, mean.x, mean.y, mag.x, mag.y))}
}

head(spi6)
head(spi6, 100)
tail(spi6, 50)
summary(spi6$intensity)
summary(spi6$peak)
summary(spi6$months)

### make table of start and end dates for each drought intensity event
debug_print("87, make table of start and end dates for each drought intensity event")

Cell.DataSPI <- data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI) <- c("Start","End", "Duration","A Intensity", "P Intensity", "Magnitude")
Cell.DataSPI$Start <- as.Date(Cell.DataSPI$Start)
Cell.DataSPI$End <- as.Date(Cell.DataSPI$End)
Cell.DataSPI
colnames(Cell.DataSPI)

### loop through rows and fill table

head(spi6, 50)
tail(spi6, 50)

for (y in 1:nrow(spi6[which(!is.na(spi6$SPI)), ])) {
  
  # get row
  row <- spi6[y, ]
  
  # get previous and next row
  prev <- spi6[y - 1, ]
  nex <- spi6[y + 1, ]
  
  # get event number
  number <- row$event_ct
  
  if(!is.na(row$intensity) && is.na(prev$intensity)) {Cell.DataSPI[number,]$Start<-as.Date(nex$date)}
  #if(!is.na(row$intensity) && is.na(prev$intensity)) {Cell.DataSPI[number,]$intensity<-row$intensity}
  if(!is.na(row$intensity) && is.na(nex$intensity)) {Cell.DataSPI[number,]$End<-as.Date(nex$date)}
  if(y==nrow(spi6[which(!is.na(spi6$SPI)),]) && !is.na(row$intensity)) {Cell.DataSPI[number,]$End<-NA}
  if(!is.na(row$intensity)) {Cell.DataSPI[number,]$Duration<-row$months}
  #if(!is.na(row$intensity)) {Cell.DataSPI[number,]$event_ct<-row$event_ct}
  if(!is.na(row$intensity)) {Cell.DataSPI[number,]$`P Intensity`<-row$peak}
  if(!is.na(row$intensity)) {Cell.DataSPI[number,]$`A Intensity`<-row$mean}
  if(!is.na(row$intensity)) {Cell.DataSPI[number,]$Magnitude<-row$mag}
}

# fix date formats
Cell.DataSPI$Start2 <- as.yearmon(sub(" .*", "", Cell.DataSPI$Start), "%Y-%m-%d")
Cell.DataSPI$End2 <- as.yearmon(sub(" .*", "", Cell.DataSPI$End), "%Y-%m-%d")

# order by date
Cell.DataSPI <- Cell.DataSPI[order(Cell.DataSPI$Start), ]

# df <- df[order(df$date), ]
Cell.DataSPI

###### Back to Ryan's code #####

write.csv(
  Cell.DataSPI,
  paste0(PATH_WITH_PROJECT_NAME, "Drought History.csv"),
  row.names = F
)

##########   Count Droughts For Figure - Derek's edits
debug_print("88, Count Droughts For Figure")

Cell.DataSPI$`P Intensity` <- as.numeric(Cell.DataSPI$`P Intensity`)

EX_Cnt <- sum(Cell.DataSPI[5] > 2)
SV_Cnt <- sum(Cell.DataSPI[5] > 1.5)
MO_Cnt <- sum(Cell.DataSPI[5] > 1)

D_Cnt <- MO_Cnt
SV_Cnt2 <- SV_Cnt - EX_Cnt
MO_Cnt2 <- D_Cnt -  SV_Cnt2 - EX_Cnt

# create SPI count dataframe
Cell.SPICNT <- data.frame(matrix(ncol = 4, nrow = 4))
colnames(Cell.SPICNT) <- c("Total", "SPI_12_Long", "SPI_3_Short", "SPI_12_Short")
Cell.SPICNT[1:4, 1] <- c("Drought Events", "Moderate", "Severe", "Extreme")

Cell.SPICNT[1:4, 2] <- c(D_Cnt, MO_Cnt2, SV_Cnt2, EX_Cnt)
Cell.SPICNT

##########   Remove Positive SPI and make absolute values
debug_print("89, Remove Postive SPI and make absloute values")

### Derek's code
SPIVEC <- SPI_ALL[which(SPI_ALL$m.scale == 12), ]$SPI
# SPIVEC[SPIVEC > 0] <- 0
# SPIVEC_Abs <- as.vector(abs(SPIVEC))

M_SPI.ts <- ts(SPIVEC, c(1920, 1), end = c(ey, em), frequency = 12)
myts66 <- as.vector(window(M_SPI.ts, start = c(1920, 1), end = c(ey, em)))
DateT1 <- as.Date(seq(as.Date("1920-01-01"), as.Date(ed), by = "months"))
#short.date_M = strftime(MonthlyRF$Ndate, "%Y-%m")

xx <- data.frame(DateT1, myts66)
colnames(xx) <- c("DT", "SP")
xx$DT <- as.Date(xx$DT)

write.csv(xx,
  paste0(PATH_WITH_PROJECT_NAME, "SPI_12.csv"),
  row.names = F)

# Calculate monthly climatology SPI-12 values (min, mean, max)
spi12a <- xx
head(spi12a, 20)

min <- min(spi12a$SP, na.rm = T)
mean <- mean(spi12a$SP, na.rm = T)
max <- max(spi12a$SP, na.rm = T)

### Make monthly min, mean, max SPI value dataset

# get month and year columns
spi12a$year <- substr(spi12a$DT, 1, 4)
spi12a$month <- as.numeric(substr(spi12a$DT, 6, 7))

# remove NA rows
spi12a <- spi12a[which(!is.na(spi12a$SP)), ]

# aggregate over months
spi12amean <- sapply(split(spi12a$SP, spi12a$month), mean)
spi12amin <- sapply(split(spi12a$SP, spi12a$month), min)
spi12amax <- sapply(split(spi12a$SP, spi12a$month), max)
spi12amean

spi12am <- rbind(spi12amean, spi12amin, spi12amax)
rownames(spi12am) <- c("mean", "min", "max")
spi12am

write.csv(spi12am, paste0(PATH_WITH_PROJECT_NAME, "SPI_12 monthly.csv"))

# change positive SPI values to 0 and make drought figure
SPIVEC <- SPI_ALL[which(SPI_ALL$m.scale == 12), ]$SPI
SPIVEC[SPIVEC > 0] <- 0
SPIVEC_Abs <- as.vector(abs(SPIVEC))

M_SPI.ts <- ts(SPIVEC_Abs, c(1920, 1), end = c(ey, em), frequency = 12)
myts66 <- as.vector(window(M_SPI.ts, start = c(1920, 1), end = c(ey, em)))
DateT1 <- as.Date(seq(as.Date("1920-01-01"), as.Date(ed), by = "months"))
#short.date_M = strftime(MonthlyRF$Ndate, "%Y-%m")

xx <- data.frame(DateT1, myts66)
colnames(xx) <- c("DT", "SP")
xx$DT <- as.Date(xx$DT)
xx
tail(xx)

dpi = 300
png(
  paste0(PATH_WITH_PROJECT_NAME, "Drought_History.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

print

ggplot(xx, aes(x = DT, y = SP)) +
  geom_area(fill="darkorange", color="black") +
  xlab("") +
  
  labs(
    title = paste0("SPI-12 Drought Events 1920 -",ey,": ",UNIT_Ns[u]),
    x = "",
    y = "Drought Intensity") +
  geom_hline(yintercept=2, linetype="dashed", color = "darkred", size = 1) + 
  geom_hline(yintercept=1.5, linetype="dashed", color = "red", size = 1) +
  geom_hline(yintercept=1, linetype="dashed", color = "orange", size = 1)  +
  
  annotate("text", x = xx$DT[350], y = 3.4, label = paste0("Drought Events = ",  D_Cnt),size=5, fontface="bold") +
  annotate("text", x = xx$DT[350], y = 3.1, label = paste0("Moderate Droughts = ", MO_Cnt2),size=5, fontface="bold",colour = "orange") +
  annotate("text", x = xx$DT[350], y = 2.8, label = paste0("Severe Droughts = ", SV_Cnt2),size=5, fontface="bold",colour = "red") +
  annotate("text", x = xx$DT[350], y = 2.5, label = paste0("Extreme Droughts = ", EX_Cnt),size=5, fontface="bold",colour = "darkred") +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )



dev.off()





##########################################################################################################################################
debug_print("90, Short (_S) timescales")

#Short (_S) timescales - Derek's version

#### SPI-3

head(SPI_ALL)
spi2 <- SPI_ALL[which(SPI_ALL$m.scale == 3), ]
spi2 <- spi2[which(spi2$date >= as.Date("1990-01-01")), ]
head(spi2)

# create binary column of drought (SPI>=1) yes or no
spi2$drought <- ifelse(spi2$spi_negs > 0, 1, 0)
head(spi2, 30)

# make column of consecutive drought months
spi2$DRGHT_ct<-with(spi2, (drought == 1)*
    ave(drought, rleid(drought == 1),
      FUN=seq_along))

# subset data for drought event months only
spi3 <- spi2[which(spi2$drought >= 1), ]
head(spi3, 20)

# create empty event count column
spi3$event_ct <- 0

### fill in event_ct column
for (r in 1:nrow(spi3)) {
  
  SPI_M<-spi3[r,]
  SPI_M
  
  if(r == 1) {spi3[r,]$event_ct<-r} 
  if(SPI_M$DRGHT_ct > 1) {spi3[r,]$event_ct<-spi3[r-1,]$event_ct}
  if(SPI_M$DRGHT_ct == 1 && r > 1) {spi3[r,]$event_ct<-spi3[r-1,]$event_ct + 1}
}

head(spi3, 50)
summary(spi3$event_ct)

### keep only events with peak spi_negs >=1
event <- unique(spi3$event_ct)

for (x in event) {
  
  sub <- spi3[which(spi3$event_ct == x), ]
  
  if (max(sub$spi_negs) < 1) { spi3[which(spi3$event_ct == x), ]$event_ct <- 0 }
}

head(spi3, 100)

### make final event_ct column
# subset for events only
spi4 <- spi3[which(spi3$event_ct > 0), ]
spi4$event_ct2 <- NA
spi4

event <- data.frame(unique(spi4$event_ct))
event

# fill event_Ct2 column with event number
for (x in 1:nrow(event)) {
  
  # get event number
  e <- event[x, ]
  
  # for each event, add event_ct2 value
  spi4[which(spi4$event_ct == e), ]$event_ct2 <- x
}

spi4

# merge event_ct back into full SPI dataset
spi5 <- merge(spi2, spi4, all = TRUE)

# sort by date
spi5 <- spi5[order(as.Date(spi5$date)), ]
head(spi5, 10)
summary(spi5$event_ct2)

# get rid of wrong event column
spi5 <- subset(spi5, select = -c(event_ct))
colnames(spi5)[which(names(spi5) == "event_ct2")] <- "event_ct"

### create column which labels each drought event by intensity
# create new copy of dataframe
spi6 <- spi5

# get total events count as list
events <- as.list(1:max(spi6$event_ct, na.rm = T))

# loop through each event and label by maximum SPI value
for (x in events) {
  
  SPI_I <- spi5[which(spi5$event_ct == x), ]
  SPI_I$months <- NA
  SPI_I$intensity <- NA
  SPI_I$peak <- NA
  SPI_I$mean <- NA
  SPI_I$mag <- NA
  SPI_I
  
  SPI_I$months <- nrow(SPI_I)
  
  if(max(SPI_I$spi_negs) <= 1.5) {SPI_I$intensity<-1}
  if(max(SPI_I$spi_negs) > 1.5 && max(SPI_I$spi_negs) <= 2) {SPI_I$intensity<-2}
  if(max(SPI_I$spi_negs) > 2) {SPI_I$intensity<-3}
  
  SPI_I$peak <- max(SPI_I$spi_negs)
  SPI_I$mean <- mean(SPI_I$spi_negs)
  SPI_I$mag <- sum(SPI_I$spi_negs)
  
  # keep only these columns
  SPI_I <- SPI_I[c("date", "months", "intensity", "peak", "mean", "mag")]
  SPI_I
  
  spi6 <- merge(spi6, SPI_I, by = "date", all.x = TRUE)
  #
  # # sort by date
  # spi5<-spi5[order(as.Date(spi5$date)),]
  # head(spi5, 100)
  
  # merge intensity columns together
  if (x>1) {spi6$months<-ifelse(is.na(spi6$months.x), spi6$months.y, spi6$months.x)}
  if (x>1) {spi6$intensity<-ifelse(is.na(spi6$intensity.x), spi6$intensity.y, spi6$intensity.x)}
  if (x>1) {spi6$peak<-ifelse(is.na(spi6$peak.x), spi6$peak.y, spi6$peak.x)}
  if (x>1) {spi6$mean<-ifelse(is.na(spi6$mean.x), spi6$mean.y, spi6$mean.x)}
  if (x>1) {spi6$mag<-ifelse(is.na(spi6$mag.x), spi6$mag.y, spi6$mag.x)}
  
  if (x>1) {spi6<-subset(spi6, select=-c(months.x, months.y, intensity.x, 
    intensity.y, peak.x, peak.y, mean.x, mean.y, mag.x, mag.y))}
}

head(spi6, 100)
tail(spi6, 50)
summary(spi6$intensity)
summary(spi6$peak)
summary(spi6$months)

### make table of start and end dates for each drought intensity event
debug_print("91, make table of start and end dates for each drought intensity event")

Cell.DataSPI3_S <- data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI3_S) <- c("Start","End", "Duration","A Intensity", "P Intensity", "Magnitude")

Cell.DataSPI3_S$Start <- as.POSIXct(Cell.DataSPI3_S$Start)
Cell.DataSPI3_S$End <- as.POSIXct(Cell.DataSPI3_S$End)
Cell.DataSPI3_S
colnames(Cell.DataSPI3_S)

### loop through rows and fill table

for (y in 1:nrow(spi6[which(!is.na(spi6$SPI)), ])) {

  # get row
  row <- spi6[y, ]
  
  # get previous and next row
  prev <- spi6[y - 1, ]
  
  nex <- spi6[y + 1, ]
  
  # get event number
  number <- row$event_ct
  
  if(!is.na(row$intensity) && is.na(prev$intensity)) {Cell.DataSPI3_S[number,]$Start<-as.Date(nex$date)}
  # if(!is.na(row$intensity) && is.na(prev$intensity)) {event.t[number,]$intensity<-row$intensity}
  if(!is.na(row$intensity) && is.na(nex$intensity)) {Cell.DataSPI3_S[number,]$End<-as.Date(nex$date)}
  if(y==nrow(spi6[which(!is.na(spi6$SPI)),]) && !is.na(row$intensity)) {Cell.DataSPI3_S[number,]$End<-NA}
  if(!is.na(row$intensity)) {Cell.DataSPI3_S[number,]$Duration<-row$months}
  # if(!is.na(row$intensity)) {event.t[number,]$event_ct<-row$event_ct}
  if(!is.na(row$intensity)) {Cell.DataSPI3_S[number,]$`P Intensity`<-row$peak}
  if(!is.na(row$intensity)) {Cell.DataSPI3_S[number,]$`A Intensity`<-row$mean}
  if(!is.na(row$intensity)) {Cell.DataSPI3_S[number,]$Magnitude<-row$mag}
}

# fix date formats
Cell.DataSPI3_S$Start <- as.yearmon(sub(" .*", "", Cell.DataSPI3_S$Start), "%Y-%m-%d")
Cell.DataSPI3_S$End <- as.yearmon(sub(" .*", "", Cell.DataSPI3_S$End), "%Y-%m-%d")

Cell.DataSPI3_S

write.csv(
  Cell.DataSPI3_S,
  paste0(PATH_WITH_PROJECT_NAME, "Drought History SPI_3.csv"),
  row.names = F
)

##########   Count Droughts For Figure - Derek's edits
debug_print("92, Count Droughts For Figure")

Cell.DataSPI3_S$`P Intensity` <- as.numeric(Cell.DataSPI3_S$`P Intensity`)

EX_Cnt <- sum(Cell.DataSPI3_S[5] > 2)
SV_Cnt <- sum(Cell.DataSPI3_S[5] > 1.5)
MO_Cnt <- sum(Cell.DataSPI3_S[5] > 1)

D_Cnt <- MO_Cnt
SV_Cnt2 <- SV_Cnt - EX_Cnt
MO_Cnt2 <- D_Cnt -  SV_Cnt2 - EX_Cnt

Cell.SPICNT[1:4, 3] <- c(D_Cnt, MO_Cnt2, SV_Cnt2, EX_Cnt)
Cell.SPICNT

##########   Remove Positive SPI and make absolute values
debug_print("93, Remove Positive SPI and make absolute values")

### Derek's Version ###

# use full SPI (3 and 12) dataframe
head(SPI_ALL)

# subset for SPI3 and 1990 - 2020 time period
SPIVEC <- SPI_ALL[which(SPI_ALL$m.scale == 3), ]
SPIVEC <- SPIVEC[which(SPIVEC$date >= as.Date("1990-01-01")), ]$SPI

# change positive SPI values to 0
SPIVEC[SPIVEC > 0] <- 0

# convert negatives to absolute
SPIVEC_Abs <- as.vector(abs(SPIVEC))

# make dataframe of date and drought (absolute SPI) value
M_SPI.ts <- ts(SPIVEC_Abs, c(1990, 1), end = c(ey, em), frequency = 12)
myts66 <- as.vector(window(M_SPI.ts, start = c(1990, 1), end = c(ey, em)))
DateT1 <- as.Date(seq(as.Date("1990-01-01"), as.Date(ed), by = "months"))
xx <- data.frame(DateT1, myts66)
colnames(xx) <- c("DT", "SP")
xx$DT <- as.Date(xx$DT)
xx

write.csv(xx,
  paste0(PATH_WITH_PROJECT_NAME, "SPI_3.csv"),
  row.names = F)

# plot and write figure
dpi = 300
png(
  paste0(PATH_WITH_PROJECT_NAME, "Drought_HistoryS_3.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

print(
  ggplot(xx, aes(x = DT, y = SP)) +
    geom_area(fill = "darkorange", color = "black") +
    xlab("") +
    
    labs(title = paste0("SPI-3 Drought Events 1990-",ey,": ",UNIT_Ns[u]),
      x = "",
      y = "Drought Intensity") +
    geom_hline(yintercept=2, linetype="dashed", color = "darkred", size = 1) + 
    geom_hline(yintercept=1.5, linetype="dashed", color = "red", size = 1) +
    geom_hline(yintercept=1, linetype="dashed", color = "orange", size = 1)  +
    annotate("text", x = xx$DT[70], y = 3.4, label = paste0("Drought Events = ",  D_Cnt),size=5, fontface="bold") +
    annotate("text", x = xx$DT[70], y = 3.1, label = paste0("Moderate Droughts = ", MO_Cnt2),size=5, fontface="bold",colour = "orange") +
    annotate("text", x = xx$DT[70], y = 2.8, label = paste0("Severe Droughts = ", SV_Cnt2),size=5, fontface="bold",colour = "red") +
    annotate("text", x = xx$DT[70], y = 2.5, label = paste0("Extreme Droughts = ", EX_Cnt),size=5, fontface="bold",colour = "darkred")+
    annotate("text", x = xx$DT[300], y = 3.1, label = "SPI-3",size=12, fontface="bold",colour = "black") +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16, face = "bold"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  
)

dev.off()


######################################################################################################################################



debug_print("94, SPI-12 short")

#### SPI-12 short

head(SPI_ALL)
spi2 <- SPI_ALL[which(SPI_ALL$m.scale == 12), ]
spi2 <- spi2[which(spi2$date >= as.Date("1990-01-01")), ]
head(spi2)

# create binary column of drought (SPI>=1) yes or no
spi2$drought <- ifelse(spi2$spi_negs > 0, 1, 0)
head(spi2, 30)

# make column of consecutive drought months
spi2$DRGHT_ct<-with(spi2, (drought == 1)*
    ave(drought, rleid(drought == 1),
      FUN=seq_along))

# subset data for drought event months only
spi3 <- spi2[which(spi2$DRGHT_ct >= 1), ]
head(spi3, 20)

# create empty event count column
spi3$event_ct <- 0

### fill in event_ct column
for (r in 1:nrow(spi3)) {
  
  SPI_M<-spi3[r,]
  SPI_M
  
  if(r == 1) {spi3[r,]$event_ct<-r}
  if(SPI_M$DRGHT_ct > 1) {spi3[r,]$event_ct<-spi3[r-1,]$event_ct}
  if(SPI_M$DRGHT_ct == 1 && r > 1) {spi3[r,]$event_ct<-spi3[r-1,]$event_ct + 1}
}

head(spi3, 100)
summary(spi3$event_ct)

### keep only events with peak spi_negs >=1
event <- unique(spi3$event_ct)

for (x in event) {
  
  sub <- spi3[which(spi3$event_ct == x), ]
  
  if (max(sub$spi_negs) < 1) { spi3[which(spi3$event_ct == x), ]$event_ct <- 0 }
}

head(spi3, 100)

### make final event_ct column
# subset for events only
spi4 <- spi3[which(spi3$event_ct > 0), ]
spi4$event_ct2 <- NA
spi4

event <- data.frame(unique(spi4$event_ct))
event

# fill event_Ct2 column with event number
for (x in 1:nrow(event)) {
  
  # get event number
  e <- event[x, ]
  
  # for each event, add event_ct2 value
  spi4[which(spi4$event_ct == e), ]$event_ct2 <- x
}

spi4

# merge event_ct back into full SPI dataset
spi5 <- merge(spi2, spi4, all = TRUE)

# sort by date
spi5 <- spi5[order(as.Date(spi5$date)), ]
head(spi5, 10)
summary(spi5$event_ct2)

# get rid of wrong event column
spi5 <- subset(spi5, select = -c(event_ct))
colnames(spi5)[which(names(spi5) == "event_ct2")] <- "event_ct"

### create column which labels each drought event by intensity
# create new copy of dataframe
spi6 <- spi5

# get total events count as list
events <- as.list(1:max(spi6$event_ct, na.rm = T))

# loop through each event and label by maximum SPI value
for (x in events) {
  
  SPI_I<-spi5[which(spi5$event_ct == x),]
  SPI_I$months<-NA
  SPI_I$intensity<-NA
  SPI_I$peak<-NA
  SPI_I$mean<-NA
  SPI_I$mag<-NA
  SPI_I
  
  SPI_I$months<-nrow(SPI_I)
  
  if(max(SPI_I$spi_negs) <= 1.5) {SPI_I$intensity<-1}
  if(max(SPI_I$spi_negs) > 1.5 && max(SPI_I$spi_negs) <= 2) {SPI_I$intensity<-2}
  if(max(SPI_I$spi_negs) > 2) {SPI_I$intensity<-3}
  
  SPI_I$peak<-max(SPI_I$spi_negs)
  SPI_I$mean<-mean(SPI_I$spi_negs)
  SPI_I$mag<-sum(SPI_I$spi_negs)
  
  # keep only these columns
  SPI_I<-SPI_I[c("date","months","intensity","peak","mean","mag")]
  SPI_I
  
  spi6<-merge(spi6, SPI_I, by="date", all.x=T)
  
  # merge intensity columns together
  if (x>1) {spi6$months<-ifelse(is.na(spi6$months.x), spi6$months.y, spi6$months.x)}
  if (x>1) {spi6$intensity<-ifelse(is.na(spi6$intensity.x), spi6$intensity.y, spi6$intensity.x)}
  if (x>1) {spi6$peak<-ifelse(is.na(spi6$peak.x), spi6$peak.y, spi6$peak.x)}
  if (x>1) {spi6$mean<-ifelse(is.na(spi6$mean.x), spi6$mean.y, spi6$mean.x)}
  if (x>1) {spi6$mag<-ifelse(is.na(spi6$mag.x), spi6$mag.y, spi6$mag.x)}
  
  if (x>1) {spi6<-subset(spi6, select=-c(months.x, months.y, intensity.x, 
    intensity.y, peak.x, peak.y, mean.x, mean.y, mag.x, mag.y))}
}

head(spi6, 100)
summary(spi6$intensity)
summary(spi6$peak)
summary(spi6$months)

### make table of start and end dates for each drought intensity event
debug_print("95, make table of start and end dates for each drought intensity event")

Cell.DataSPI12_S <- data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI12_S) <- c("Start","End", "Duration","A Intensity", "P Intensity", "Magnitude")

Cell.DataSPI12_S$Start <- as.POSIXct(Cell.DataSPI12_S$Start)
Cell.DataSPI12_S$End <- as.POSIXct(Cell.DataSPI12_S$End)
Cell.DataSPI12_S
colnames(Cell.DataSPI12_S)

### loop through rows and fill table

for (y in 1:nrow(spi6[which(!is.na(spi6$SPI)),])) {
  
  # get row
  row<-spi6[y,]
  row
  
  # get previous and next row
  prev<-spi6[y-1,]
  
  nex<-spi6[y+1,]
  
  # get event number
  number<-row$event_ct
  
  if(!is.na(row$intensity) && is.na(prev$intensity)) {Cell.DataSPI12_S[number,]$Start<-as.Date(nex$date)}
  # if(!is.na(row$intensity) && is.na(prev$intensity)) {event.t[number,]$intensity<-row$intensity}
  if(!is.na(row$intensity) && is.na(nex$intensity)) {Cell.DataSPI12_S[number,]$End<-as.Date(nex$date)}
  if(y==nrow(spi6[which(!is.na(spi6$SPI)),]) && !is.na(row$intensity)) {Cell.DataSPI12_S[number,]$End<-NA}
  if(!is.na(row$intensity)) {Cell.DataSPI12_S[number,]$Duration<-row$months}
  # if(!is.na(row$intensity)) {event.t[number,]$event_ct<-row$event_ct}
  if(!is.na(row$intensity)) {Cell.DataSPI12_S[number,]$`P Intensity`<-row$peak}
  if(!is.na(row$intensity)) {Cell.DataSPI12_S[number,]$`A Intensity`<-row$mean}
  if(!is.na(row$intensity)) {Cell.DataSPI12_S[number,]$Magnitude<-row$mag}
}

# fix date formats
Cell.DataSPI12_S$Start <- as.yearmon(sub(" .*", "", Cell.DataSPI12_S$Start), "%Y-%m-%d")
Cell.DataSPI12_S$End <- as.yearmon(sub(" .*", "", Cell.DataSPI12_S$End), "%Y-%m-%d")

Cell.DataSPI12_S
Cell.DataSPI12_S[which(Cell.DataSPI12_S$`P Intensity` > 1.5), ]


##########   Count Droughts For Figure - Derek's edits
debug_print("96, Count Droughts For Figure")

Cell.DataSPI12_S$`P Intensity` <- as.numeric(Cell.DataSPI12_S$`P Intensity`)

EX_Cnt <- sum(Cell.DataSPI12_S[5] > 2)
SV_Cnt <- sum(Cell.DataSPI12_S[5] > 1.5)
MO_Cnt <- sum(Cell.DataSPI12_S[5] > 1)

D_Cnt <- MO_Cnt
SV_Cnt2 <- SV_Cnt - EX_Cnt
MO_Cnt2 <- D_Cnt -  SV_Cnt2 - EX_Cnt

Cell.SPICNT[1:4, 4] <- c(D_Cnt, MO_Cnt2, SV_Cnt2, EX_Cnt)
Cell.SPICNT

write.csv(
  Cell.SPICNT,
  paste0(PATH_WITH_PROJECT_NAME, "Drought Count.csv"),
  row.names = F
)

##########   Remove Positive SPI and make absolute values
debug_print("97, Remove Positive SPI and make absolute values")

### Derek's Version ###

# use full SPI (3 and 12) dataframe
head(SPI_ALL)

# subset for SPI12 and 1990 - 2020 time period
SPIVEC <- SPI_ALL[which(SPI_ALL$m.scale == 12), ]
SPIVEC <- SPIVEC[which(SPIVEC$date >= as.Date("1990-01-01")), ]$SPI

# change positive SPI values to 0
SPIVEC[SPIVEC > 0] <- 0

# convert negatives to absolute
SPIVEC_Abs <- as.vector(abs(SPIVEC))

# make dataframe of date and drought (absolute SPI) value
M_SPI.ts <- ts(SPIVEC_Abs, c(1990, 1), end = c(ey, em), frequency = 12)
myts66 <- as.vector(window(M_SPI.ts, start = c(1990, 1), end = c(ey, em)))
DateT1 <- as.Date(seq(as.Date("1990-01-01"), as.Date(ed), by = "months"))
xx <- data.frame(DateT1, myts66)
colnames(xx) <- c("DT", "SP")
xx$DT <- as.Date(xx$DT)
xx

# plot and write figure
dpi = 300
png(
  paste0(PATH_WITH_PROJECT_NAME, "Drought_HistoryS_12.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

print(
  ggplot(xx, aes(x = DT, y = SP)) +
    geom_area(fill="darkorange", color="black") +
    xlab("") +
    
    labs(
      title = paste0("SPI-12 Drought Events 1990-",ey,": ",UNIT_Ns[u]),
      x = "",
      y = "Drought Intensity") +
    geom_hline(yintercept=2, linetype="dashed", color = "darkred", size = 1) + 
    geom_hline(yintercept=1.5, linetype="dashed", color = "red", size = 1) +
    geom_hline(yintercept=1, linetype="dashed", color = "orange", size = 1) +
    annotate("text", x = xx$DT[70], y = 3.4, label = paste0("Drought Events = ",  D_Cnt),size=5, fontface="bold") +
    annotate("text", x = xx$DT[70], y = 3.1, label = paste0("Moderate Droughts = ", MO_Cnt2),size=5, fontface="bold",colour = "orange") +
    annotate("text", x = xx$DT[70], y = 2.8, label = paste0("Severe Droughts = ", SV_Cnt2),size=5, fontface="bold",colour = "red") +
    annotate("text", x = xx$DT[70], y = 2.5, label = paste0("Extreme Droughts = ", EX_Cnt),size=5, fontface="bold",colour = "darkred")+
    annotate("text", x = xx$DT[300], y = 3.1, label = "SPI-12",size=12, fontface="bold",colour = "black") +
      theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )
)

dev.off()


################################################################################################################################################

##########  MEI

# load ENSO phase dataset (using same ONI dataset from Guam analysis)
enso<-ONI_raw
head(enso)

# format Date column
enso$Date <- as.Date(enso$Date)

### plot ENSO phases over time

# ---- Define color bands ----
data_breaks <- data.frame(
  start = c(-2.8, -1.5, -0.5, 0.5, 1.5),
  end   = c(-1.5, -0.5, 0.5, 1.5, 2.8),
  colors = c("blue3", "lightskyblue1", "white", "indianred1", "red3"),
  label = c("Strong La Niña", "Weak La Niña", "Neutral", "Weak El Niño", "Strong El Niño")
)

# ---- Define positions ----
y_positions <- (data_breaks$start + data_breaks$end) / 2
x_min <- min(enso$Date, na.rm = TRUE)
x_max <- max(enso$Date, na.rm = TRUE)

# ---- Base plot ----
p <- ggplot() +
  # ENSO phase color bands
  geom_rect(
    data = data_breaks,
    aes(ymin = start, ymax = end,
        xmin = x_min, xmax = x_max),
    fill = data_breaks$colors, alpha = 0.5
  ) +
  # ONI line
  geom_line(
    data = enso,
    aes(x = Date, y = ONI, group = 1),
    linewidth = 1
  ) +
  # Black horizontal line at 0
  geom_hline(yintercept = 0, color = "black") +
  # Axes
  scale_x_date(
    limits = c(x_min, x_max),
    date_breaks = "5 years",
    labels = date_format("%Y"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title = "ENSO Phases based on Sea Surface Temperature (SST)",
    y = "Change in SST (°F)",
    x = "Year"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.65, size = 16),  # increased from 12 → 16
    axis.text.y = element_text(size = 16),                            # increased from 12 → 16
    axis.title.x = element_text(size = 18),                           # optional: larger axis titles
    axis.title.y = element_text(size = 18),
    plot.title  = element_text(size = 24, face = "bold"),             # optional: larger title
    plot.margin = margin(10, 160, 10, 10)
  ) +
  coord_cartesian(clip = "off")

# ---- Add labels outside the plot ----
for (i in seq_along(data_breaks$label)) {
  p <- p + annotation_custom(
    grid::textGrob(
      label = data_breaks$label[i],
      gp = gpar(fontsize = 18, fontface = "bold"),
      just = "left"
    ),
    xmin = x_max + 100,  # nudged ~50 days to right of plot
    xmax = x_max + 100,
    ymin = y_positions[i],
    ymax = y_positions[i]
  )
}

p

# ---- Save plot as PNG ----

# Define the output dimensions
width_px <- 3902
height_px <- 1929

png(
  filename = paste0(PATH_WITH_PROJECT_NAME, "ENSO_timeseries.png"),
  width = width_px,
  height = height_px,
  res = dpi,
  units = "px",
  bg = "white"
)

# Print or draw your ggplot here
print(p)

# Close device
dev.off()

##########   UNIT
debug_print("98, MEI")

Cell.MEI <- data.frame(matrix(nrow = 5, ncol = 9))
colnames(Cell.MEI) <- c("Phase","W-Mean","D-Mean","W-Max","D-Max","W-Min","D-Min","W-Count","D-Count")
Cell.MEI[1:5, 1] <- c("Strong EL", "Weak EL", "Neutral", "Weak LA", "Strong LA")

##########   All MEI
head(MEI)
tail(MEI)

# get last year
ly <- 1949 + nrow(MEI)
# ly<-max(MRF100$Year)
# MEI$Year <- seq(1950, ly, 1)

#Cell.MEI_All[u:1] <- UNIT_N[u]
#Cell.MEI_All[u,8] <- W_MRF_MEAN
#Cell.MEI_All[UNIT_C+u,8] <- D_MRF_MEAN
# MEI_W <- MEI_W[2:nrow(MEI_W), ]
MEI_D <- subset(MEI, select = -c(MEI_W))
# remove first row from wet season (1950 wet season not available because it starts in 1949)
MEI_W <- subset(MEI, select = -c(MEI_D))

debug_print("98.4, MEI")

# ## Can read in the csv from above to start script here
# MRF100<-read.csv(paste0(OUTPUTS_FOLDER,UNIT_N[u],"/",UNIT_N[u],"_Monthly Rainfall_in.csv"))

head(MRF100)
tail(MRF100)

# # If starting this section from monthly rainfall csv and only want up to 2012
# Cell.AF_Maps<-MRF100[which(MRF100$Year<2013),]

# # Using MEI and Abby's rainfall dataset
# MRF  =  Cell.AF_Maps[c(1:1116),]

# Using ONI and the combo dataset 1950 (row 361) - current
# MRF = MRF100[c(361:(nrow(MRF100)-12)),]
MRF100[361,]
MRF = MRF100[361:nrow(MRF100),]
head(MRF)
tail(MRF)

##########  Need to remove Rows to get seasons correct
# first row always starts on May 1950
MRF2 =  MRF[-c(1:4),]
head(MRF2)

# last row changes depending on dataset end point
# get last row in dataset
lr<-MRF2[nrow(MRF2),]
lr

# get count of how many months in the last year
yc<-max(MRF2$Year)
yc<-nrow(MRF2[which(MRF2$Year == yc),])
yc

# get row number of April or October of last year
mrn<-data.frame(which(MRF2$Year == lr$Year))
mrn

if(nrow(mrn>4)) {arn<-mrn[4,]}
if(nrow(mrn>10)) {orn<-mrn[10,]}

# make dataset end at the end of the last complete season
if(yc<4){MRF2<-MRF2[which(MRF2$Year < lr$Year),]} 
if(yc>4 && yc<10) {MRF2<-MRF2[1:arn,]}
if(yc>10) {MRF2<-MRF2[1:orn,]}

tail(MRF2)

MRF3 <- MRF2
MRF3$Month <- as.numeric(MRF3$Month)
debug_print("98.9, MEI")

### for each consecutive 6 months, aggregate ANOM by season and keep maximum
# make list of values 6 values apart
rows <- seq(from = 1, to = nrow(MRF3), by = 6)
rows

# create empty dataframe
seasons <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), 
                    c("Year", "RF", "season"))
head(seasons)
head(MRF3)

# loop through consecutive months and aggregate season oni2 values

n <- 1

for (y in rows) {

  b<-MRF3[c(y:(y+5)),]
  if(!is.na(mean(b$RF))) {
    # calculate average seasonal rainfall value and set season
    seasons[n,]$RF<-mean(b$RF)
    seasons[n,]$Year<-max(b$Year)
    if(as.numeric(min(b$Month) == 5)) {seasons[n,]$season<-"dry"}
    if(as.numeric(min(b$Month) == 1)) {seasons[n,]$season<-"wet"}

    n<-n+1
  }
}

head(seasons, 10)
tail(seasons)
summary(seasons$RF)

# # selects every other value since they alternate (dry, wet, dry, etc. Starts with dry (May))
# DryRF <- SixMo[c(TRUE, FALSE)]
# WetRF <- SixMo[c(FALSE, TRUE)]

# DryRF<-seasons[which(seasons$season == "dry"),]$RF
# WetRF<-seasons[which(seasons$season == "wet"),]$RF


##########   Bind MEI and Data
debug_print("98.12, Bind MEI and Data")

MEI_D
DryRF <- seasons[which(seasons$season == "dry"), ]
MEI_Dly <- max(MEI_D$Year)

DryRF
# DryRFly <- DryRF[which(DryRF$Year <= MEI_Dly), ]

if(nrow(MEI_D[which(!is.na(MEI_D$MEI_D)), , drop = FALSE]) != nrow(DryRF))
{DryRF<-DryRF[1:nrow(MEI_D[which(!is.na(MEI_D$MEI_D)), , drop = FALSE]), ]}

L0_D<-cbind(DryRF,MEI_D[which(!is.na(MEI_D$MEI_D)),])
L0_D

MEI_W
MEI_W[which(!is.na(MEI_W$MEI_W)), ] 

WetRF <- seasons[which(seasons$season == "wet"), ]
WetRF

if(nrow(MEI_W[which(!is.na(MEI_W$MEI_W)), , drop = FALSE]) != nrow(WetRF))
{WetRF<-WetRF[1:nrow(MEI_W[which(!is.na(MEI_W$MEI_W)), , drop = FALSE]), ]}

L0_W <- cbind(WetRF, MEI_W[which(!is.na(MEI_W$MEI_W)), ])
L0_W

##########   Wet Season
##########   Separate by ENSO Phase
debug_print("98.13, Wet Season")

ELN_W  <- subset(L0_W, MEI_W > 0.5)
LAN_W  <- subset(L0_W, MEI_W < -0.5)
NUT_Wx <- subset(L0_W, MEI_W > -0.5)
NUT_W  <- subset(NUT_Wx, MEI_W < 0.5)

##########   Strong and Weak EL Wet Season
debug_print("98.14, Strong and Weak EL Wet Season")

ELN_W_Weak <- subset(ELN_W , MEI_W  <= 1.5)
ELN_W_Strong <- subset(ELN_W , MEI_W  > 1.5)

##########   Strong and Weak La Wet Season
debug_print("98.15, Strong and Weak La Wet Season")

LAN_W_Weak <- subset(LAN_W , MEI_W  >= -1.5)
LAN_W_Strong <- subset(LAN_W , MEI_W  < -1.5)

##########   DRY Season
debug_print("98.16, Dry Season")

ELN_D <- subset(L0_D, MEI_D > 0.5)
LAN_D <- subset(L0_D, MEI_D < -0.5)
NUT_Dx <- subset(L0_D, MEI_D > -0.5)
NUT_D  <- subset(NUT_Dx, MEI_D < 0.5)

##########   Strong and Weak EL Dry Season
debug_print("98.17, Strong and Weak EL Dry Season")

ELN_D_Weak <- subset(ELN_D , MEI_D  <= 1.5)
ELN_D_Weak
ELN_D_Strong <- subset(ELN_D , MEI_D  > 1.5)
ELN_D_Strong

##########   Strong and Weak La Dry Season
debug_print("98.18, Strong and Weak La Dry Season")

LAN_D_Weak <- subset(LAN_D , MEI_D  >= -1.5)
LAN_D_Strong <- subset(LAN_D , MEI_D  < -1.5)

##########   Extract RF values
debug_print("98.19, Extract RF values")

EL_W_S <- ELN_W_Strong[, 2]
EL_W_W <- ELN_W_Weak[, 2]
LA_W_S <- LAN_W_Strong[, 2]
LA_W_W <- LAN_W_Weak[, 2]
NU_W <- NUT_W[, 2]

EL_D_S <- ELN_D_Strong[, 2]
EL_D_W <- ELN_D_Weak[, 2]
LA_D_S <- LAN_D_Strong[, 2]
LA_D_W <- LAN_D_Weak[, 2]
NU_D <- NUT_D[, 2]

##########   Counting the number in each phase
debug_print("98.20, Counting the number in each phase")

C_EL_W_S <- sum(!is.na(EL_W_S))
C_EL_W_W <- sum(!is.na(EL_W_W))
C_LA_W_S <- sum(!is.na(LA_W_S))
C_LA_W_W <- sum(!is.na(LA_W_W))
C_NU_W <-   sum(!is.na(NU_W))
C_EL_D_S <- sum(!is.na(EL_D_S))
C_EL_D_W <- sum(!is.na(EL_D_W))
C_LA_D_S <- sum(!is.na(LA_D_S))
C_LA_D_W <- sum(!is.na(LA_D_W))
C_NU_D <-   sum(!is.na(NU_D))

Cell.MEI[1:5, 8] <- c(C_EL_W_S, C_EL_W_W, C_NU_W, C_LA_W_W, C_LA_W_S)
Cell.MEI[1:5, 9] <- c(C_EL_D_S, C_EL_D_W, C_NU_D, C_LA_D_W, C_LA_D_S)
Cell.MEI
##########   Mean RF values for each season-phase
debug_print("98.21, Mean RF values for each season-phase")

Me_EL_W_S <- round(mean(EL_W_S, na.rm = T), 1)
Me_EL_W_W <- round(mean(EL_W_W, na.rm = T), 1)
Me_LA_W_S <- round(mean(LA_W_S, na.rm = T), 1)
Me_LA_W_W <- round(mean(LA_W_W, na.rm = T), 1)
Me_NU_W <-   round(mean(NU_W, na.rm = T), 1)
Me_EL_D_S <- round(mean(EL_D_S, na.rm = T), 1)
Me_EL_D_W <- round(mean(EL_D_W, na.rm = T), 1)
Me_LA_D_S <- round(mean(LA_D_S, na.rm = T), 1)
Me_LA_D_W <- round(mean(LA_D_W, na.rm = T), 1)
Me_NU_D <-   round(mean(NU_D, na.rm = T), 1)

Cell.MEI[1:5, 2] <- c(Me_EL_W_S, Me_EL_W_W, Me_NU_W, Me_LA_W_W, Me_LA_W_S)
Cell.MEI[1:5, 3] <- c(Me_EL_D_S, Me_EL_D_W, Me_NU_D, Me_LA_D_W, Me_LA_D_S)

##########   MAX
debug_print("98.22, MAX")

Mx_EL_W_S <- round(max(EL_W_S, na.rm = T), 1)
Mx_EL_W_W <- round(max(EL_W_W, na.rm = T), 1)
Mx_LA_W_S <- round(max(LA_W_S, na.rm = T), 1)
Mx_LA_W_W <- round(max(LA_W_W, na.rm = T), 1)
Mx_NU_W <-   round(max(NU_W, na.rm = T), 1)
Mx_EL_D_S <- round(max(EL_D_S, na.rm = T), 1)
Mx_EL_D_W <- round(max(EL_D_W, na.rm = T), 1)
Mx_LA_D_S <- round(max(LA_D_S, na.rm = T), 1)
Mx_LA_D_W <- round(max(LA_D_W, na.rm = T), 1)
Mx_NU_D <-   round(max(NU_D, na.rm = T), 1)

Cell.MEI[1:5, 4] <- c(Mx_EL_W_S, Mx_EL_W_W, Mx_NU_W, Me_LA_W_W, Mx_LA_W_S)
Cell.MEI[1:5, 5] <- c(Mx_EL_D_S, Mx_EL_D_W, Mx_NU_D, Me_LA_D_W, Mx_LA_D_S)

##########   MIN
debug_print("98.23, MIN")

Mn_EL_W_S <- round(min(EL_W_S, na.rm = T), 1)
Mn_EL_W_W <- round(min(EL_W_W, na.rm = T), 1)
Mn_LA_W_S <- round(min(LA_W_S, na.rm = T), 1)
Mn_LA_W_W <- round(min(LA_W_W, na.rm = T), 1)
Mn_NU_W <-   round(min(NU_W, na.rm = T), 1)
Mn_EL_D_S <- round(min(EL_D_S, na.rm = T), 1)
Mn_EL_D_W <- round(min(EL_D_W, na.rm = T), 1)
Mn_LA_D_S <- round(min(LA_D_S, na.rm = T), 1)
Mn_LA_D_W <- round(min(LA_D_W, na.rm = T), 1)
Mn_NU_D <-   round(min(NU_D, na.rm = T), 1)

Cell.MEI[1:5, 6] <- c(Mn_EL_W_S, Mn_EL_W_W, Mn_NU_W, Mn_LA_W_W, Mn_LA_W_S)
Cell.MEI[1:5, 7] <- c(Mn_EL_D_S, Mn_EL_D_W, Mn_NU_D, Mn_LA_D_W, Mn_LA_D_S)

# ##########   SUM
# Su_EL_W_S <- round(sum(EL_W_S,na.rm=T),1)
# Su_EL_W_W <- round(sum(EL_W_W,na.rm=T),1)
# Su_LA_W_S <- round(sum(LA_W_S,na.rm=T),1)
# Su_LA_W_W <- round(sum(LA_W_W,na.rm=T),1)
# Su_NU_W <-   round(sum(NU_W,na.rm=T),1)
# Su_EL_D_S <- round(sum(EL_D_S,na.rm=T),1)
# Su_EL_D_W <- round(sum(EL_D_W,na.rm=T),1)
# Su_LA_D_S <- round(sum(LA_D_S,na.rm=T),1)
# Su_LA_D_W <- round(sum(LA_D_W,na.rm=T),1)
# Su_NU_D <-   round(sum(NU_D,na.rm=T),1)
#
# Cell.MEI[1:5,8] <- c( Su_EL_W_S, Su_EL_W_W, Su_NU_W, Su_LA_W_W, Su_LA_W_S)
# Cell.MEI[1:5,9] <- c( Su_EL_D_S, Su_EL_D_S, Su_NU_D, Su_LA_D_W, Su_LA_D_S)

Cell.MEI

##########   DRY SEASON
debug_print("99, MEI, DRY SEASON")

##########   Scaler
MAXA <- max(c(Mx_EL_D_S, Mx_EL_D_W, Mx_LA_D_S, Mx_LA_D_W, Mx_NU_D))
MAXB <- MAXA * 0.05
MAXB2 <- MAXA * 0.10
MAXC <- MAXA * 0.15
MAXD <- (MAXC * 1.2) + MAXA
MAXD2 <- MAXD - MAXB2
MAXD22 <- MAXD - MAXB
MAXD3 <- MAXD - MAXC

dpi = 300

png(
  paste0(PATH_WITH_PROJECT_NAME, "MEI_DRY.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

DRY_ENSO <- c(EL_D_S, EL_D_W, NU_D, LA_D_S, LA_D_W)

boxplot(c(EL_D_S),(EL_D_W),(NU_D),(LA_D_W),(LA_D_S),col=c("darkred","red","grey","blue","darkblue"),
  names = c("SEL" ,"WEL","NUT","WLA", "SLA"),ylab = paste0("Avg. Monthly Rainfall (",RFUnit2,")"),
  ylim=c(0,MAXD))


title(paste0("Dry Sea. RF: ",UNIT_Ns[u]), line = 1,cex.main = 1.1,col.main="gold4")

text(1, MAXD, paste0("Count = ", C_EL_D_S), cex = 0.7)
text(1, MAXD2, paste0("Mean = ", Me_EL_D_S), cex = 0.7)
text(1, MAXD3, paste0("Max = ", Mx_EL_D_S), cex = 0.7)
text(1, MAXD22, paste0("Min = ", Mn_EL_D_S), cex = 0.7)

text(2, MAXD, paste0("Count = ", C_EL_D_W), cex = 0.7)
text(2, MAXD2, paste0("Mean = ", Me_EL_D_W), cex = 0.7)
text(2, MAXD3, paste0("Max = ", Mx_EL_D_W), cex = 0.7)
text(2, MAXD22, paste0("Min = ", Mn_EL_D_W), cex = 0.7)

text(3, MAXD, paste0("Count = ", C_NU_D), cex = 0.7)
text(3, MAXD2, paste0("Mean = ", Me_NU_D), cex = 0.7)
text(3, MAXD3, paste0("Max = ", Mx_NU_D), cex = 0.7)
text(3, MAXD22, paste0("Min = ", Mn_NU_D), cex = 0.7)

text(4, MAXD, paste0("Count = ", C_LA_D_W), cex = 0.7)
text(4, MAXD2, paste0("Mean = ", Me_LA_D_W), cex = 0.7)
text(4, MAXD3, paste0("Max = ", Mx_LA_D_W), cex = 0.7)
text(4, MAXD22, paste0("Min = ", Mn_LA_D_W), cex = 0.7)

text(5, MAXD, paste0("Count = ", C_LA_D_S), cex = 0.7)
text(5, MAXD2, paste0("Mean = ", Me_LA_D_S), cex = 0.7)
text(5, MAXD3, paste0("Max = ", Mx_LA_D_S), cex = 0.7)
text(5, MAXD22, paste0("Min = ", Mn_LA_D_S), cex = 0.7)

dev.off()

##########   Wet Season
##########   Scaler
debug_print("100, MEI, Wet SEASON")

MAXA <- max(c(Mx_EL_W_S, Mx_EL_W_W, Mx_LA_W_S, Mx_LA_W_W, Mx_NU_W))
MAXB <- MAXA * 0.05
MAXB2 <- MAXA * 0.10
MAXC <- MAXA * 0.15
MAXD <- (MAXC * 1.2) + MAXA
MAXD2 <- MAXD - MAXB2
MAXD22 <- MAXD - MAXB
MAXD3 <- MAXD - MAXC

png(
  paste0(PATH_WITH_PROJECT_NAME, "MEI_WET.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

boxplot(c(EL_W_S),(EL_W_W),(NU_W),(LA_W_W),(LA_W_S),col=c("darkred","red","grey","blue","darkblue"),
  names = c("SEL" ,"WEL","NUT","WLA", "SLA"),ylab=paste0("Avg. Monthly Rainfall (",RFUnit2,")"),
  ylim=c(0,MAXD))
title(paste0("Wet Sea RF: ",UNIT_Ns[u]), line = 1,cex.main = 1.1,col.main="darkgreen")

text(1, MAXD, paste0("Count = ", C_EL_W_S), cex = 0.7)
text(1, MAXD2, paste0("Mean = ", Me_EL_W_S), cex = 0.7)
text(1, MAXD3, paste0("Max = ", Mx_EL_W_S), cex = 0.7)
text(1, MAXD22, paste0("Min = ", Mn_EL_W_S), cex = 0.7)

text(2, MAXD, paste0("Count = ", C_EL_W_W), cex = 0.7)
text(2, MAXD2, paste0("Mean = ", Me_EL_W_W), cex = 0.7)
text(2, MAXD3, paste0("Max = ", Mx_EL_W_W), cex = 0.7)
text(2, MAXD22, paste0("Min = ", Mn_EL_W_W), cex = 0.7)

text(3, MAXD, paste0("Count = ", C_NU_W), cex = 0.7)
text(3, MAXD2, paste0("Mean = ", Me_NU_W), cex = 0.7)
text(3, MAXD3, paste0("Max = ", Mx_NU_W), cex = 0.7)
text(3, MAXD22, paste0("Min = ", Mn_NU_W), cex = 0.7)

text(4, MAXD, paste0("Count = ", C_LA_W_W), cex = 0.7)
text(4, MAXD2, paste0("Mean = ", Me_LA_W_W), cex = 0.7)
text(4, MAXD3, paste0("Max = ", Mx_LA_W_W), cex = 0.7)
text(4, MAXD22, paste0("Min = ", Mn_LA_W_W), cex = 0.7)

text(5, MAXD, paste0("Count = ", C_LA_W_S), cex = 0.7)
text(5, MAXD2, paste0("Mean = ", Me_LA_W_S), cex = 0.7)
text(5, MAXD3, paste0("Max = ", Mx_LA_W_S), cex = 0.7)
text(5, MAXD22, paste0("Min = ", Mn_LA_W_S), cex = 0.7)

dev.off()

Cell.MEI
# I don't think this is used...
write.csv(Cell.MEI,
  paste0(PATH_WITH_PROJECT_NAME, "MEI_A.csv"),
  row.names = F)



#####################################################
##### ENSO rainfall barplots (Derek's addition)
debug_print("101, ENSO rainfall barplots")

# Average Monthly Rainfall by ENSO Phase and Season barplot
# spell out the ENSO phases
Cell.MEI2 <- Cell.MEI
Cell.MEI2$Phase2 <- NA
Cell.MEI2

p <- as.list(unique(Cell.MEI2$Phase))
p

for (i in 1:nrow(Cell.MEI2)) {
  y<-Cell.MEI2[i,]
  if(y$Phase == p[[1]]) {Cell.MEI2[i,]$Phase2 <- "Strong El Nino"}
  if(y$Phase == p[[2]]) {Cell.MEI2[i,]$Phase2 <- "Weak El Nino"}
  if(y$Phase == p[[3]]) {Cell.MEI2[i,]$Phase2 <- "Neutral"}
  if(y$Phase == p[[4]]) {Cell.MEI2[i,]$Phase2 <- "Weak La Nina"}
  if(y$Phase == p[[5]]) {Cell.MEI2[i,]$Phase2 <- "Strong La Nina"}
}

Cell.MEI2

# reformat table for barplots

colnames(Cell.MEI2)
cd <- subset(Cell.MEI2, select = c("D-Mean", "D-Count", "Phase2"))
colnames(cd) <- c("Mean", "Count", "Phase")
cd$Season <- "dry"
cw <- subset(Cell.MEI2, select = c("W-Mean", "W-Count", "Phase2"))
colnames(cw) <- c("Mean", "Count", "Phase")
cw$Season <- "wet"

cc <- rbind(cd, cw)
cc

### Plot
# set order of seasons for plotting
cc$Season <- factor(cc$Season, levels = c("wet", "dry"))

# set order of ENSO phases for plotting
cc$Phase<-factor(cc$Phase, levels=c("Strong El Nino","Weak El Nino","Neutral","Weak La Nina","Strong La Nina"))

cc

# set ylim
ylim <- max(cc$Mean) * 1.2

png(
  paste0(PATH_WITH_PROJECT_NAME, "ENSO_season_barplot.png"),
  width = 6 * dpi,
  height = 2.5 * dpi,
  res = dpi
)

ggplot(data=cc, 
  aes(x=Season, y=Mean, group=Phase)) +
  geom_bar(aes(fill=Phase), position = position_dodge(width=0.7), stat="identity", color="black", 
    alpha=.7, width=.55) +
  labs(title="Average Monthly Rainfall by Season and ENSO Phase",
    y = "Rainfall (inches)", x= "Season") +
  scale_fill_manual(values=c("darkred","red","grey","lightskyblue1","darkblue"),
    limits=c("Strong El Nino","Weak El Nino","Neutral","Weak La Nina","Strong La Nina")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  guides(fill=guide_legend(title="ENSO Phase")) +
  geom_text(aes(label=Count), position=position_dodge(width=0.7), vjust=-0.8) +
  theme_bw()+
  theme(axis.text.x=element_text(size=13),
    axis.text.y=element_text(size=13),
    axis.title.x=element_text(size=13),
    axis.title.y=element_text(size=13))

dev.off()

cc
### Export table of average monthly rainfall by season and ENSO phase
write.csv(cc, paste0(PATH_WITH_PROJECT_NAME, "MEI_S.csv"))

#####################################################
##### Air temperature graph
debug_print("102, Air temperature graph")

########## Month Year Air Temperature Maps
#AT_Map_Path_A <- ("D:/PDKE/CCVD/CCVD INPUTS/air_temp/data_map_newer/")
AT_Map_Path_A <- paste0(INPUTS_FOLDER, "/air_temp/data_map_newer/")
maps <- AT_Map_Path_A

# make list of years
years <- as.list(list.files(maps))

# create empty table for raster stats
table <- data.frame(matrix(ncol = 5, nrow = 0, 
                           dimnames = list(NULL, c("year", "month", "min", "max", "mean"))))

# set row
no <- 1

### for each year, loop through months and extract min, max, mean air temp value
### from within the study area

for (i in years) {
  # i="2019"
  setwd(paste0(maps, i, "/"))
  # make list of monthly files
  f <- as.list(list.files(), pattern = ".tif")
  # loop through monthly files and fill in table
  for (x in f) {
    # x=f[[1]]
    # set year and month
    y <- substr(x, nchar(x) - 10, nchar(x))
    table[no, ]$year <- substr(y, 1, 4)
    m <- substr(x, nchar(x) - 5, nchar(x))
    table[no, ]$month <- substr(m, 1, 2)
    # load raster
    map <- raster(x)
    # crop to study area
    m2 <- crop(map, extent(HALE))
    m3 <- mask(m2, HALE)
    m3 <- reclassify(m3, cbind(-Inf, 0, NA), right = TRUE)
    summary(m3)
    
    # calculate min, max, mean air temp
    table[no, ]$min <- round(cellStats(m3, stat = "min"), digits = 2)
    table[no, ]$max <- round(cellStats(m3, stat = "max"), digits = 2)
    table[no, ]$mean <- round(cellStats(m3, stat = "mean"), digits = 2)
    
    no <- no + 1
  }
}

head(table)
summary(table$min)

###### Monthly air temperature time series
debug_print("103, Monthly air temperature time series")

dat <- table

head(dat)
tail(dat)

dat[which(dat$max == max(dat$max)), ]

# convert all celcius to farenheit
dat$min <- (dat$min * (9 / 5)) + 32
dat$max <- (dat$max * (9 / 5)) + 32
dat$mean <- (dat$mean * (9 / 5)) + 32

# date column
dat$date <- as.Date(paste0(dat$year, "-", dat$month, "-01"))

# only keep complete years
dat$full <- c("N")
head(dat)

y <- as.list(unique(dat$year))

for (i in y) {
  if (length(which(dat$year == i)) == 12) { dat[which(dat$year == i), ]$full <- "Y" }
}

tail(dat, 20)

ic <- unique(dat[which(dat$full == "N"), ]$year)


# dat<-dat[which(dat$full == "Y"),]
head(dat)
tail(dat)
#############################
### Annual air temp trend
debug_print("104, Annual air temp trend")

# aggregate monthly to annual air temp
rm(min)
rm(max)
rm(mean)

dat.y <- aggregate(mean ~ year, dat, mean)
dat.y$mean <- round(dat.y$mean, digit = 2)
dat.y$meanmin <- round((aggregate(mean ~ year, dat, min))$mean, digits = 2)
dat.y$meanmax <- round((aggregate(mean ~ year, dat, max))$mean, digits = 2)
dat.y$min <- round((aggregate(min ~ year, dat, min))$min, digits = 2)
dat.y$max <- round((aggregate(max ~ year, dat, max))$max, digits = 2)

dat.y

# remove stats from incomplete year
if (length(ic) != 0) {
  dat.y <- dat.y[which(dat.y$year < ic), ]
  dat.y[(nrow(dat.y) + 1), ]$year <- ic
}

# make date column
dat.y$date <- as.Date(paste0(dat.y$year, "-01-01"))

# write to csv table of annual air temp values from monthly means (min, max, mean)
write.csv(dat.y, paste0(PATH_WITH_PROJECT_NAME, "monthly_airtemp.csv"))

slope <- formatC((coef(lm(dat.y$mean ~ dat.y$date))[2]), format = "e", digits = 2)
slope

dpi = 300

tail(dat.y)
png(
  paste0(PATH_WITH_PROJECT_NAME, "annual_airtemp.png"),
  width = 6 * dpi,
  height = 3 * dpi,
  res = dpi
)

ggplot(dat.y, aes(x=date, y=mean)) +
  geom_line(size=1.2, color="orange") +
  geom_smooth(method=lm, se=F, size=1, color="black") +
  stat_cor(method="pearson", label.x=as.Date("1991-01-01"), label.y=58) +
  scale_x_date(date_breaks = "5 years", labels = date_format(format="%Y")) +
  labs(title="Annual Air Temperature and Extremes",
    y = "Temperature (F)", x = "Year") +
  geom_line(aes(x=date,y=min), size=1.2, color="blue") +
  geom_line(aes(x=date,y=max), size=1.2, color="red") +
  annotate("text", x=as.Date("1994-07-01"), y=54,
    label=paste0("Slope = ",slope)) +
  theme(panel.background=element_rect(fill=NA, color="black"),
    panel.grid.major=element_line(color="grey90"),
    panel.grid.minor=element_blank())

dev.off()

##### add annual mean values
dat2 <- merge(dat, dat.y[, c("year", "mean", "min", "max")], by = "year", all.x = TRUE)
head(dat2)
tail(dat2)

# write to csv table of monthly air temp values (min, max, mean)
write.csv(dat2, paste0(PATH_WITH_PROJECT_NAME, "daily_airtemp.csv"))

# set y-axis limits
ylow <- min(dat.y$min, na.rm = T) * .5
yhi <- max(dat.y$max, na.rm = T) * 1.0005

# set slope
slope <- formatC((coef(lm(dat2$mean.x ~ dat2$date))[2]), format = "e", digits = 2)
slope

# set location for linear trend values
yl <- min(dat.y$min, na.rm = T) * .75

dpi = 300

png(
  paste0(PATH_WITH_PROJECT_NAME, "monthly_airtemp.png"),
  width = 6 * dpi,
  height = 3 * dpi,
  res = dpi
)

ggplot(dat2, aes(x=date,y=mean.x)) +
  geom_line(color="grey60") +
  geom_smooth(span=0.2, se=F, size=1.3, color="orange") +
  geom_smooth(method=lm, se=F, color="black") +
  stat_cor(method="pearson", label.x=as.Date("1993-01-01"), label.y=yl-5) +
  scale_x_date(date_breaks = "4 years", labels = date_format(format="%Y")) +
  ylim(ylow,yhi) +
  labs(title="Air Temperature Trends",
    y="Temperature (F)", x="") +
  geom_hline(yintercept=0) +
  annotate("text", x=as.Date("1996-07-01"), y=yl, 
    label=paste0("Slope = ",slope)) +
  geom_line(data = dat.y, aes(x=date, y=min), size=1.2, color="blue") +
  geom_line(data = dat.y, aes(x=date,y=max), size=1.2, color="red") +
  theme(panel.background=element_rect(fill=NA, color="black"),
    panel.grid.major=element_line(color="grey90"),
    panel.grid.minor=element_blank())

dev.off()

######## Air Temp Anomalies ########
#setwd("D:/PDKE/CCVD/MINI_Phase2/")
debug_print("105, Air Temp Anomalies")

setwd(WORKING_DIR)

# make list of climatology files
list <- list.files(AT_CLIM_Path_A, pattern = "\\.tif$")

### Make dataframe of average monthly climatology values within study area
# create empty table for mean AT values
clim <- data.frame(matrix(ncol = 2, nrow = 0,
  dimnames = list(NULL, c("month", "mean"))))
# set row
no <- 1
### for each monthly climatology calculate mean value within study area
for (x in list) {
  # set month
  m <- substring(x, 15, nchar(x))
  m <- as.numeric(sub(".tif.*", "", m))
  clim[no, ]$month <- m
  # load raster
  map <- raster(paste0(AT_CLIM_Path_A, x))
  # crop to study area
  m2 <- crop(map, extent(HALE))
  m3 <- mask(m2, HALE)
  # calculate mean value
  clim[no, ]$mean <- round(cellStats(m3, stat = "mean"), digits = 2)
  
  no <- no + 1
}
clim <- clim[order(clim$month), ]
clim$mean <- (clim$mean * (9 / 5)) + 32
clim

### Make dataframe of difference between month-year average air temp and that month's climatology
head(table)
anom <- data.frame()

# for each row (month-year) calculate the difference from climatology (month - climatology)
for (i in 1:nrow(table)) {
  r <- table[i, ]
  r$mean_f <- (r$mean * (9 / 5)) + 32
  m <- as.numeric(r$month)
  c <- clim[which(clim$month == m), ]$mean
  r$clim_f <- c
  r$anom_f <- r$mean_f - c
  anom <- rbind(anom, r)
}

head(anom)
tail(anom)

# make date column
anom$date <- as.Date(paste0(anom$year, "/", anom$month, "/", "01"), format = "%Y/%m/%d")

# write to csv table of monthly air temp anomaly values
write.csv(anom, paste0(PATH_WITH_PROJECT_NAME, "anomaly_airtemp.csv"))

# set y-axis limits
ylow <- min(anom$anom_f) * .95
yhi <- max(anom$anom_f) * 1.0005

# set slope
slope <- formatC((coef(lm(anom$mean_f ~ anom$date))[2]), format = "e", digits = 2)
slope

# set location for linear trend values
yl <- (((min(anom$anom_f) - max(anom$anom_f)) / 2) + max(anom$anom_f))

dpi = 300

png(
  paste0(PATH_WITH_PROJECT_NAME, "monthly_airtemp_anomaly.png"),
  width = 6 * dpi,
  height = 3 * dpi,
  res = dpi
)

ggplot(anom, aes(x=date,y=anom_f)) +
  geom_line(color="grey60") +
  geom_smooth(span=0.2, se=F, size=1.3, color="orange") +
  geom_smooth(method=lm, se=F, color="black") +
  # stat_cor(method="pearson", label.x=as.Date("1993-01-01"), label.y=yl-3) +
  scale_x_date(date_breaks = "4 years", labels = date_format(format="%Y")) +
  ylim(ylow,yhi) +
  labs(title="Air Temperature Anomaly Trends",
    y="Temperature (F)", x="") +
  # geom_hline(yintercept=0) +
  # annotate("text", x=as.Date("1996-07-01"), y=yl,
  #          label=paste0("Slope = ",slope)) +
  # geom_line(data = dat.y, aes(x=date, y=min), size=1.2, color="blue") +
  # geom_line(data = dat.y, aes(x=date,y=max), size=1.2, color="red") +
  theme(panel.background=element_rect(fill=NA, color="black"),
    panel.grid.major=element_line(color="grey90"),
    panel.grid.minor=element_blank())

dev.off()

end_time <- Sys.time()
debug_print(paste0("start: ", format(start_time, "%Y-%m-%d_%H-%M-%S")))
debug_print(paste0("end: ", format(end_time, "%Y-%m-%d_%H-%M-%S")))
debug_print(paste0("Execution time: ", end_time - start_time))

if (ENV_TYPE == "linux") {
  run_string <- paste0(RSCRIPT_PATH, " ", MYSCRIPT_PATH, " ",
   shQuote(email), " ",
   shQuote(PATH), " ",
   shQuote(NM), " ",
   shQuote(NM_s), " ",
   shQuote(NP_FILE))
  debug_print(paste0("runString: ", run_string))
  system(run_string, wait = FALSE)
}
if (ENV_TYPE == "windows") {
  # not actually used, but useful to print out in case we need to re-run via command line
  run_string <- paste0(RSCRIPT_PATH, " ", MYSCRIPT_PATH, " ",
    shQuote(email), " ",
    shQuote(PATH), " ",
    shQuote(NM), " ",
    shQuote(NM_s), " ",
    shQuote(NP_FILE))
  debug_print(paste0("runString: ", run_string))
  
  args <- c(
    MYSCRIPT_PATH,
    email,
    file.path(PATH),
    NM,
    NM_s,
    NP_FILE
  )
  
  print("run ppt")
  result <- processx::run(
    command = RSCRIPT_PATH,  # The Rscript or other command
    args = args,            # Arguments as a character vector
    echo = TRUE
  )
  
  cat(result$stdout)
  if (result$stderr != "") {
    cat("Error:", result$stderr)
  }
}

### END! ###