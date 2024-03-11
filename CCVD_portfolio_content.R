

Sys.setenv(PROJ_LIB = "/Users/jgeis/opt/anaconda3/lib/python3.8/site-packages/fiona/proj_data")
PROJ_DEBUG = 3
#install.packages("maptools", repos = "https://packagemanager.posit.co/cran/2023-10-13")
#install.packages("rgdal", repos = "https://packagemanager.posit.co/cran/2023-10-13")
#install.packages("rgeos", repos = "https://packagemanager.posit.co/cran/2023-10-13")

packages <-
  c(
    "gstat",
    "raster",
    "sp" ,
    "maptools",
    "rgdal",
    "RColorBrewer",
    "gridExtra",
    "ggplot2",
    "grid",
    "data.table",
    "devtools",
    "DescTools",
    "lubridate",
    "rgeos",
    "latticeExtra",
    "rasterVis",
    "plotrix",
    "plyr",
    "dplyr",
    "xts",
    "timeSeries",
    "ggfortify",
    "changepoint",
    "scales",
    "reshape2",
    "hydroTSM",
    "tiff",
    "lmomco",
    "parallel",
    "SPEI",
    "sf",
    "ggpubr",
    "terrainr",
    "ggmap",
    "ggthemes",
    "zoo"
  )
for (package in packages) {
  print(paste("pkgTest: ", package))
  # require() returns TRUE if the packages is loaded, FALSE if it isn't.
  if (!require(package, character.only = TRUE)) {
    print(paste("  ", package, "1"))
    install.packages(package, dep = TRUE, repos = 'http://cran.us.r-project.org')
  }
  else if (!require(package, character.only = TRUE)) {
    print(paste("  ", package, "2"))
    install.packages(package, dep = TRUE, repos = 'https://packagemanager.posit.co/cran/2023-10-13')
  }
  else if (!require(package, character.only = TRUE)) {
    print(paste("  ", package, "3"))
    install.packages(package, dep = TRUE, repos = 'http://cran.rstudio.com/')
  }
  else if (!require(package, character.only = TRUE)) {
    print(paste("Package: ", package, " not found"))
  }
  if (require(package, character.only = TRUE)) {
    print(paste("Package: ", package, " installed"))
    library(package, character.only = TRUE, quietly = TRUE)
    
  }
}

##########   This code will the generate the inputs for a CCVD portfolio
##########   Required Inputs: 1. Folder destinations 2. Measurement units 3. Shapefiles 4. List of files and output names
##########   Searchable sections of the Code:
##########   Maps - Elevation - Mean Climate - Downscaling - Rainfall Extract - SPI 1 - SPI 2 - MEI

##########################################################################################################################

BASE_DIR <- "/Users/jgeis/Work/PDKE"
WORKING_DIR <- paste0(BASE_DIR, "/CCVD/MINI_Phase2/")
setwd(WORKING_DIR)               # WORKING DIRECTORY
INPUTS_FOLDER <-
  paste0(BASE_DIR, "/CCVD/CCVD_INPUTS/")       # INPUT FOLDER
OUTPUTS_FOLDER <-
  paste0(BASE_DIR, "/CCVD/CCVD_OUTPUTS/")       # OUTPUT FOLDER
print(paste("PDKE: 1,INPUTS_FOLDER: ", INPUTS))
print(paste("PDKE: 1,OUTPUTS_FOLDER: ", OUTPUTS_FOLDER))



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

COAST_PATH <-
  paste0(INPUTS_FOLDER, "COAST/Coast_2/coast_geo_shp.dbf")
Coast <- readOGR(COAST_PATH)
print(paste("PDKE: 1,COAST_PATH", COAST_PATH))
plot(Coast[which(Coast$COAST_GEO_ == 10), ])

# ALL ISLANDS
Coast_Crop <-
  crop(Coast , extent(-160.0,-154.8066, 18.91069, 22.23))
Coast_Crop_T <- spTransform(Coast_Crop, crs(EXAMP))
# BIG ISLAND
Coast_Crop <- Coast[which(Coast$COAST_GEO_ == 14), ]
Coast_BI <- spTransform(Coast_Crop, crs(EXAMP))
# MAUI
Coast_Crop <- Coast[which(Coast$COAST_GEO_ == 10), ]
Coast_MN <- spTransform(Coast_Crop, crs(EXAMP))
#Molokai
Coast_Crop <- Coast[which(Coast$COAST_GEO_ == 9), ]
Coast_MO <- spTransform(Coast_Crop, crs(EXAMP))
#Lanai
Coast_Crop <- Coast[which(Coast$COAST_GEO_ == 11), ]
Coast_LA <- spTransform(Coast_Crop, crs(EXAMP))
plot(Coast_LA)
#Oahu
Coast_Crop <- Coast[which(Coast$COAST_GEO_ == 5), ]
Coast_OA <- spTransform(Coast_Crop, crs(EXAMP))
#KAUI
Coast_Crop <- Coast[which(Coast$COAST_GEO_ == 2), ]
Coast_KA <- spTransform(Coast_Crop, crs(EXAMP))
#Kahoolawe
Coast_Crop <- Coast[which(Coast$COAST_GEO_ == 13), ]
Coast_KO <-
  spTransform(Coast_Crop, crs(EXAMP))

############################################################
#NP_ALL <- readOGR("F:/PDKE/CCVD/sites/HAVO_boundary_2022_Pohue.shp")
NP_DIR <- paste0(INPUTS_FOLDER, "waikiki_watershed/")
NP_FILE <- paste0(NP_DIR, "waikiki_watershed.shp")

print(paste("PDKE: 1,NP_FILE", NP_FILE))
print(file.exists(NP_FILE))
# these all work, once I had all the files, not just the shp file
NP_ALL <- readOGR(NP_FILE)
#NP_ALL <- readOGR("/Users/jgeis/Work/PDKE/CCVD/CCVD_INPUTS/waikiki_watershed/waikiki_watershed.shp")
#NP_ALL <- readOGR(dsn = NP_DIR, #Provide the directory of shp file
#                  layer = "waikiki_watershed", #Provide the name of the shp file without extension (.shp)
#                  verbose = TRUE)
#NP_ALL <- readOGR(NP_FILE, layer="waikiki_watershed", verbose = TRUE)

print("PDKE: 1")

HALE <- NP_ALL
NM <- "Waikiki Watershed"
NM_s <- "Waikiki"
ILE <- "Oahu"
ILE_s <- "OA"

plot(HALE)
print("PDKE: 2")

if (ILE_s == "BI") {
  plot(Coast_BI, main = ILE) + plot(HALE , add = T, col = "red")
}
if (ILE_s == "MN") {
  plot(Coast_MN, main = ILE) + plot(HALE , add = T, col = "red")
}
if (ILE_s == "MO") {
  plot(Coast_MO, main = ILE) + plot(HALE , add = T, col = "red")
}
if (ILE_s == "LA") {
  plot(Coast_LA, main = ILE) + plot(HALE , add = T, col = "red")
}
if (ILE_s == "KO") {
  plot(Coast_KO, main = ILE) + plot(HALE , add = T, col = "red")
}
if (ILE_s == "OA") {
  plot(Coast_OA, main = ILE) + plot(HALE , add = T, col = "red")
}
if (ILE_s == "KA") {
  plot(Coast_KA, main = ILE) + plot(HALE , add = T, col = "red")
}

print("PDKE: 3")

##################################################
##########   Define Units

##########   Assign AOI shapefile coordinates
HALE <- spTransform(HALE, crs(EXAMP))

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
print("PDKE: 4")


##########   COUNT INPUTS

#UNIT_Shape <- sum(!is.na(UNIT_X))
UNIT_C <- sum(!is.na(UNIT_N))  # This will mark the end of the Loop
#UNIT_Island <- sum(!is.na(UNIT_I))
#ShortNm <- sum(!is.na(UNIT_Ns))

##########   CONFIRM INPUTS ARE THE SAME

# UNIT_C %>%
# ShortNm
# UNIT_Shape
# UNIT_Island

##########  Mean Rainfall DATA

MeanRF_ALL <-
  dir(
    paste0(INPUTS_FOLDER, "Mean_RF_Data/StateMaps/"),
    pattern = "*x.adf",
    recursive = T,
    full.names = T
  )
#print(paste0("PDKE: 5:", MeanRF_ALL))

# ### Find state-wide mean monthly rainfall stats (min, mean, max, percentiles)
# # loop through month maps and calculate stats, also make average monthly rainfall map
# n<-1
# fq<-data.frame(quantile(c(0:0)))
# print("PDKE: 6")

# sm<-stack()
# print("PDKE: 7")

# for (i in MeanRF_ALL[1:12]) {
#   m<-raster(i)
#   print(paste0("PDKE: 8", i))
#   print(paste0("PDKE: m", m))

#   # convert mm to inches
#   if(RFUnit == " in") {m  <- m  * 0.0393701}
#   print("PDKE: 9")

#   # calculate quantiles and put in dataframe
#   quant<-quantile(m)
#   print("PDKE: 10")

#   mq<-data.frame(quant)
#   print("PDKE: 11")

#   colnames(mq)<-n
#   fq<-cbind(fq,mq)
#   n<-n+1

#   # add raster to stack
#   sm<-stack(sm, m)
#   print("PDKE: 12")
# }
# sm
# print("PDKE: 13")

# rm(mean)
# print("PDKE: 14")

#fq<-round(fq[2:ncol(fq)], 2)
#stateRF<-rowMeans(fq)
# stateRFM<-calc(sm, mean, na.rm=TRUE)
# print("PDKE: 15")
stateRFM <-
  read.csv(paste0(INPUTS_FOLDER, "Mean_RF_Data/StateMaps/stateRFMd_input.csv"))
#stateRFM<-read.csv(paste0(INPUTS_FOLDER,"Mean_RF_Data/StateMaps/stateRFM_quantiles_input.csv"))


# make dataframe from average monthly rainfall map
stateRFMd <- as.data.frame(stateRFM, na.rm = T)
stateRFMd <- data.frame(stateRFMd[1])

u <- 1
print("PDKE: 16")

##########   FIRE OCCURRENCE AND RISK
#Fire Occurrence Shape 2022 (From Clay)
#FIRE_Shape <- readOGR(paste0(INPUTS_FOLDER,"StateFire_1999/fires_1999_2022/fires_1999_2022.shp"))
FIRE_Shape_File <-
  paste0(INPUTS_FOLDER, "fires_1999_2022/fires_1999_2022.shp")
print(paste("PDKE: 17,FIRE_Shape_File", FIRE_Shape_File))
FIRE_Shape <- readOGR(FIRE_Shape_File)
print("PDKE: 18")
FIRE_Shape_T <- spTransform(FIRE_Shape, crs(EXAMP))
print("PDKE: 19")
crs(FIRE_Shape_T) <- CRS(proj4string(EXAMP))
print("PDKE: 20")

# Fire_Mask <- mask(x = FIRE_Shape_T, mask = UNIT_X[[u]])
# Fire_Crop <- crop(x = FIRE_Shape_T, y = extent(UNIT_X[[u]]))

#Fire risk categories 2023 (From Clay)
FRISK <-
  raster(
    paste0(
      INPUTS_FOLDER,
      "FireRisk/Avg_Landscape_Fire_Risk_Hawaii_2023/Avg_Landscape_Fire_Risk_Hawaii_fire_risk_categories_2023.tif"
    )
  )
crs(FRISK) <- "+proj=longlat +datum=WGS84 +no_defs"
print("PDKE: 21")

FRISK_Shape <-
  readOGR(
    paste0(
      INPUTS_FOLDER,
      "FireRisk/Avg_Landscape_Fire_Risk_Hawaii_2023/fire_risk_polygons.shp"
    )
  )
FRISK_Shape_T <- spTransform(FRISK_Shape, crs(EXAMP))
crs(FRISK_Shape_T) <- CRS(proj4string(EXAMP))
print("PDKE: 22")

##########   Forest Roads
# F_Roads <- readOGR(paste0(INPUTS_FOLDER,"Forestry Roads/Forestry_Roads.shp"))
F_Roads_File <-
  paste0(INPUTS_FOLDER, "Forestry Roads/Forestry_Roads.shp")
print(paste("PDKE: 23,F_Roads_File", F_Roads_File))
F_Roads <- readOGR(F_Roads_File)
F_Roads <- spTransform(F_Roads, crs(EXAMP))
plot(F_Roads)

##########   Trail Inventory
# I_Trails <- readOGR(paste0(INPUTS_FOLDER,"Inventory Trails/InventoryTrails.shp"))
I_Trails_File <-
  paste0(INPUTS_FOLDER, "Inventory Trails/InventoryTrails.shp")
print(paste("PDKE: 24,I_Trails_File", I_Trails_File))
I_Trails <- readOGR(I_Trails_File)
I_Trails <- spTransform(I_Trails, crs(EXAMP))
plot(I_Trails)

##########   Digital Elevation Model

ELEV2 <-
  dir(
    paste0(INPUTS_FOLDER, "ned_dem/"),
    pattern = "*x.adf",
    recursive = T,
    full.names = T
  )
# ### if doing portfolio for whole island (needed for Big Island at least)
#   ELEV2 = dir(paste0(INPUTS_FOLDER), pattern="*ned_dem_crs.tif", recursive=T, full.names=T)
ELEV <- raster(ELEV2[1])
if (ELUnit == " ft") {
  ELEV <- ELEV * 3.28084
}
crs(ELEV) <- "+proj=longlat +datum=WGS84 +no_defs"
plot(ELEV)

##########  Hillshade
HS <- raster(paste0(INPUTS_FOLDER, "hawaii_hillshade.tif"))
crs(HS) <- "+proj=longlat +datum=WGS84 +no_defs"
plot(HS)

# # convert to dataframe for plotting
# HSd<-as.data.frame(HS, xy=T)
# HSd<-na.omit(HSd)

##########  Landcover

LC2 <-
  dir(
    paste0(INPUTS_FOLDER, "Landcover/"),
    pattern = "*LCMAP_HI_2020_V10_LCPRI_crs.tif$",
    recursive = T,
    full.names = T
  )
LC <- raster(LC2)
crs(LC) <- "+proj=longlat +datum=WGS84 +no_defs"
plot(LC)

##########  Moku and Ahupuaa
MOKU <- readOGR(paste0(INPUTS_FOLDER, "Moku.shp"))
MOKU <- spTransform(MOKU, crs(EXAMP))
plot(MOKU)

AHU <- readOGR(paste0(INPUTS_FOLDER, "Ahupuaa2.shp"))
AHU <- spTransform(AHU, crs(EXAMP))
plot(AHU)

##########  Streams and Aquifers
STRM <-
  readOGR(paste0(
    INPUTS_FOLDER,
    "NHD_H_Hawaii_State_Shape/NHD_Flowlines2.shp"
  ))
STRM <- spTransform(STRM, crs(EXAMP))

AQU <-
  readOGR(paste0(INPUTS_FOLDER, "Aquifers/DOH_Aquifers_type_status.shp"))
AQU <- spTransform(AQU, crs(EXAMP))

##########  Downscaling

# should be 132 DyDs_files. If TIF's were opened in arcgis, extra files may have been created
DyDs_files <-
  dir(
    paste0(INPUTS_FOLDER, "HI_Downscaling_Final/"),
    recursive = T,
    full.names = T,
    pattern = "\\.tif$"
  )
D2_files <-
  dir(
    paste0(INPUTS_FOLDER, "Lulin_Downscaling/"),
    recursive = T,
    full.names = T
  )
print(paste("PDKE: 25,DyDs_files: ", DyDs_files))
print(paste("PDKE: 26,D2_files: ", D2_files))


##########   Month Year Rainfall Maps (Frazier et al., 2016)

# RF_Map_Path_A <- ("F:/PDKE/CCVD/data/production/rainfall/legacy/month/statewide/data_map/")
RF_Map_Path_A_DIR <-
  paste0(BASE_DIR,
         "/CCVD/data/production/rainfall/legacy/month/statewide/data_map/")
print(paste("PDKE: 27,RF_Map_Path_A_DIR: ", RF_Map_Path_A_DIR))
RF_Map_Path_A <- (RF_Map_Path_A_DIR)

##########   Month-Year Rainfall Maps New Lucas et al., In Review

# Need to download maps from HCDP
# Legacy maps 1920-1989
# Matty maps 1990-NRT
# Put in INput folder
#For now Read this path in seperate as not in input folder
#RF_Map_Path <- ("F:/PDKE/CCVD/NEW RF MAPS/statewide_rf_mm/rf_mm/")
RF_Map_Path_DIR <-
  paste0(BASE_DIR, "/CCVD/NEW_RF_MAPS/statewide_rf_mm/rf_mm/")
RF_Map_Path <- (RF_Map_Path_DIR)
print(paste("PDKE: 28,RF_Map_Path: ", RF_Map_Path))

##########  Multi-Variate ENSO Index

#Replace this with the ONI Dataset at some point
MEI <- read.csv(paste0(INPUTS_FOLDER, "ONI_Season.csv"), sep = ",")



########## Monthly Climatology Air Temperature Maps
#AT_CLIM_Path_A <- ("F:/PDKE/CCVD/CCVD INPUTS/air_temp/air_temp_climatology/")
AT_CLIM_Path_A <-
  paste0(INPUTS_FOLDER, "/air_temp/air_temp_climatology/")

#################################################################################################################################

##########   DataFrames For Output

#Set up Group matrix for all climate Variables
#Climate For all Ranches
Cell.DataCL <- data.frame(matrix(ncol = 15, nrow = UNIT_C))
colnames(Cell.DataCL) <-
  c(
    "Unit",
    "Elev Min",
    "Elev Mean",
    "Elev Max",
    "RF",
    "Tavg",
    "Tmax",
    "Tmin",
    "RH",
    "SM",
    "KD",
    "ET",
    "CF",
    "Island",
    "Short Name"
  )
Cell.DataCL[1] <- UNIT_N
Cell.DataCL[16] <- UNIT_Ns

#Climate For Individual Ranches
Cell.DataCLR <- data.frame(matrix(ncol = 16, nrow = 3))
colnames(Cell.DataCLR) <-
  c(
    "Unit",
    "Elev",
    "RF",
    "Tavg",
    "Tmax",
    "Tmin",
    "RH",
    "SM",
    "KD",
    "ET",
    "CF",
    "WS",
    "DryRF",
    "WetRF",
    "Island",
    "Short Name"
  )

#Create a matrix For a site Comparision of rainfall
Cell.RF_Year <- data.frame(matrix(ncol = 15, nrow = UNIT_C))
colnames(Cell.RF_Year) <-
  c(
    "Unit",
    "Mean ELEV",
    "JAN",
    "FEB",
    "MAR",
    "APR",
    "MAY",
    "JUN",
    "JUL",
    "AUG",
    "SEP",
    "OCT",
    "NOV",
    "DEC",
    "ANN"
  )
Cell.RF_Year[1] <- UNIT_N

u <- 1
# for (u in 1:1) {

print("Analysis UNIT ")  # Unit Analyzed
print(u)                 # Number In Loop
print(UNIT_N[u])         # Name Of Unit
print("Maps")

#Create A Directory for Output
path <- paste0(OUTPUTS_FOLDER, UNIT_N[u])
print(paste("PDKE 29, PATH:", path))
dir.create(path, showWarnings = TRUE, recursive = FALSE)
#
# TODO: the above got the following, so add error handling:
# Warning message:
# In dir.create(path, showWarnings = TRUE, recursive = FALSE) :
#  '/Users/jgeis/Work/PDKE/CCVD/MINI_Phase2/Waikiki Watershed' already exists


IS <- UNIT_I[u]       # Island
SH <- UNIT_X[[u]]     # Shape File

##########   Get Coast for Island
print("PDKE: 30, Get Coast for Island")

if (IS == "BI") {
  CoastM <- Coast_BI
  Iname <- "Hawaii"
}
if (IS == "MN") {
  CoastM <- Coast_MN
  Iname <- "Maui"
}
if (IS == "MO") {
  CoastM <- Coast_MO
  Iname <- "Molokai"
}
if (IS == "LA") {
  CoastM <- Coast_LA
  Iname <- "Lanai"
}
if (IS == "OA") {
  CoastM <- Coast_OA
  Iname <- "Oahu"
}
if (IS == "KA") {
  CoastM <- Coast_KA
  Iname <- "Kauai"
}
if (IS == "KO") {
  CoastM <- Coast_KO
  Iname <- "Kahoolawe"
}
if (IS == "ALL") {
  CoastM <- Coast_Crop_T
  Iname <- "Hawaiian Islands"
}

##########   GET Central Latitude and Longitude for shapefile

Xmin <- extent(SH)[1]
Xmax <- extent(SH)[2]
ymin <- extent(SH)[3]
ymax <- extent(SH)[4]

LONG <- round(mean(Xmin, Xmax), 3)
LAT <- round(mean(ymin, ymax), 3)


##########   Make a Plot of Shapefile with Island Coast and LAT and LON
print("PDKE: 31, Make a Plot of Shapefile with Island Coast and LAT and LON")

TIT <- paste0(UNIT_Ns[u], " (", Iname, ")")  # Title For Figure
dpi <- 300

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

if (sb <= 3) {
  sb2 <- 1
}
if (sb > 3 && sb <= 5) {
  sb2 <- 5
}
if (sb > 5 && sb <= 10) {
  sb2 <- 10
}
if (sb > 10 && sb <= 20) {
  sb2 <- 20
}

sb2 <- sb2 * 1.60934

############# plot Mokupuni (island)
print("PDKE: 32, Mokupuni")

dpi <- 300
png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Mokupuni.png"),
  width = 5 * dpi,
  height = 4 * dpi,
  res = dpi
)

ggplot() +
  geom_polygon(
    data = Coast,
    aes(x = long, y = lat, group = group),
    col = "black",
    fill = "aliceblue",
    size = 0.5
  ) +
  geom_polygon(
    data = CoastM,
    aes(x = long, y = lat, group = group),
    col = "brown",
    fill = NA,
    size = 0.5
  ) +
  geom_polygon(
    data = HALE,
    aes(x = long, y = lat, group = group),
    col = "blue",
    fill = NA,
    size = 0.25
  ) +
  coord_fixed() +
  theme_void()

dev.off()

# par(mfrow=c(1,1))
# plot(Coast_Crop_T, lwd=1,lty = "solid",axes=F,las=2)
# plot(HALE,add=T,col='blue')
# plot(CoastM, add=T, border="brown")

############# plot Moku
print("PDKE: 33, Plot Moku")


# find which moku the study area is in
# MI<-as.data.frame(gIntersects(MOKU,HALE,byid=T))
# MI
# MI2<-as.data.frame(colnames(MI[which(MI == TRUE)]))
# colnames(MI2)<-c("objectID")

MIt <- as.data.frame(raster::intersect(MOKU, HALE))
MI2 <- data.frame(MIt$objectid)
MI2
colnames(MI2) <- "ObjectID"


# # add 1 to the objectID's because they're off by 1...
# MI2$objectID<-as.numeric(MI2$objectID) + 1
# MI2
#
# # iterate through moku objectID's to make a list of them
# r2<-c()
# for (i in 1:nrow(MI2)) {
#
#   m<-MI2$objectID[i]
#   r<-c(m)
#   r2<-c(r2,r)
# }
# r2

# subset moku for study area
MI3 <- MOKU[MI2$ObjectID, ]
plot(MI3)

# get subsetted moku as a dataframe for text labels
MId <- as(MI3, "data.frame")
write.csv(
  data.frame(MId),
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Moku.csv"),
  row.names = F
)

# # get basemap at moku extent
# mapBound <- c(xmin(MI3), ymin(MI3), xmax(MI3), ymax(MI3))
# basemap<-ggmap::get_stamenmap(bbox=mapBound, maptype = "terrain")
# dpi=300

# get extent of Moku
xm <- xmin(MI3)
xma <- xmax(MI3)
ym <- ymin(MI3) * 0.999
yma <- ymax(MI3) * 1.0001

# make map
png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Moku.png"),
  width = 5 * dpi,
  height = 4 * dpi,
  res = dpi
)

# par(mfrow=c(1,1))
# plot(Coast_Crop_T, lwd=1,lty = "solid",axes=F,las=2)
# plot(HALE,add=T,col='blue')
# plot(CoastM, add=T, border="brown")
# par(mfrow=c(1,1))

ggplot() +
  # geom_tile(data = dplyr::filter(HSd, !is.na(value)),
  #           aes(x=x, y=y, fill=value), show.legend = F) +
  # scale_fill_gradient("Elevation", low = "black", high="white", na.value = 0) +
  geom_polygon(
    data = Coast,
    aes(x = long, y = lat, group = group),
    col = "black",
    fill = "aliceblue",
    alpha = 50,
    size = 1
  ) +
  geom_polygon(
    data = MI3,
    aes(x = long, y = lat, group = group),
    col = "brown",
    fill = NA,
    size = 1
  ) +
  geom_polygon(
    data = HALE,
    aes(x = long, y = lat, group = group),
    col = "blue",
    fill = NA,
    size = 1
  ) +
  geom_text(
    data = MId,
    aes(label = moku, x = cent_long, y = cent_lat),
    colour = "brown",
    size = 4,
    fontface = "bold"
  ) +
  coord_sf(xlim = c(xm, xma), ylim = c(ym, yma)) +
  theme_void()

dev.off()

########## Plot Ahupuaa

print("PDKE: 34, Plot Ahupuaa")

# Clip ahupuaa polygons to HALE polygons (removes all attribute data, use objectid
# to re-assign attribute data afterwards)
# AI<-gIntersection(AHU, HALE, byid=T, id = as.character(AHU$objectid))
# plot(AI)
plot(HALE)
plot(AHU)
AIt <- as.data.frame(raster::intersect(AHU, HALE))
AI2 <- data.frame(AIt$objectid)
AI2
colnames(AI2) <- "ObjectID"

# subset ahupuaa for study area
AI3 <- AHU[AI2$ObjectID, ]
plot(AI3)

# get subsetted ahupuaa as a dataframe for text labels
AId <- as(AI3, "data.frame")
write.csv(
  data.frame(AId),
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Ahupuaa.csv"),
  row.names = F
)

# get extent of Ahupuaa
xm <- xmin(AI3)
xma <- xmax(AI3)
ym <- ymin(AI3) * 0.999
yma <- ymax(AI3) * 1.0001

dpi <- 300

# make map
png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Ahupuaa.png"),
  width = 5 * dpi,
  height = 4 * dpi,
  res = dpi
)

ggplot() +
  geom_polygon(
    data = Coast,
    aes(x = long, y = lat, group = group),
    col = "black",
    fill = "aliceblue",
    alpha = 50,
    size = 1
  ) +
  geom_polygon(
    data = AI3,
    aes(x = long, y = lat, group = group),
    col = "brown",
    fill = NA,
    size = 1
  ) +
  geom_polygon(
    data = HALE,
    aes(x = long, y = lat, group = group),
    col = "blue",
    fill = NA,
    size = 1
  ) +
  geom_text(
    data = AId,
    aes(label = ahupuaa2, x = cent_long, y = cent_lat),
    size = 4,
    fontface = "bold"
  ) +
  coord_sf(xlim = c(xm, xma), ylim = c(ym, yma)) +
  theme_void()

dev.off()

##########  Fire Maps

print("PDKE: 35, Fire Maps")

# Fire history map
FIRE_CropI <- FIRE_Shape_T
TITF <- "Fire Occurrence 1999-2022"

dpi <- 300
png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Fire.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

par(mfrow = c(1, 1))
plot(
  CoastM,
  lwd = 1,
  lty = "solid",
  axes = F,
  las = 2
)
title(TITF , line = 0.5, cex.main = 1.5)
plot(SH,
     add = TRUE,
     lwd = 2,
     col = "cyan")
plot(FIRE_CropI, add = TRUE, col = "darkred")


legend(
  "topright",
  legend = c("Fire Occurrence", UNIT_Ns[[u]]),
  pch = c(15, 15),
  col = c("darkred", "cyan"),
  bty = "n",
  horiz = FALSE,
  inset = c(0.05, 0.05),
  cex = 0.8
)

dev.off()

# Fire risk map

# get extent of AOI
xm <- xmin(HALE)
xma <- xmax(HALE)
ym <- ymin(HALE) * 0.9999
yma <- ymax(HALE) * 1.0001

# separate fire risk categories
FR1 <- FRISK_Shape_T[which(FRISK_Shape_T$gridcode == "101"), ]
FR2 <- FRISK_Shape_T[which(FRISK_Shape_T$gridcode == "102"), ]
FR3 <- FRISK_Shape_T[which(FRISK_Shape_T$gridcode == "103"), ]

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " FireRisk.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

par(mfrow = c(1, 1))
ggplot() +
  geom_polygon(data = FR1,
               aes(x = long, y = lat, group = group),
               fill = "yellow") +
  geom_polygon(data = FR2,
               aes(x = long, y = lat, group = group),
               fill = "orange") +
  geom_polygon(data = FR3,
               aes(x = long, y = lat, group = group),
               fill = "red") +
  geom_polygon(
    data = Coast,
    aes(x = long, y = lat, group = group),
    col = "black",
    fill = NA,
    size = 1
  ) +
  geom_polygon(
    data = HALE,
    aes(x = long, y = lat, group = group),
    col = "blue",
    fill = NA,
    size = 1
  ) +
  coord_sf(xlim = c(xm, xma), ylim = c(ym, yma)) +
  theme_void()

dev.off()

########### Elevation Extract and Map

print("PDKE: 36, Elevation")

ELEV_MaskI <- mask(x = ELEV, mask = CoastM)
ELEV_CropI <- crop(x = ELEV_MaskI, y = extent(CoastM))

#Mask For Geography
ELEV_Mask <- mask(x = ELEV, mask = UNIT_X[[u]])
ELEV_Crop <- crop(x = ELEV_Mask, y = extent(UNIT_X[[u]]))
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

SHAPE <- UNIT_X[[u]]

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " ELMap.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

print(
  spplot(
    ELEV_CropI,
    col.regions = colfuncEL,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    main = list(
      label = paste0("Elevation", " ", Iname, " (", ELUnit2, ")"),
      cex = 0.8
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    )
  ) +
    layer(sp.polygons(SHAPE, lwd = 2)) +
    layer(sp.polygons(CoastM, lwd = 1))
)

dev.off()

########## Landcover Extraction and Map
print("PDKE: 37, Landcover")

LC_MaskI <- mask(x = LC, mask = CoastM)
LC_CropI <- crop(x = LC_MaskI, y = extent(CoastM))
LC_MaskI
plot(LC)
plot(CoastM)

#Mask For Geography
LC_Mask <- mask(x = LC, mask = UNIT_X[[u]])
LC_Crop <- crop(x = LC_Mask, y = extent(UNIT_X[[u]]))

plot(LC_Crop)

#Convert raster to dataframe for stats
LC_Crop2 <- as.data.frame(LC_Crop, xy = T)

# remove LC "0" (no data)
LC_Crop3 <-
  LC_Crop2[which(LC_Crop2$LCMAP_HI_2020_V10_LCPRI_crs != 0), ]

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
  if (LC_ct2[x, ]$pct == 0) {
    LC_ct2[x, ]$pct <- "<1"
  }
}

LC_ct2$pct <- paste0(LC_ct2$pct, "%")
LC_ct2

#define class names and colors
LC_ct2$class_name <- as.character(NA)
LC_ct2$color <- as.character(NA)

for (x in 1:nrow(LC_ct2)) {
  if (LC_ct2[x, ]$LC == 1) {
    LC_ct2[x, ]$class_name <- c("Developed")
  }
  if (LC_ct2[x, ]$LC == 2) {
    LC_ct2[x, ]$class_name <- c("Cropland")
  }
  if (LC_ct2[x, ]$LC == 3) {
    LC_ct2[x, ]$class_name <- c("Grass/Shrub")
  }
  if (LC_ct2[x, ]$LC == 4) {
    LC_ct2[x, ]$class_name <- c("Tree Cover")
  }
  if (LC_ct2[x, ]$LC == 5) {
    LC_ct2[x, ]$class_name <- c("Water")
  }
  if (LC_ct2[x, ]$LC == 6) {
    LC_ct2[x, ]$class_name <- c("Wetland")
  }
  if (LC_ct2[x, ]$LC == 8) {
    LC_ct2[x, ]$class_name <- c("Barren")
  }
  
}

LC_ct2

# write to table
write.csv(LC_ct2,
          paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Landcover.csv"),
          row.names = F)

# set class names as leveled factors
LC_ct2$class_name <- factor(
  LC_ct2$class_name,
  levels = c(
    "Developed",
    "Cropland",
    "Grass/Shrub",
    "Tree Cover",
    "Water",
    "Wetland",
    "Barren"
  )
)

# set first thru third landcover classes
LC1 <- paste0(LC_ct2[1, ]$class_name, "(", LC_ct2[1, ]$pct, ")")
LC2 <- paste0(LC_ct2[2, ]$class_name, "(", LC_ct2[2, ]$pct, ")")
LC3 <- paste0(LC_ct2[3, ]$class_name, "(", LC_ct2[3, ]$pct, ")")

### Plot Bar Graph
# set bargraph ylim
ylim <- max(LC_ct2$acres) + (max(LC_ct2$acres) * 0.15)

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " LC_barchart.png"),
  width = 5 * dpi,
  height = 3 * dpi,
  res = dpi
)

print(
  ggplot(LC_ct2, aes(
    x = class_name, y = acres, fill = class_name
  )) +
    geom_col(
      width = 0.9,
      position = position_dodge(0.7),
      color = "black"
    ) +
    theme_bw() +
    ylim(0, ylim) +
    scale_fill_manual(
      values = c(
        "Developed" = "black",
        "Cropland" = "darkgoldenrod1",
        "Grass/Shrub" = "darkolivegreen1",
        "Tree Cover" = "darkgreen",
        "Water" = "blue",
        "Wetland" = "aquamarine",
        "Barren" = "chocolate4"
      )
    ) +
    geom_text(aes(label = pct), vjust = -0.8, size = 4.5) +
    theme(axis.text.x = element_text(
      angle = 35,
      vjust = .65,
      size = 14
    )) +
    labs(y = "Area (acres)", x = "",) +
    theme(legend.position = "none") +
    theme(
      axis.title.y = element_text(size = 15),
      axis.text.y = element_text(size = 14)
    )
)

dev.off()

#
# Cell.DataCL[u,14:16] <- c(LC1,LC2,LC3)
# Cell.DataCLR[1,12:14] <- c(LC1,LC2,LC3)

SHAPE <- UNIT_X[[u]]

### remove "0" values from raster
LC_CropI2 <- reclassify(LC_CropI, cbind(0, NA))
summary(LC_CropI2)

# export map
TITF <- paste0("Landcover", " ", Iname)

dpi = 300

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " LCMap.png"),
  width = 7 * dpi,
  height = 5 * dpi,
  res = dpi
)

par(mar = c(1, 1, 2, 0.5))
plot(
  LC_CropI2,
  legend = F,
  col = c(
    "black",
    "darkgoldenrod1",
    "darkolivegreen1",
    "darkgreen",
    "aquamarine",
    "blue",
    "chocolate4"
  ),
  xaxt = 'n',
  yaxt = 'n'
)
plot(CoastM, add = T)
plot(
  SHAPE,
  lwd = 1,
  border = "red",
  lwd = 2,
  add = T
)
title(TITF , line = 0.5, cex.main = 1.5)

dev.off()


################################################################################
########## Water Sources
print("PDKE: 38, Water Sources")

### Aquifers
# AQ1<-as.data.frame(gIntersects(AQU, HALE, byid=T))
# AQ2<-as.data.frame(colnames(AQ1[which(AQ1 == TRUE)]))
# colnames(AQ2)<-c("objectID")
# head(AQ2)

AIt <- as.data.frame(raster::intersect(AQU, HALE))
AI2 <- data.frame(AIt$objectid)
colnames(AI2) <- "ObjectID"

# subset aquifers for study area
AQ3 <- AQU[AI2$ObjectID, ]
plot(AQ3)


AQ3 <- raster::intersect(AQU, HALE)

# subset for hydrology
AQb <- AQ3[which(AQ3$typea1 == 1), ]
AQh <- AQ3[which(AQ3$typea1 == 2), ]

# make data frame and only keep relevant columns
AQd <- as(AQ3, "data.frame")
AQd <-
  subset(
    AQd,
    select = c(
      "objectid",
      "doh_aquife",
      "typea1",
      "typea3",
      "stata1",
      "strata3",
      "cent_lat",
      "cent_long"
    )
  )

# translate type and status codes to text
AQd$Hydrology <- NA
AQd$Geology <- NA
AQd$Salinity <- NA
AQd$Use <- NA


for (i in 1:nrow(AQd)) {
  if (AQd[i, ]$typea1 == 1) {
    AQd[i, ]$Hydrology <- c("Basal")
  }
  if (AQd[i, ]$typea1 == 2) {
    AQd[i, ]$Hydrology <- c("High Level")
  }
  if (AQd[i, ]$typea3 == 1) {
    AQd[i, ]$Geology <- c("Flank")
  }
  if (AQd[i, ]$typea3 == 2) {
    AQd[i, ]$Geology <- c("Dike")
  }
  if (AQd[i, ]$typea3 == 3) {
    AQd[i, ]$Geology <- c("Flank/Dike")
  }
  if (AQd[i, ]$typea3 == 4) {
    AQd[i, ]$Geology <- c("Perched")
  }
  if (AQd[i, ]$typea3 == 5) {
    AQd[i, ]$Geology <- c("Dike/Perched")
  }
  if (AQd[i, ]$typea3 == 6) {
    AQd[i, ]$Geology <- c("Sedimentary")
  }
  if (AQd[i, ]$stata1 == 1) {
    AQd[i, ]$Use <- c("Currently used")
  }
  if (AQd[i, ]$stata1 == 2) {
    AQd[i, ]$Use <- c("Potential use")
  }
  if (AQd[i, ]$stata1 == 3) {
    AQd[i, ]$Use <- c("No potential use")
  }
  if (AQd[i, ]$strata3 == 1) {
    AQd[i, ]$Salinity <- c("Fresh")
  }
  if (AQd[i, ]$strata3 == 2) {
    AQd[i, ]$Salinity <- c("Low")
  }
  if (AQd[i, ]$strata3 == 3) {
    AQd[i, ]$Salinity <- c("Moderate")
  }
  if (AQd[i, ]$strata3 == 4) {
    AQd[i, ]$Salinity <- c("High")
  }
  if (AQd[i, ]$strata3 == 5) {
    AQd[i, ]$Salinity <- c("Seawater")
  }
  
}

AQd

# trim it down to the final table
AQe <-
  subset(AQd,
         select = c("doh_aquife", "Hydrology", "Geology", "Salinity", "Use"))
colnames(AQe)[1] = "DOH Aquifer"
write.csv(
  data.frame(AQe),
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Aquifer.csv"),
  row.names = F
)
AQe


# get extent of Aquifrs
xm <- xmin(AQ3)
xma <- xmax(AQ3)
ym <- ymin(AQ3) * 0.999
yma <- ymax(AQ3) * 1.0001

# make map
dpi = 300

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Aquifers.png"),
  width = 4 * dpi,
  height = 4 * dpi,
  res = dpi
)
par(mfrow = c(1, 1))

# has conditionals based on which aquifer types are present
if (nrow(data.frame(AQb)) > 0 && nrow(data.frame(AQh)) > 0) {
  # print(ggmap(basemap, extent="device") +
  ggplot() +
    geom_polygon(
      data = CoastM,
      aes(x = long, y = lat, group = group),
      fill = "gray90",
      col = "black"
    ) +
    geom_polygon(
      data = AQb,
      aes(x = long, y = lat, group = group),
      fill = "lightblue",
      col = "black"
    ) +
    geom_polygon(
      data = AQh,
      aes(x = long, y = lat, group = group),
      fill = "grey",
      col = "black"
    ) +
    geom_polygon(
      data = HALE,
      aes(x = long, y = lat, group = group),
      col = "red",
      fill = NA,
      size = 1
    ) +
    geom_text(
      data = AQd,
      aes(label = doh_aquife, x = cent_long, y = cent_lat),
      size = 3,
      fontface = "bold"
    ) +
    theme_void() +
    coord_sf(xlim = c(xm, xma), ylim = c(ym, yma))
}

if (nrow(data.frame(AQb)) == 0 && nrow(data.frame(AQh)) > 0) {
  # print(ggmap(basemap, extent="device") +
  ggplot() +
    geom_polygon(
      data = CoastM,
      aes(x = long, y = lat, group = group),
      fill = "gray90",
      col = "black"
    ) +
    geom_polygon(
      data = AQh,
      aes(x = long, y = lat, group = group),
      fill = "grey",
      col = "black"
    ) +
    geom_polygon(
      data = HALE,
      aes(x = long, y = lat, group = group),
      col = "red",
      fill = NA,
      size = 1
    ) +
    geom_text(
      data = AQd,
      aes(label = doh_aquife, x = cent_long, y = cent_lat),
      size = 3,
      fontface = "bold"
    ) +
    theme_void() +
    coord_sf(xlim = c(xm, xma), ylim = c(ym, yma))
  
}

if (nrow(data.frame(AQb)) > 0 && nrow(data.frame(AQh)) == 0) {
  # print(ggmap(basemap, extent="device") +
  ggplot() +
    geom_polygon(
      data = CoastM,
      aes(x = long, y = lat, group = group),
      fill = "gray90",
      col = "black"
    ) +
    geom_polygon(
      data = AQb,
      aes(x = long, y = lat, group = group),
      fill = "lightblue",
      col = "black"
    ) +
    geom_polygon(
      data = HALE,
      aes(x = long, y = lat, group = group),
      col = "red",
      fill = NA,
      size = 1
    ) +
    geom_text(
      data = AQd,
      aes(label = doh_aquife, x = cent_long, y = cent_lat),
      size = 3,
      fontface = "bold"
    ) +
    theme_void() +
    coord_sf(xlim = c(xm, xma), ylim = c(ym, yma))
  
}

dev.off()

################################################################################
### Streams, hydrologic features

# # get basemap at buffered AOI extent
# mapBound1 <- c((xmin(HALE)*1.00001), (ymin(HALE)*0.9999), (xmax(HALE)*0.99999), (ymax(HALE)*1.0001))
# basemap<-ggmap::get_stamenmap(bbox=mapBound1, maptype = "terrain")
print("PDKE: 39, Streams")

# subset for feature types
STRMi <- STRM[which(STRM$fcode == 46003), ]
STRMp <- STRM[which(STRM$fcode == 46006), ]
STRMa <- STRM[which(STRM$fcode == 55800), ]
STRMc <- STRM[which(STRM$Feature_Ty == "CANAL/DITCH"), ]
STRMpi <- STRM[which(STRM$Feature_Ty == "PIPELINE"), ]
STRMco <- STRM[which(STRM$Feature_Ty == "CONNECTOR"), ]

### prepare hillshade for plotting
#clip to island
HSmoku <- mask(x = HS, mask = CoastM)

HSmoku <- crop(x = HS, y = extent(CoastM))
HSmoku[HSmoku == 0] <- NA

# use this function I found online (https://stackoverflow.com/questions/47116217/overlay-raster-layer-on-map-in-ggplot2-in-r)
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

# get extent of AOI again
xm <- xmin(HALE)
xma <- xmax(HALE)
ym <- ymin(HALE) * 0.9999
yma <- ymax(HALE) * 1.0001

# make map
dpi = 300
png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Streams.png"),
  width = 4 * dpi,
  height = 4 * dpi,
  res = dpi
)

par(mfrow = c(1, 1))
ggplot() +
  geom_tile(
    data = dplyr::filter(HSd,!is.na(value)),
    aes(x = x, y = y, fill = value),
    show.legend = F
  ) +
  scale_fill_gradient(
    "Elevation",
    low = "black",
    high = "white",
    na.value = 0
  ) +
  geom_polygon(
    data = Coast,
    aes(x = long, y = lat, group = group),
    fill = NA,
    col = "black",
    size = 0.5
  ) +
  geom_path(
    data = STRMi,
    aes(x = long, y = lat, group = group),
    col = "deepskyblue1",
    size = 1
  ) +
  geom_path(
    data = STRMp,
    aes(x = long, y = lat, group = group),
    col = "darkblue",
    size = 1
  ) +
  geom_path(
    data = STRMa,
    aes(x = long, y = lat, group = group),
    col = "green",
    size = 1
  ) +
  geom_path(
    data = STRMc,
    aes(x = long, y = lat, group = group),
    col = "orange",
    size = 1
  ) +
  geom_path(
    data = STRMpi,
    aes(x = long, y = lat, group = group),
    col = "lightgrey",
    size = 1
  ) +
  geom_path(
    data = STRMco,
    aes(x = long, y = lat, group = group),
    col = "yellow",
    size = 1
  ) +
  geom_polygon(
    data = HALE,
    aes(x = long, y = lat, group = group),
    col = "red",
    fill = NA,
    size = 1
  ) +
  coord_sf(xlim = c(xm, xma), ylim = c(ym, yma)) +
  theme_void()

dev.off()

# clip features to study area
ST <- crop(STRM, HALE)
plot(ST)

# make dataframe of features
STdf <- as(ST, "data.frame")

# get all hydrologic types present
ht <- unique(STdf$Feature_Ty)
for (i in 1:length(ht)) {
  if (ht[i] == "ARTIFICIAL PATH") {
    ht[i] <- "Managed Waterway"
  }
  if (ht[i] == "CANAL/DITCH") {
    ht[i] <- "Canal/Ditch"
  }
  if (ht[i] == "STREAM/RIVER") {
    ht[i] <- "Stream"
  }
  if (ht[i] == "PIPELINE") {
    ht[i] <- "Pipeline"
  }
}

write.csv(
  data.frame(ht),
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Hydro_features.csv"),
  row.names = F
)


########## Rainfall Stations near AOI
print("PDKE: 40, Rain Station Locations")

# load csv of station networks and website links
sn <- read.csv(paste0(INPUTS_FOLDER, "rain_stations_links.csv"))

# load stations point shapefile, tranform points and polygons to planar projection
sl <- readOGR(paste0(INPUTS_FOLDER, "rain_stations_2023.shp"))
# sl <- spTransform(sl, crs(EXAMP))

utmStr <-
  "+proj=utm +zone=%d +datum=NAD83 +units=m +no_defs +ellps=GRS80"
crs <- CRS(sprintf(utmStr, 4))
sl <- spTransform(sl, crs)
HALE1 <- spTransform(HALE, crs)

summary(sl)
unique(sl$Network)

# calculate distance from each point to the AOI polygon
sd <- as.data.frame(apply(gDistance(sl, HALE1, byid = TRUE), 2, min))
sd <- cbind(rownames(sd), data.frame(sd, row.names = NULL))
colnames(sd) <- c("ID2", "distance")

# subtract 1 from each FID because they're all off by one (started with 0)
sd$ID2 <- as.numeric(sd$ID2) - 1
sd <- sd[order(sd$distance), ]
head(sd)

# convert station points shp into dataframe
summary(sl)
sld <- as.data.frame(sl)
head(sld)

# make table of 3 closest stations and join network links
s1 <- sld[which(sld$ID2 == sd[1, ]$ID2), ]
s2 <- sld[which(sld$ID2 == sd[2, ]$ID2), ]
s3 <- sld[which(sld$ID2 == sd[3, ]$ID2), ]
st2 <- rbind(s1, s2, s3)
st2
st <- st2[c("Station_Na", "Network")]
st <- join(st, sn)
colnames(st) <- c("Station Name", "Network", "Website")
st

# export table
write.csv(st,
          paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " rain_stations.csv"),
          row.names = F)

# ### Make map of island with station and AOI
# # get basemap at island extent
# mapBound <- c(xmin(CoastM), ymin(CoastM), xmax(CoastM), ymax(CoastM))
# basemap<-ggmap::get_stamenmap(bbox=mapBound, maptype = "terrain-background")

# get extent of island
xm <- xmin(CoastM)
xma <- xmax(CoastM)
ym <- ymin(CoastM)
yma <- ymax(CoastM)

# make map
png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " rf_stations.png"),
  width = 5 * dpi,
  height = 3.5 * dpi,
  res = dpi
)

par(mfrow = c(1, 1))

ggplot() +
  geom_tile(
    data = dplyr::filter(HSd,!is.na(value)),
    aes(x = x, y = y, fill = value),
    show.legend = F
  ) +
  scale_fill_gradient(
    "Elevation",
    low = "black",
    high = "white",
    na.value = 0
  ) +
  geom_polygon(
    data = HALE,
    aes(x = long, y = lat, group = group),
    col = "blue",
    fill = NA,
    size = 1
  ) +
  geom_point(
    data = sld,
    aes(x = LON, y = LAT),
    col = "black",
    size = 1.5
  ) +
  geom_point(
    data = st2,
    aes(x = LON, y = LAT),
    color = "black",
    fill = "orange",
    shape = 21,
    size = 2
  ) +
  coord_sf(xlim = c(xm, xma), ylim = c(ym, yma)) +
  theme_void()

dev.off()

########## Mean Climate Extract
print("PDKE: 41, Mean Climate")

#Create a climate matrix
Cell.matrix5 <- matrix(ncol = 14, nrow = 9)
Cell.CL_Year <- data.frame(Cell.matrix5)
colnames(Cell.CL_Year) <-
  c(
    "Variable",
    "JAN",
    "FEB",
    "MAR",
    "APR",
    "MAY",
    "JUN",
    "JUL",
    "AUG",
    "SEP",
    "OCT",
    "NOV",
    "DEC",
    "ANN"
  )
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
print("PDKE: 42, Solar Radiation")

Jan <- raster(Mean_CLIM[58])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropKD <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
JanMKD   <- round(cellStats(Jan_CropKD, 'mean'), 0)
JanMKDx   <- round(cellStats(Jan_CropKD, 'max'), 0)
JanMKDn   <- round(cellStats(Jan_CropKD, 'min'), 0)

Feb <- raster(Mean_CLIM[57])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropKD <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
FebMKD   <- round(cellStats(Feb_CropKD, 'mean'), 0)
FebMKDx   <- round(cellStats(Feb_CropKD, 'max'), 0)
FebMKDn   <- round(cellStats(Feb_CropKD, 'min'), 0)

Mar <- raster(Mean_CLIM[61])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropKD <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
MarMKD   <- round(cellStats(Mar_CropKD, 'mean'), 0)
MarMKDx   <- round(cellStats(Mar_CropKD, 'max'), 0)
MarMKDn   <- round(cellStats(Mar_CropKD, 'min'), 0)

Apr <- raster(Mean_CLIM[54])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropKD <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
AprMKD   <- round(cellStats(Apr_CropKD, 'mean'), 0)
AprMKDx   <- round(cellStats(Apr_CropKD, 'max'), 0)
AprMKDn   <- round(cellStats(Apr_CropKD, 'min'), 0)

May <- raster(Mean_CLIM[62])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropKD <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
MayMKD   <- round(cellStats(May_CropKD, 'mean'), 0)
MayMKDx   <- round(cellStats(May_CropKD, 'max'), 0)
MayMKDn   <- round(cellStats(May_CropKD, 'min'), 0)

Jun <- raster(Mean_CLIM[60])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropKD <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
JunMKD   <- round(cellStats(Jun_CropKD, 'mean'), 0)
JunMKDx   <- round(cellStats(Jun_CropKD, 'max'), 0)
JunMKDn   <- round(cellStats(Jun_CropKD, 'min'), 0)

Jul <- raster(Mean_CLIM[59])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropKD <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
JulMKD   <- round(cellStats(Jul_CropKD, 'mean'), 0)
JulMKDx   <- round(cellStats(Jul_CropKD, 'max'), 0)
JulMKDn   <- round(cellStats(Jul_CropKD, 'min'), 0)

Aug <- raster(Mean_CLIM[55])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropKD <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
AugMKD   <- round(cellStats(Aug_CropKD, 'mean'), 0)
AugMKDx   <- round(cellStats(Aug_CropKD, 'max'), 0)
AugMKDn   <- round(cellStats(Aug_CropKD, 'min'), 0)

Sep <- raster(Mean_CLIM[65])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropKD <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
SepMKD   <- round(cellStats(Sep_CropKD, 'mean'), 0)
SepMKDx   <- round(cellStats(Sep_CropKD, 'max'), 0)
SepMKDn   <- round(cellStats(Sep_CropKD, 'min'), 0)

Oct <- raster(Mean_CLIM[64])
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropKD <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
OctMKD   <- round(cellStats(Oct_CropKD, 'mean'), 0)
OctMKDx   <- round(cellStats(Oct_CropKD, 'max'), 0)
OctMKDn   <- round(cellStats(Oct_CropKD, 'min'), 0)

Nov <- raster(Mean_CLIM[63])
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropKD <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
NovMKD   <- round(cellStats(Nov_CropKD, 'mean'), 0)
NovMKDx   <- round(cellStats(Nov_CropKD, 'max'), 0)
NovMKDn   <- round(cellStats(Nov_CropKD, 'min'), 0)

Dec <- raster(Mean_CLIM[56])
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropKD <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
DecMKD   <- round(cellStats(Dec_CropKD, 'mean'), 0)
DecMKDx   <- round(cellStats(Dec_CropKD, 'max'), 0)
DecMKDn   <- round(cellStats(Dec_CropKD, 'min'), 0)

Ann <- raster(Mean_CLIM[53])
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropKD <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
AnnMKD   <- round(cellStats(Ann_CropKD, 'mean'), 0)
AnnMKDx   <- round(cellStats(Ann_CropKD, 'max'), 0)
AnnMKDn   <- round(cellStats(Ann_CropKD, 'min'), 0)

MEANKD2 <-
  c(
    JanMKD,
    FebMKD,
    MarMKD,
    AprMKD,
    MayMKD,
    JunMKD,
    JulMKD,
    AugMKD,
    SepMKD,
    OctMKD,
    NovMKD,
    DecMKD,
    AnnMKD
  )
MEANKD3 <- MEANKD2 * 100

#GKD thresholds for scale Bar
KDUP <-
  max(
    JanMKDx,
    FebMKDx,
    MarMKDx,
    AprMKDx,
    MayMKDx,
    JunMKDx,
    JulMKDx,
    AugMKDx,
    SepMKDx,
    OctMKDx,
    NovMKDx,
    DecMKDx,
    AnnMKDx
  )
KDLO <-
  min(
    JanMKDn,
    FebMKDn,
    MarMKDn,
    AprMKDn,
    MayMKDn,
    JunMKDn,
    JulMKDn,
    AugMKDn,
    SepMKDn,
    OctMKDn,
    NovMKDn,
    DecMKDn,
    AnnMKDn
  )

Cell.CL_Year[9, 2:14] <- MEANKD2


##########   Soil Moisture
print("PDKE: 43, Soil Moisture")

Jan <- raster(Mean_CLIM[45])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropSM <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
JanMSM   <- round(cellStats(Jan_CropSM, 'mean'), 2)
JanMSMx   <- round(cellStats(Jan_CropSM, 'max'), 2)
JanMSMn   <- round(cellStats(Jan_CropSM, 'min'), 2)

Feb <- raster(Mean_CLIM[44])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropSM <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
FebMSM   <- round(cellStats(Feb_CropSM, 'mean'), 2)
FebMSMx   <- round(cellStats(Feb_CropSM, 'max'), 2)
FebMSMn   <- round(cellStats(Feb_CropSM, 'min'), 2)

Mar <- raster(Mean_CLIM[48])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropSM <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
MarMSM   <- round(cellStats(Mar_CropSM, 'mean'), 2)
MarMSMx   <- round(cellStats(Mar_CropSM, 'max'), 2)
MarMSMn   <- round(cellStats(Mar_CropSM, 'min'), 2)

Apr <- raster(Mean_CLIM[41])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropSM <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
AprMSM   <- round(cellStats(Apr_CropSM, 'mean'), 2)
AprMSMx   <- round(cellStats(Apr_CropSM, 'max'), 2)
AprMSMn   <- round(cellStats(Apr_CropSM, 'min'), 2)

May <- raster(Mean_CLIM[49])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropSM <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
MayMSM   <- round(cellStats(May_CropSM, 'mean'), 2)
MayMSMx   <- round(cellStats(May_CropSM, 'max'), 2)
MayMSMn   <- round(cellStats(May_CropSM, 'min'), 2)

Jun <- raster(Mean_CLIM[47])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropSM <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
JunMSM   <- round(cellStats(Jun_CropSM, 'mean'), 2)
JunMSMx   <- round(cellStats(Jun_CropSM, 'max'), 2)
JunMSMn   <- round(cellStats(Jun_CropSM, 'min'), 2)

Jul <- raster(Mean_CLIM[46])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropSM <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
JulMSM   <- round(cellStats(Jul_CropSM, 'mean'), 2)
JulMSMx   <- round(cellStats(Jul_CropSM, 'max'), 2)
JulMSMn   <- round(cellStats(Jul_CropSM, 'min'), 2)

Aug <- raster(Mean_CLIM[42])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropSM <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
AugMSM   <- round(cellStats(Aug_CropSM, 'mean'), 2)
AugMSMx   <- round(cellStats(Aug_CropSM, 'max'), 2)
AugMSMn   <- round(cellStats(Aug_CropSM, 'min'), 2)

Sep <- raster(Mean_CLIM[52])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropSM <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
SepMSM   <- round(cellStats(Sep_CropSM, 'mean'), 2)
SepMSMx   <- round(cellStats(Sep_CropSM, 'max'), 2)
SepMSMn   <- round(cellStats(Sep_CropSM, 'min'), 2)

Oct <- raster(Mean_CLIM[51])
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropSM <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
OctMSM   <- round(cellStats(Oct_CropSM, 'mean'), 2)
OctMSMx   <- round(cellStats(Oct_CropSM, 'max'), 2)
OctMSMn   <- round(cellStats(Oct_CropSM, 'min'), 2)

Nov <- raster(Mean_CLIM[50])
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropSM <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
NovMSM   <- round(cellStats(Nov_CropSM, 'mean'), 2)
NovMSMx   <- round(cellStats(Nov_CropSM, 'max'), 2)
NovMSMn   <- round(cellStats(Nov_CropSM, 'min'), 2)

Dec <- raster(Mean_CLIM[43])
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropSM <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
DecMSM   <- round(cellStats(Dec_CropSM, 'mean'), 2)
DecMSMx   <- round(cellStats(Dec_CropSM, 'max'), 2)
DecMSMn   <- round(cellStats(Dec_CropSM, 'min'), 2)

Ann <- raster(Mean_CLIM[40])
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropSM <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
AnnMSM   <- round(cellStats(Ann_CropSM, 'mean'), 2)
AnnMSMx   <- round(cellStats(Ann_CropSM, 'max'), 2)
AnnMSMn   <- round(cellStats(Ann_CropSM, 'min'), 2)

MEANSM2 <-
  c(
    JanMSM,
    FebMSM,
    MarMSM,
    AprMSM,
    MayMSM,
    JunMSM,
    JulMSM,
    AugMSM,
    SepMSM,
    OctMSM,
    NovMSM,
    DecMSM,
    AnnMSM
  )
MEANSM3 <- MEANSM2 * 100

#GSM thresholds for scale Bar
SMUP <-
  max(
    JanMSMx,
    FebMSMx,
    MarMSMx,
    AprMSMx,
    MayMSMx,
    JunMSMx,
    JulMSMx,
    AugMSMx,
    SepMSMx,
    OctMSMx,
    NovMSMx,
    DecMSMx,
    AnnMSMx
  )
SMLo <-
  min(
    JanMSMn,
    FebMSMn,
    MarMSMn,
    AprMSMn,
    MayMSMn,
    JunMSMn,
    JulMSMn,
    AugMSMn,
    SepMSMn,
    OctMSMn,
    NovMSMn,
    DecMSMn,
    AnnMSMn
  )

Cell.CL_Year[8, 2:14] <- MEANSM3




##########   EVAPOTRANSIPRATION
print("PDKE: 44, EVAPOTRANSIPRATION")

Jan <- raster(Mean_CLIM[19])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropET <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Jan_CropET <- Jan_CropET  * 0.0393701
}
JanMET   <- round(cellStats(Jan_CropET, 'mean'), 0)
JanMETx   <- round(cellStats(Jan_CropET, 'max'), 0)
JanMETn   <- round(cellStats(Jan_CropET, 'min'), 0)

Feb <- raster(Mean_CLIM[18])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropET <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Feb_CropET <- Feb_CropET  * 0.0393701
}
FebMET   <- round(cellStats(Feb_CropET, 'mean'), 0)
FebMETx   <- round(cellStats(Feb_CropET, 'max'), 0)
FebMETn   <- round(cellStats(Feb_CropET, 'min'), 0)

Mar <- raster(Mean_CLIM[22])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropET <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Mar_CropET <- Mar_CropET  * 0.0393701
}
MarMET   <- round(cellStats(Mar_CropET, 'mean'), 0)
MarMETx   <- round(cellStats(Mar_CropET, 'max'), 0)
MarMETn   <- round(cellStats(Mar_CropET, 'min'), 0)

Apr <- raster(Mean_CLIM[15])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropET <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Apr_CropET <- Apr_CropET  * 0.0393701
}
AprMET   <- round(cellStats(Apr_CropET, 'mean'), 0)
AprMETx   <- round(cellStats(Apr_CropET, 'max'), 0)
AprMETn   <- round(cellStats(Apr_CropET, 'min'), 0)

May <- raster(Mean_CLIM[23])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropET <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  May_CropET <- May_CropET  * 0.0393701
}
MayMET   <- round(cellStats(May_CropET, 'mean'), 0)
MayMETx   <- round(cellStats(May_CropET, 'max'), 0)
MayMETn   <- round(cellStats(May_CropET, 'min'), 0)

Jun <- raster(Mean_CLIM[21])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropET <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Jun_CropET <- Jun_CropET  * 0.0393701
}
JunMET   <- round(cellStats(Jun_CropET, 'mean'), 0)
JunMETx   <- round(cellStats(Jun_CropET, 'max'), 0)
JunMETn   <- round(cellStats(Jun_CropET, 'min'), 0)

Jul <- raster(Mean_CLIM[20])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropET <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Jul_CropET <- Jul_CropET  * 0.0393701
}
JulMET   <- round(cellStats(Jul_CropET, 'mean'), 0)
JulMETx   <- round(cellStats(Jul_CropET, 'max'), 0)
JulMETn   <- round(cellStats(Jul_CropET, 'min'), 0)

Aug <- raster(Mean_CLIM[16])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropET <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Aug_CropET <- Aug_CropET  * 0.0393701
}
AugMET   <- round(cellStats(Aug_CropET, 'mean'), 0)
AugMETx   <- round(cellStats(Aug_CropET, 'max'), 0)
AugMETn   <- round(cellStats(Aug_CropET, 'min'), 0)

Sep <- raster(Mean_CLIM[26])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropET <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Sep_CropET <- Sep_CropET  * 0.0393701
}
SepMET   <- round(cellStats(Sep_CropET, 'mean'), 0)
SepMETx   <- round(cellStats(Sep_CropET, 'max'), 0)
SepMETn   <- round(cellStats(Sep_CropET, 'min'), 0)

Oct <- raster(Mean_CLIM[25])
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropET <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Oct_CropET <- Oct_CropET  * 0.0393701
}
OctMET   <- round(cellStats(Oct_CropET, 'mean'), 0)
OctMETx   <- round(cellStats(Oct_CropET, 'max'), 0)
OctMETn   <- round(cellStats(Oct_CropET, 'min'), 0)

Nov <- raster(Mean_CLIM[24])
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropET <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Nov_CropET <- Nov_CropET  * 0.0393701
}
NovMET   <- round(cellStats(Nov_CropET, 'mean'), 0)
NovMETx   <- round(cellStats(Nov_CropET, 'max'), 0)
NovMETn   <- round(cellStats(Nov_CropET, 'min'), 0)

Dec <- raster(Mean_CLIM[17])
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropET <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Dec_CropET <- Dec_CropET  * 0.0393701
}
DecMET   <- round(cellStats(Dec_CropET, 'mean'), 0)
DecMETx   <- round(cellStats(Dec_CropET, 'max'), 0)
DecMETn   <- round(cellStats(Dec_CropET, 'min'), 0)

Ann <- raster(Mean_CLIM[14])
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropET <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Ann_CropET <- Ann_CropET  * 0.0393701
}
AnnMET   <- round(cellStats(Ann_CropET, 'mean'), 0)
AnnMETx   <- round(cellStats(Ann_CropET, 'max'), 0)
AnnMETn   <- round(cellStats(Ann_CropET, 'min'), 0)

MEANET2 <-
  c(
    JanMET,
    FebMET,
    MarMET,
    AprMET,
    MayMET,
    JunMET,
    JulMET,
    AugMET,
    SepMET,
    OctMET,
    NovMET,
    DecMET,
    AnnMET
  )


#Get thresholds for scale Bar
ETUP <-
  max(
    JanMETx,
    FebMETx,
    MarMETx,
    AprMETx,
    MayMETx,
    JunMETx,
    JulMETx,
    AugMETx,
    SepMETx,
    OctMETx,
    NovMETx,
    DecMETx,
    AnnMETx
  )
ETLo <-
  min(
    JanMETn,
    FebMETn,
    MarMETn,
    AprMETn,
    MayMETn,
    JunMETn,
    JulMETn,
    AugMETn,
    SepMETn,
    OctMETn,
    NovMETn,
    DecMETn,
    AnnMETn
  )

Cell.CL_Year[7, 2:14] <- MEANET2

plot(Aug)
##########   Cloud Frequency
print("PDKE: 45, Cloud Frequency")

Jan <- raster(Mean_CLIM[6])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropCF <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
JanMCF   <- round(cellStats(Jan_CropCF, 'mean'), 2)
JanMCFx   <- round(cellStats(Jan_CropCF, 'max'), 2)
JanMCFn   <- round(cellStats(Jan_CropCF, 'min'), 2)

Feb <- raster(Mean_CLIM[5])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropCF <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
FebMCF   <- round(cellStats(Feb_CropCF, 'mean'), 2)
FebMCFx   <- round(cellStats(Feb_CropCF, 'max'), 2)
FebMCFn   <- round(cellStats(Feb_CropCF, 'min'), 2)

Mar <- raster(Mean_CLIM[9])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropCF <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
MarMCF   <- round(cellStats(Mar_CropCF, 'mean'), 2)
MarMCFx   <- round(cellStats(Mar_CropCF, 'max'), 2)
MarMCFn   <- round(cellStats(Mar_CropCF, 'min'), 2)

Apr <- raster(Mean_CLIM[2])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropCF <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
AprMCF   <- round(cellStats(Apr_CropCF, 'mean'), 2)
AprMCFx   <- round(cellStats(Apr_CropCF, 'max'), 2)
AprMCFn   <- round(cellStats(Apr_CropCF, 'min'), 2)

May <- raster(Mean_CLIM[10])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropCF <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
MayMCF   <- round(cellStats(May_CropCF, 'mean'), 2)
MayMCFx   <- round(cellStats(May_CropCF, 'max'), 2)
MayMCFn   <- round(cellStats(May_CropCF, 'min'), 2)

Jun <- raster(Mean_CLIM[8])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropCF <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
JunMCF   <- round(cellStats(Jun_CropCF, 'mean'), 2)
JunMCFx   <- round(cellStats(Jun_CropCF, 'max'), 2)
JunMCFn   <- round(cellStats(Jun_CropCF, 'min'), 2)

Jul <- raster(Mean_CLIM[7])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropCF <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
JulMCF   <- round(cellStats(Jul_CropCF, 'mean'), 2)
JulMCFx   <- round(cellStats(Jul_CropCF, 'max'), 2)
JulMCFn   <- round(cellStats(Jul_CropCF, 'min'), 2)

Aug <- raster(Mean_CLIM[3])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropCF <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
AugMCF   <- round(cellStats(Aug_CropCF, 'mean'), 2)
AugMCFx   <- round(cellStats(Aug_CropCF, 'max'), 2)
AugMCFn   <- round(cellStats(Aug_CropCF, 'min'), 2)

Sep <- raster(Mean_CLIM[13])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropCF <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
SepMCF   <- round(cellStats(Sep_CropCF, 'mean'), 2)
SepMCFx   <- round(cellStats(Sep_CropCF, 'max'), 2)
SepMCFn   <- round(cellStats(Sep_CropCF, 'min'), 2)

Oct <- raster(Mean_CLIM[12])
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropCF <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
OctMCF   <- round(cellStats(Oct_CropCF, 'mean'), 2)
OctMCFx   <- round(cellStats(Oct_CropCF, 'max'), 2)
OctMCFn   <- round(cellStats(Oct_CropCF, 'min'), 2)

Nov <- raster(Mean_CLIM[11])
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropCF <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
NovMCF   <- round(cellStats(Nov_CropCF, 'mean'), 2)
NovMCFx   <- round(cellStats(Nov_CropCF, 'max'), 2)
NovMCFn   <- round(cellStats(Nov_CropCF, 'min'), 2)

Dec <- raster(Mean_CLIM[4])
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropCF <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
DecMCF   <- round(cellStats(Dec_CropCF, 'mean'), 2)
DecMCFx   <- round(cellStats(Dec_CropCF, 'max'), 2)
DecMCFn   <- round(cellStats(Dec_CropCF, 'min'), 2)

Ann <- raster(Mean_CLIM[1])
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropCF <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
AnnMCF   <- round(cellStats(Ann_CropCF, 'mean'), 2)
AnnMCFx   <- round(cellStats(Ann_CropCF, 'max'), 2)
AnnMCFn   <- round(cellStats(Ann_CropCF, 'min'), 2)

MEANCF2 <-
  c(
    JanMCF,
    FebMCF,
    MarMCF,
    AprMCF,
    MayMCF,
    JunMCF,
    JulMCF,
    AugMCF,
    SepMCF,
    OctMCF,
    NovMCF,
    DecMCF,
    AnnMCF
  )
MEANCF3 <- MEANCF2 * 100

#Get thresholds for scale Bar
CFUP <-
  max(
    JanMCFx,
    FebMCFx,
    MarMCFx,
    AprMCFx,
    MayMCFx,
    JunMCFx,
    JulMCFx,
    AugMCFx,
    SepMCFx,
    OctMCFx,
    NovMCFx,
    DecMCFx,
    AnnMCFx
  )
CFLO <-
  min(
    JanMCFn,
    FebMCFn,
    MarMCFn,
    AprMCFn,
    MayMCFn,
    JunMCFn,
    JulMCFn,
    AugMCFn,
    SepMCFn,
    OctMCFn,
    NovMCFn,
    DecMCFn,
    AnnMCFn
  )

Cell.CL_Year[6, 2:14] <- MEANCF3

##########

##########   MEAN TEMPERATURE
print("PDKE: 46, MEAN TEMPERATURE")

Jan <- raster(Mean_CLIM[71])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropTA <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Jan_CropTA <- (Jan_CropTA * 1.8) + 32
}
JanMTA   <- round(cellStats(Jan_CropTA, 'mean'), 1)
JanMTAx   <- round(cellStats(Jan_CropTA, 'max'), 1)
JanMTAn   <- round(cellStats(Jan_CropTA, 'min'), 1)

Feb <- raster(Mean_CLIM[70])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropTA <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Feb_CropTA <- (Feb_CropTA * 1.8) + 32
}
FebMTA   <- round(cellStats(Feb_CropTA, 'mean'), 1)
FebMTAx   <- round(cellStats(Feb_CropTA, 'max'), 1)
FebMTAn   <- round(cellStats(Feb_CropTA, 'min'), 1)

Mar <- raster(Mean_CLIM[74])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropTA <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Mar_CropTA <- (Mar_CropTA * 1.8) + 32
}
MarMTA   <- round(cellStats(Mar_CropTA, 'mean'), 1)
MarMTAx   <- round(cellStats(Mar_CropTA, 'max'), 1)
MarMTAn   <- round(cellStats(Mar_CropTA, 'min'), 1)

Apr <- raster(Mean_CLIM[67])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropTA <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Apr_CropTA <- (Apr_CropTA * 1.8) + 32
}
AprMTA   <- round(cellStats(Apr_CropTA, 'mean'), 1)
AprMTAx   <- round(cellStats(Apr_CropTA, 'max'), 1)
AprMTAn   <- round(cellStats(Apr_CropTA, 'min'), 1)

May <- raster(Mean_CLIM[75])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropTA <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  May_CropTA <- (May_CropTA * 1.8) + 32
}
MayMTA   <- round(cellStats(May_CropTA, 'mean'), 1)
MayMTAx   <- round(cellStats(May_CropTA, 'max'), 1)
MayMTAn   <- round(cellStats(May_CropTA, 'min'), 1)

Jun <- raster(Mean_CLIM[73])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropTA <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Jun_CropTA <- (Jun_CropTA * 1.8) + 32
}
JunMTA   <- round(cellStats(Jun_CropTA, 'mean'), 1)
JunMTAx   <- round(cellStats(Jun_CropTA, 'max'), 1)
JunMTAn   <- round(cellStats(Jun_CropTA, 'min'), 1)

Jul <- raster(Mean_CLIM[72])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropTA <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Jul_CropTA <- (Jul_CropTA * 1.8) + 32
}
JulMTA   <- round(cellStats(Jul_CropTA, 'mean'), 1)
JulMTAx   <- round(cellStats(Jul_CropTA, 'max'), 1)
JulMTAn   <- round(cellStats(Jul_CropTA, 'min'), 1)

Aug <- raster(Mean_CLIM[68])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropTA <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Aug_CropTA <- (Aug_CropTA * 1.8) + 32
}
AugMTA   <- round(cellStats(Aug_CropTA, 'mean'), 1)
AugMTAx   <- round(cellStats(Aug_CropTA, 'max'), 1)
AugMTAn   <- round(cellStats(Aug_CropTA, 'min'), 1)

Sep <- raster(Mean_CLIM[78])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropTA <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Sep_CropTA <- (Sep_CropTA * 1.8) + 32
}
SepMTA   <- round(cellStats(Sep_CropTA, 'mean'), 1)
SepMTAx   <- round(cellStats(Sep_CropTA, 'max'), 1)
SepMTAn   <- round(cellStats(Sep_CropTA, 'min'), 1)

Oct <- raster(Mean_CLIM[77])
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropTA <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Oct_CropTA <- (Oct_CropTA * 1.8) + 32
}
OctMTA   <- round(cellStats(Oct_CropTA, 'mean'), 1)
OctMTAx   <- round(cellStats(Oct_CropTA, 'max'), 1)
OctMTAn   <- round(cellStats(Oct_CropTA, 'min'), 1)

Nov <- raster(Mean_CLIM[76])
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropTA <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Nov_CropTA <- (Nov_CropTA * 1.8) + 32
}
NovMTA   <- round(cellStats(Nov_CropTA, 'mean'), 1)
NovMTAx   <- round(cellStats(Nov_CropTA, 'max'), 1)
NovMTAn   <- round(cellStats(Nov_CropTA, 'min'), 1)

Dec <- raster(Mean_CLIM[69])
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropTA <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Dec_CropTA <- (Dec_CropTA * 1.8) + 32
}
DecMTA   <- round(cellStats(Dec_CropTA, 'mean'), 1)
DecMTAx   <- round(cellStats(Dec_CropTA, 'max'), 1)
DecMTAn   <- round(cellStats(Dec_CropTA, 'min'), 1)

Ann <- raster(Mean_CLIM[66])
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropTA <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Ann_CropTA <- (Ann_CropTA * 1.8) + 32
}
AnnMTA   <- round(cellStats(Ann_CropTA, 'mean'), 1)
AnnMTAx   <- round(cellStats(Ann_CropTA, 'max'), 1)
AnnMTAn   <- round(cellStats(Ann_CropTA, 'min'), 1)

MEANTA2 <-
  c(
    JanMTA,
    FebMTA,
    MarMTA,
    AprMTA,
    MayMTA,
    JunMTA,
    JulMTA,
    AugMTA,
    SepMTA,
    OctMTA,
    NovMTA,
    DecMTA,
    AnnMTA
  )

#Get thresholds for scale Bar
TaUP <-
  max(
    JanMTAx,
    FebMTAx,
    MarMTAx,
    AprMTAx,
    MayMTAx,
    JunMTAx,
    JulMTAx,
    AugMTAx,
    SepMTAx,
    OctMTAx,
    NovMTAx,
    DecMTAx,
    AnnMTAx
  )
TaLo <-
  min(
    JanMTAn,
    FebMTAn,
    MarMTAn,
    AprMTAn,
    MayMTAn,
    JunMTAn,
    JulMTAn,
    AugMTAn,
    SepMTAn,
    OctMTAn,
    NovMTAn,
    DecMTAn,
    AnnMTAn
  )

Cell.CL_Year[3, 2:14] <- MEANTA2

########## MAX TEMPERATURE
print("PDKE: 47, MAX TEMPERATURE")

Jan <- raster(Mean_CLIM[84])
crs(Jan) <- crs(EXAMP)
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropTX <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Jan_CropTX <- (Jan_CropTX * 1.8) + 32
}
JanMTX  <- round(cellStats(Jan_CropTX , 'mean'), 1)
JanMTXx   <- round(cellStats(Jan_CropTX, 'max'), 1)
JanMTXn   <- round(cellStats(Jan_CropTX, 'min'), 1)

Feb <- raster(Mean_CLIM[83])
crs(Feb) <- crs(EXAMP)
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropTX  <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Feb_CropTX <- (Feb_CropTX * 1.8) + 32
}
FebMTX   <- round(cellStats(Feb_CropTX , 'mean'), 1)
FebMTXx   <- round(cellStats(Feb_CropTX, 'max'), 1)
FebMTXn   <- round(cellStats(Feb_CropTX, 'min'), 1)

Mar <- raster(Mean_CLIM[87])
crs(Mar) <- crs(EXAMP)
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropTX  <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Mar_CropTX <- (Mar_CropTX * 1.8) + 32
}
MarMTX   <- round(cellStats(Mar_CropTX , 'mean'), 1)
MarMTXx   <- round(cellStats(Mar_CropTX, 'max'), 1)
MarMTXn   <- round(cellStats(Mar_CropTX, 'min'), 1)

Apr <- raster(Mean_CLIM[80])
crs(Apr) <- crs(EXAMP)
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropTX  <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Apr_CropTX <- (Apr_CropTX * 1.8) + 32
}
AprMTX   <- round(cellStats(Apr_CropTX , 'mean'), 1)
AprMTXx   <- round(cellStats(Apr_CropTX, 'max'), 1)
AprMTXn   <- round(cellStats(Apr_CropTX, 'min'), 1)

May <- raster(Mean_CLIM[88])
crs(May) <- crs(EXAMP)
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropTX  <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  May_CropTX <- (May_CropTX * 1.8) + 32
}
MayMTX    <- round(cellStats(May_CropTX, 'mean'), 1)
MayMTXx   <- round(cellStats(May_CropTX, 'max'), 1)
MayMTXn   <- round(cellStats(May_CropTX, 'min'), 1)

Jun <- raster(Mean_CLIM[86])
crs(Jun) <- crs(EXAMP)
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropTX  <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Jun_CropTX <- (Jun_CropTX * 1.8) + 32
}
JunMTX   <- round(cellStats(Jun_CropTX , 'mean'), 1)
JunMTXx   <- round(cellStats(Jun_CropTX, 'max'), 1)
JunMTXn   <- round(cellStats(Jun_CropTX, 'min'), 1)

Jul <- raster(Mean_CLIM[85])
crs(Jul) <- crs(EXAMP)
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropTX  <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Jul_CropTX <- (Jul_CropTX * 1.8) + 32
}
JulMTX    <- round(cellStats(Jul_CropTX, 'mean'), 1)
JulMTXx   <- round(cellStats(Jul_CropTX, 'max'), 1)
JulMTXn   <- round(cellStats(Jul_CropTX, 'min'), 1)

Aug <- raster(Mean_CLIM[81])
crs(Aug) <- crs(EXAMP)
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropTX  <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Aug_CropTX <- (Aug_CropTX * 1.8) + 32
}
AugMTX   <- round(cellStats(Aug_CropTX, 'mean'), 1)
AugMTXx   <- round(cellStats(Aug_CropTX, 'max'), 1)
AugMTXn   <- round(cellStats(Aug_CropTX, 'min'), 1)

Sep <- raster(Mean_CLIM[91])
crs(Sep) <- crs(EXAMP)
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropTX  <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Sep_CropTX <- (Sep_CropTX * 1.8) + 32
}
SepMTX   <- round(cellStats(Sep_CropTX, 'mean'), 1)
SepMTXx   <- round(cellStats(Sep_CropTX, 'max'), 1)
SepMTXn   <- round(cellStats(Sep_CropTX, 'min'), 1)

Oct <- raster(Mean_CLIM[90])
crs(Oct) <- crs(EXAMP)
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropTX  <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Oct_CropTX <- (Oct_CropTX * 1.8) + 32
}
OctMTX   <- round(cellStats(Oct_CropTX, 'mean'), 1)
OctMTXx   <- round(cellStats(Oct_CropTX, 'max'), 1)
OctMTXn   <- round(cellStats(Oct_CropTX, 'min'), 1)

Nov <- raster(Mean_CLIM[89])
crs(Nov) <- crs(EXAMP)
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropTX  <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Nov_CropTX <- (Nov_CropTX * 1.8) + 32
}
NovMTX   <- round(cellStats(Nov_CropTX, 'mean'), 1)
NovMTXx   <- round(cellStats(Nov_CropTX, 'max'), 1)
NovMTXn   <- round(cellStats(Nov_CropTX, 'min'), 1)

Dec <- raster(Mean_CLIM[82])
crs(Dec) <- crs(EXAMP)
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropTX  <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Dec_CropTX <- (Dec_CropTX * 1.8) + 32
}
DecMTX   <- round(cellStats(Dec_CropTX, 'mean'), 1)
DecMTXx   <- round(cellStats(Dec_CropTX, 'max'), 1)
DecMTXn   <- round(cellStats(Dec_CropTX, 'min'), 1)

Ann <- raster(Mean_CLIM[79])
crs(Ann) <- crs(EXAMP)
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropTX <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Ann_CropTX <- (Ann_CropTX * 1.8) + 32
}
AnnMTX   <- round(cellStats(Ann_CropTX, 'mean'), 1)
AnnMTXx   <- round(cellStats(Ann_CropTX, 'max'), 1)
AnnMTXn   <- round(cellStats(Ann_CropTX, 'min'), 1)

MEANTX2 <-
  c(
    JanMTX,
    FebMTX,
    MarMTX,
    AprMTX,
    MayMTX,
    JunMTX,
    JulMTX,
    AugMTX,
    SepMTX,
    OctMTX,
    NovMTX,
    DecMTX,
    AnnMTX
  )

##########   Add to dataframe

Cell.CL_Year[4, 2:14] <- MEANTX2

##########   Get thresholds for figures

TxUP <-
  max(
    JanMTXx,
    FebMTXx,
    MarMTXx,
    AprMTXx,
    MayMTXx,
    JunMTXx,
    JulMTXx,
    AugMTXx,
    SepMTXx,
    OctMTXx,
    NovMTXx,
    DecMTXx,
    AnnMTXx
  )
TxLO <-
  min(
    JanMTXn,
    FebMTXn,
    MarMTXn,
    AprMTXn,
    MayMTXn,
    JunMTXn,
    JulMTXn,
    AugMTXn,
    SepMTXn,
    OctMTXn,
    NovMTXn,
    DecMTXn,
    AnnMTXn
  )

##########    MIN TEMPERATURE
print("PDKE: 48, MIN TEMPERATURE")

Jan <- raster(Mean_CLIM[97])
crs(Jan) <- crs(EXAMP)
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropTN <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Jan_CropTN <- (Jan_CropTN * 1.8) + 32
}
JanMTN   <- round(cellStats(Jan_CropTN, 'mean'), 1)
JanMTNx   <- round(cellStats(Jan_CropTN , 'max'), 1)
JanMTNn   <- round(cellStats(Jan_CropTN , 'min'), 1)

Feb <- raster(Mean_CLIM[96])
crs(Feb) <- crs(EXAMP)
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropTN <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Feb_CropTN <- (Feb_CropTN * 1.8) + 32
}
FebMTN  <- round(cellStats(Feb_CropTN, 'mean'), 1)
FebMTNx   <- round(cellStats(Feb_CropTN , 'max'), 1)
FebMTNn   <- round(cellStats(Feb_CropTN , 'min'), 1)

Mar <- raster(Mean_CLIM[100])
crs(Mar) <- crs(EXAMP)
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropTN <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Mar_CropTN <- (Mar_CropTN * 1.8) + 32
}
MarMTN   <- round(cellStats(Mar_CropTN, 'mean'), 1)
MarMTNx   <- round(cellStats(Mar_CropTN , 'max'), 1)
MarMTNn   <- round(cellStats(Mar_CropTN , 'min'), 1)

Apr <- raster(Mean_CLIM[93])
crs(Apr) <- crs(EXAMP)
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropTN <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Apr_CropTN <- (Apr_CropTN * 1.8) + 32
}
AprMTN   <- round(cellStats(Apr_CropTN, 'mean'), 1)
AprMTNx   <- round(cellStats(Apr_CropTN , 'max'), 1)
AprMTNn   <- round(cellStats(Apr_CropTN , 'min'), 1)

May <- raster(Mean_CLIM[101])
crs(May) <- crs(EXAMP)
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropTN <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  May_CropTN <- (May_CropTN * 1.8) + 32
}
MayMTN   <- round(cellStats(May_CropTN, 'mean'), 1)
MayMTNx   <- round(cellStats(May_CropTN , 'max'), 1)
MayMTNn   <- round(cellStats(May_CropTN , 'min'), 1)

Jun <- raster(Mean_CLIM[99])
crs(Jun) <- crs(EXAMP)
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropTN <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Jun_CropTN <- (Jun_CropTN * 1.8) + 32
}
JunMTN   <- round(cellStats(Jun_CropTN, 'mean'), 1)
JunMTNx   <- round(cellStats(Jun_CropTN , 'max'), 1)
JunMTNn   <- round(cellStats(Jun_CropTN , 'min'), 1)

Jul <- raster(Mean_CLIM[98])
crs(Jul) <- crs(EXAMP)
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropTN <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Jul_CropTN <- (Jul_CropTN * 1.8) + 32
}
JulMTN   <- round(cellStats(Jul_CropTN, 'mean'), 1)
JulMTNx   <- round(cellStats(Jul_CropTN , 'max'), 1)
JulMTNn   <- round(cellStats(Jul_CropTN , 'min'), 1)

Aug <- raster(Mean_CLIM[94])
crs(Aug) <- crs(EXAMP)
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropTN <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Aug_CropTN <- (Aug_CropTN * 1.8) + 32
}
AugMTN   <- round(cellStats(Aug_CropTN, 'mean'), 1)
AugMTNx   <- round(cellStats(Aug_CropTN , 'max'), 1)
AugMTNn   <- round(cellStats(Aug_CropTN , 'min'), 1)

Sep <- raster(Mean_CLIM[104])
crs(Sep) <- crs(EXAMP)
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropTN <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Sep_CropTN <- (Sep_CropTN * 1.8) + 32
}
SepMTN   <- round(cellStats(Sep_CropTN, 'mean'), 1)
SepMTNx   <- round(cellStats(Sep_CropTN , 'max'), 1)
SepMTNn   <- round(cellStats(Sep_CropTN , 'min'), 1)

Oct <- raster(Mean_CLIM[103])
crs(Oct) <- crs(EXAMP)
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropTN <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Oct_CropTN <- (Oct_CropTN * 1.8) + 32
}
OctMTN   <- round(cellStats(Oct_CropTN, 'mean'), 1)
OctMTNx   <- round(cellStats(Oct_CropTN , 'max'), 1)
OctMTNn   <- round(cellStats(Oct_CropTN , 'min'), 1)

Nov <- raster(Mean_CLIM[102])
crs(Nov) <- crs(EXAMP)
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropTN <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Nov_CropTN <- (Nov_CropTN * 1.8) + 32
}
NovMTN   <- round(cellStats(Nov_CropTN, 'mean'), 1)
NovMTNx   <- round(cellStats(Nov_CropTN , 'max'), 1)
NovMTNn   <- round(cellStats(Nov_CropTN , 'min'), 1)

Dec <- raster(Mean_CLIM[95])
crs(Dec) <- crs(EXAMP)
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropTN <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Dec_CropTN <- (Dec_CropTN * 1.8) + 32
}
DecMTN   <- round(cellStats(Dec_CropTN, 'mean'), 1)
DecMTNx   <- round(cellStats(Dec_CropTN , 'max'), 1)
DecMTNn   <- round(cellStats(Dec_CropTN , 'min'), 1)

Ann <- raster(Mean_CLIM[92])
crs(Ann) <- crs(EXAMP)
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropTN <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Ann_CropTN <- (Ann_CropTN * 1.8) + 32
}
AnnMTN   <- round(cellStats(Ann_CropTN, 'mean'), 1)
AnnMTNx   <- round(cellStats(Ann_CropTN , 'max'), 1)
AnnMTNn   <- round(cellStats(Ann_CropTN , 'min'), 1)

########## Add to Table

MEANTN2 <-
  c(
    JanMTN,
    FebMTN,
    MarMTN,
    AprMTN,
    MayMTN,
    JunMTN,
    JulMTN,
    AugMTN,
    SepMTN,
    OctMTN,
    NovMTN,
    DecMTN,
    AnnMTN
  )
Cell.CL_Year[2, 2:14] <- MEANTN2

########## Get thresholds for figures

TnUP <-
  max(
    JanMTNx,
    FebMTNx,
    MarMTNx,
    AprMTNx,
    MayMTNx,
    JunMTNx,
    JulMTNx,
    AugMTNx,
    SepMTNx,
    OctMTNx,
    NovMTNx,
    DecMTNx,
    AnnMTNx
  )
TnLO <-
  min(
    JanMTNn,
    FebMTNn,
    MarMTNn,
    AprMTNn,
    MayMTNn,
    JunMTNn,
    JulMTNn,
    AugMTNn,
    SepMTNn,
    OctMTNn,
    NovMTNn,
    DecMTNn,
    AnnMTNn
  )

##########   RELATIVE HUMIDITY
print("PDKE: 49, RELATIVE HUMIDITY")

Jan <- raster(Mean_CLIM[32])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropRH <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
JanMRH    <- round(cellStats(Jan_CropRH , 'mean'), 0)
JanMRHx   <- round(cellStats(Jan_CropRH  , 'max'), 0)
JanMRHn   <- round(cellStats(Jan_CropRH  , 'min'), 0)

Feb <- raster(Mean_CLIM[31])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropRH  <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
FebMRH    <- round(cellStats(Feb_CropRH , 'mean'), 0)
FebMRHx   <- round(cellStats(Feb_CropRH  , 'max'), 0)
FebMRHn   <- round(cellStats(Feb_CropRH  , 'min'), 0)

Mar <- raster(Mean_CLIM[35])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropRH  <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
MarMRH    <- round(cellStats(Mar_CropRH , 'mean'), 0)
MarMRHx   <- round(cellStats(Mar_CropRH  , 'max'), 0)
MarMRHn   <- round(cellStats(Mar_CropRH  , 'min'), 0)

Apr <- raster(Mean_CLIM[28])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropRH  <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
AprMRH    <- round(cellStats(Apr_CropRH , 'mean'), 0)
AprMRHx   <- round(cellStats(Apr_CropRH  , 'max'), 0)
AprMRHn   <- round(cellStats(Apr_CropRH  , 'min'), 0)

May <- raster(Mean_CLIM[36])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropRH  <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
MayMRH   <- round(cellStats(May_CropRH, 'mean'), 0)
MayMRHx   <- round(cellStats(May_CropRH  , 'max'), 0)
MayMRHn   <- round(cellStats(May_CropRH  , 'min'), 0)

Jun <- raster(Mean_CLIM[34])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropRH  <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
JunMRH    <- round(cellStats(Jun_CropRH, 'mean'), 0)
JunMRHx   <- round(cellStats(Jun_CropRH  , 'max'), 0)
JunMRHn   <- round(cellStats(Jun_CropRH  , 'min'), 0)

Jul <- raster(Mean_CLIM[33])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropRH  <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
JulMRH   <- round(cellStats(Jul_CropRH, 'mean'), 0)
JulMRHx   <- round(cellStats(Jul_CropRH  , 'max'), 0)
JulMRHn   <- round(cellStats(Jul_CropRH  , 'min'), 0)

Aug <- raster(Mean_CLIM[29])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropRH  <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
AugMRH   <- round(cellStats(Aug_CropRH, 'mean'), 0)
AugMRHx   <- round(cellStats(Aug_CropRH  , 'max'), 0)
AugMRHn   <- round(cellStats(Aug_CropRH  , 'min'), 0)

Sep <- raster(Mean_CLIM[39])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropRH  <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
SepMRH <- round(cellStats(Sep_CropRH, 'mean'), 0)
SepMRHx   <- round(cellStats(Sep_CropRH  , 'max'), 0)
SepMRHn   <- round(cellStats(Sep_CropRH  , 'min'), 0)

Oct <- raster(Mean_CLIM[38])
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropRH  <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
OctMRH <- round(cellStats(Oct_CropRH, 'mean'), 0)
OctMRHx   <- round(cellStats(Oct_CropRH  , 'max'), 0)
OctMRHn   <- round(cellStats(Oct_CropRH  , 'min'), 0)

Nov <- raster(Mean_CLIM[37])
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropRH  <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
NovMRH <- round(cellStats(Nov_CropRH , 'mean'), 0)
NovMRHx   <- round(cellStats(Nov_CropRH  , 'max'), 0)
NovMRHn   <- round(cellStats(Nov_CropRH  , 'min'), 0)

Dec <- raster(Mean_CLIM[30])
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropRH  <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
DecMRH    <- round(cellStats(Dec_CropRH, 'mean'), 0)
DecMRHx   <- round(cellStats(Dec_CropRH  , 'max'), 0)
DecMRHn   <- round(cellStats(Dec_CropRH  , 'min'), 0)

ANN <- raster(Mean_CLIM[27])
ANN_Mask <- mask(x = ANN, mask = UNIT_X[[u]])
ANN_CropRH  <- crop(x = ANN_Mask, y = extent(UNIT_X[[u]]))
ANNMRH    <- round(cellStats(ANN_CropRH, 'mean'), 0)
ANNMRHx   <- round(cellStats(ANN_CropRH  , 'max'), 0)
ANNMRHn   <- round(cellStats(ANN_CropRH  , 'min'), 0)

##########   Add to Table

MEANRH2 <-
  c(
    JanMRH,
    FebMRH ,
    MarMRH ,
    AprMRH ,
    MayMRH ,
    JunMRH ,
    JulMRH ,
    AugMRH ,
    SepMRH ,
    OctMRH ,
    NovMRH ,
    DecMRH,
    ANNMRH
  )
Cell.CL_Year[5, 2:14] <- MEANRH2
Cell.CL_Year

##########   Get thresholds for figures

RHUP <-
  max(
    JanMRHx,
    FebMRHx ,
    MarMRHx ,
    AprMRHx ,
    MayMRHx ,
    JunMRHx ,
    JulMRHx ,
    AugMRHx ,
    SepMRHx ,
    OctMRHx ,
    NovMRHx ,
    DecMRHx,
    ANNMRHx
  )
RHLO <-
  min(
    JanMRHn,
    FebMRHn ,
    MarMRHn ,
    AprMRHn ,
    MayMRHn ,
    JunMRHn ,
    JulMRHn ,
    AugMRHn ,
    SepMRHn ,
    OctMRHn ,
    NovMRHn ,
    DecMRHn,
    ANNMRHn
  )




#########   OTHER CLIMATE VARIABLES
print("PDKE: 50, OTHER CLIMATE VARIABLES")

### Read in the Annual rasters for each variable only and extract min and max to show spatial range.
# Have extracted the min max from Annual layers already This next step can be streamlined

#See Varible "ANN_CropRH" for example
RH_P <- raster(Mean_CLIM[27])
RH_P_Mask <- mask(x = RH_P, mask = UNIT_X[[u]])
RH_P_Crop <- crop(x = RH_P_Mask, y = extent(UNIT_X[[u]]))
RH_P_M   <- round(cellStats(RH_P_Crop, 'mean'), 0)
RHUPA   <- round(cellStats(RH_P_Crop  , 'max'), 1)
RHLOA   <- round(cellStats(RH_P_Crop  , 'min'), 1)

Tair_P <- raster(Mean_CLIM[66])
names(Tair_P) = "Mean Air Temp."
Tair_P_Mask <- mask(x = Tair_P, mask = UNIT_X[[u]])
Tair_P_Crop <- crop(x = Tair_P_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Tair_P_Crop <- (Tair_P_Crop * 1.8) + 32
}
Tair_P_M   <- round(cellStats(Tair_P_Crop,  'mean'), 1)
Tair_PName <- names(Tair_P)

SM_P <- raster(Mean_CLIM[40])
names(SM_P) = "Soil Moisture"
SM_P_Mask <- mask(x = SM_P, mask = UNIT_X[[u]])
SM_P_Crop <- crop(x = SM_P_Mask, y = extent(UNIT_X[[u]]))
SM_P_M   <- round(cellStats(SM_P_Crop, 'mean'), 2)
SM_PName <- names(SM_P)
SMUP   <- round(cellStats(SM_P_Crop  , 'max'), 2)
SMLO   <- round(cellStats(SM_P_Crop  , 'min'), 2)

CF_P <- raster(Mean_CLIM[1])
names(CF_P) = "Cloud Freq"
CF_P_Mask <- mask(x = CF_P, mask = UNIT_X[[u]])
CF_P_Crop <- crop(x = CF_P_Mask, y = extent(UNIT_X[[u]]))
CF_P_M   <- round(cellStats(CF_P_Crop, 'mean'), 2)
CF_PName <- names(CF_P)
CFUP   <- round(cellStats(CF_P_Crop  , 'max'), 2)
CFLO   <- round(cellStats(CF_P_Crop  , 'min'), 2)

Tmax_P <- raster(Mean_CLIM[79])
crs(Tmax_P)  <- crs(EXAMP)
names(Tmax_P) = "Max Air Temp."
Tmax_P_Mask <- mask(x = Tmax_P, mask = UNIT_X[[u]])
Tmax_P_Crop <- crop(x = Tmax_P_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Tmax_P_Crop <- (Tmax_P_Crop  * 1.8) + 32
}
Tmax_P_M   <- round(cellStats(Tmax_P_Crop, 'mean'), 1)
Tmax_PName <- names(Tmax_P)

Tmin_P <- raster(Mean_CLIM[92])
crs(Tmin_P)  <- crs(EXAMP)
names(Tmin_P) = "Mean Air Temp."
Tmin_P_Mask <- mask(x = Tmin_P, mask = UNIT_X[[u]])
Tmin_P_Crop <- crop(x = Tmin_P_Mask, y = extent(UNIT_X[[u]]))
if (TUnit == "°F") {
  Tmin_P_Crop <- (Tmin_P_Crop  * 1.8) + 32
}
Tmin_P_M   <- round(cellStats(Tmin_P_Crop, 'mean'), 1)
Tmin_PName <- names(Tmin_P)

KD_P <- raster(Mean_CLIM[53])
names(KD_P) = "Solar Radiation."
KD_P_Mask <- mask(x = KD_P, mask = UNIT_X[[u]])
KD_P_Crop <- crop(x = KD_P_Mask, y = extent(UNIT_X[[u]]))
KD_P_M   <- round(cellStats(KD_P_Crop, 'mean'), 0)
KD_PName <- names(KD_P)
KDUP   <- round(cellStats(KD_P_Crop  , 'max'), 0)
KDLO   <- round(cellStats(KD_P_Crop  , 'min'), 0)

ET_P <- raster(Mean_CLIM[14])
names(ET_P) = "Evaporation"
ET_P_Mask <- mask(x = ET_P, mask = UNIT_X[[u]])
ET_P_Crop <- crop(x = ET_P_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  ET_P_Crop <- ET_P_Crop * 0.0393701
}
ET_P_M   <- round(cellStats(ET_P_Crop , 'mean'), 0)
ET_PName <- names(ET_P)
ETUP   <- round(cellStats(ET_P_Crop  , 'max'), 0)
ETLO   <- round(cellStats(ET_P_Crop  , 'min'), 0)

VPD_P <- raster(Mean_CLIM[105])
names(VPD_P) = "VPD"
VPD_P_Mask <- mask(x = VPD_P, mask = UNIT_X[[u]])
VPD_P_Crop <- crop(x = VPD_P_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  VPD_P_Crop <- VPD_P_Crop * 0.0393701
}
VPD_P_M   <- round(cellStats(VPD_P_Crop, 'mean'), 0)
VPD_PName <- names(VPD_P)
VPDUP   <- round(cellStats(VPD_P_Crop  , 'max'), 0)
VPDLO   <- round(cellStats(VPD_P_Crop  , 'min'), 0)

WS_P <-
  raster(paste0(INPUTS_FOLDER, "Mean Climate/Wind/wind_sd_avg.tif"))
names(WS_P) = "Windspeed"
WS_P_Mask <- mask(x = WS_P, mask = UNIT_X[[u]])
WS_P_Crop <- crop(x = WS_P_Mask, y = extent(UNIT_X[[u]]))
#convert from m/sec to mph
WS_P_Crop <- WS_P_Crop * 2.237
WS_P_M   <- round(cellStats(WS_P_Crop, 'mean'), 0)
WS_PName <- names(WS_P)
WSUP   <- round(cellStats(WS_P_Crop  , 'max'), 1)
WSLO   <- round(cellStats(WS_P_Crop  , 'min'), 1)

##########   MEAN RAINFALL Data
print("PDKE: 51, MEAN RAINFALL Data")

Jan <- raster(MeanRF_ALL[1])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropRF <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Jan_CropRF  <- Jan_CropRF  * 0.0393701
}
JanMRF    <- round(cellStats(Jan_CropRF,  'mean'), 1)
JanMRFx   <- round(cellStats(Jan_CropRF  , 'max'), 1)
JanMRFn   <- round(cellStats(Jan_CropRF  , 'min'), 1)

print("PDKE: 51A")

Feb <- raster(MeanRF_ALL[2])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropRF <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Feb_CropRF  <- Feb_CropRF  * 0.0393701
}
FebMRF   <- round(cellStats(Feb_CropRF,  'mean'), 1)
FebMRFx   <- round(cellStats(Feb_CropRF  , 'max'), 1)
FebMRFn   <- round(cellStats(Feb_CropRF  , 'min'), 1)

Mar <- raster(MeanRF_ALL[3])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropRF <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Mar_CropRF  <- Mar_CropRF  * 0.0393701
}
MarMRF   <- round(cellStats(Mar_CropRF,  'mean'), 1)
MarMRFx   <- round(cellStats(Mar_CropRF  , 'max'), 1)
MarMRFn   <- round(cellStats(Mar_CropRF  , 'min'), 1)

Apr <- raster(MeanRF_ALL[4])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropRF <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Apr_CropRF  <- Apr_CropRF  * 0.0393701
}
AprMRF   <- round(cellStats(Apr_CropRF,  'mean'), 1)
AprMRFx   <- round(cellStats(Apr_CropRF  , 'max'), 1)
AprMRFn   <- round(cellStats(Apr_CropRF  , 'min'), 1)

May <- raster(MeanRF_ALL[5])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropRF <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  May_CropRF  <- May_CropRF  * 0.0393701
}
MayMRF   <- round(cellStats(May_CropRF,  'mean'), 1)
MayMRFx   <- round(cellStats(May_CropRF  , 'max'), 1)
MayMRFn   <- round(cellStats(May_CropRF  , 'min'), 1)

Jun <- raster(MeanRF_ALL[6])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropRF <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Jun_CropRF  <- Jun_CropRF  * 0.0393701
}
JunMRF   <- round(cellStats(Jun_CropRF,  'mean'), 1)
JunMRFx   <- round(cellStats(Jun_CropRF  , 'max'), 1)
JunMRFn   <- round(cellStats(Jun_CropRF  , 'min'), 1)

Jul <- raster(MeanRF_ALL[7])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropRF <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Jul_CropRF  <- Jul_CropRF  * 0.0393701
}
JulMRF   <- round(cellStats(Jul_CropRF,  'mean'), 1)
JulMRFx   <- round(cellStats(Jul_CropRF  , 'max'), 1)
JulMRFn   <- round(cellStats(Jul_CropRF  , 'min'), 1)

Aug <- raster(MeanRF_ALL[8])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropRF <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Aug_CropRF  <- Aug_CropRF  * 0.0393701
}
AugMRF   <- round(cellStats(Aug_CropRF,  'mean'), 1)
AugMRFx   <- round(cellStats(Aug_CropRF  , 'max'), 1)
AugMRFn   <- round(cellStats(Aug_CropRF  , 'min'), 1)

Sep <- raster(MeanRF_ALL[9])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropRF <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Sep_CropRF  <- Sep_CropRF  * 0.0393701
}
SepMRF   <- round(cellStats(Sep_CropRF, 'mean'), 1)
SepMRFx   <- round(cellStats(Sep_CropRF  , 'max'), 1)
SepMRFn   <- round(cellStats(Sep_CropRF  , 'min'), 1)

Oct <- raster(MeanRF_ALL[10])
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropRF <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Oct_CropRF  <- Oct_CropRF  * 0.0393701
}
OctMRF   <- round(cellStats(Oct_CropRF, 'mean'), 1)
OctMRFx   <- round(cellStats(Oct_CropRF  , 'max'), 1)
OctMRFn   <- round(cellStats(Oct_CropRF  , 'min'), 1)

Nov <- raster(MeanRF_ALL[11])
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropRF <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Nov_CropRF  <- Nov_CropRF  * 0.0393701
}
NovMRF   <- round(cellStats(Nov_CropRF,  'mean'), 1)
NovMRFx   <- round(cellStats(Nov_CropRF  , 'max'), 1)
NovMRFn   <- round(cellStats(Nov_CropRF  , 'min'), 1)

Dec <- raster(MeanRF_ALL[12])
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropRF <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  Dec_CropRF  <- Dec_CropRF  * 0.0393701
}
DecMRF   <- round(cellStats(Dec_CropRF,  'mean'), 1)
DecMRFx   <- round(cellStats(Dec_CropRF  , 'max'), 1)
DecMRFn   <- round(cellStats(Dec_CropRF  , 'min'), 1)

ANN <- raster(MeanRF_ALL[13])
ANN_Mask <- mask(x = ANN, mask = UNIT_X[[u]])
ANN_CropRF <- crop(x = ANN_Mask, y = extent(UNIT_X[[u]]))
if (RFUnit == " in") {
  ANN_CropRF  <- ANN_CropRF  * 0.0393701
}
ANNMRFM   <- round(cellStats(ANN_CropRF, 'mean'), 0)
ANNUP   <- round(cellStats(ANN_CropRF  , 'max'), 0)
ANNLO   <- round(cellStats(ANN_CropRF  , 'min'), 0)

##########   Add to Table

MEANRF2 <-
  c(
    JanMRF,
    FebMRF,
    MarMRF,
    AprMRF,
    MayMRF,
    JunMRF,
    JulMRF,
    AugMRF,
    SepMRF,
    OctMRF,
    NovMRF,
    DecMRF,
    ANNMRFM
  )
Cell.CL_Year[1, 2:14] <- MEANRF2

########## Get Threshold for figures

RFUP <-
  max(
    JanMRFx,
    FebMRFx,
    MarMRFx,
    AprMRFx,
    MayMRFx,
    JunMRFx,
    JulMRFx,
    AugMRFx,
    SepMRFx,
    OctMRFx,
    NovMRFx,
    DecMRFx
  )
RFLO <-
  min(
    JanMRFn,
    FebMRFn,
    MarMRFn,
    AprMRFn,
    MayMRFn,
    JunMRFn,
    JulMRFn,
    AugMRFn,
    SepMRFn,
    OctMRFn,
    NovMRFn,
    DecMRFn
  )

########## Write to Output Folder

write.csv(
  Cell.CL_Year,
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Annual Climate.csv"),
  row.names = F
)
Cell.CL_Year
##########   CLIMO GRAPHS
print("PDKE: 52, CLIMO GRAPHS")

#build data frame with temperature and precipitation data
df <-
  as.data.frame(c(
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec"
  ))
colnames(df) <- c("month")
Cell.CL_Year
df$month <- factor(df$month, levels = month.abb)
PREC <- as.numeric(as.vector(Cell.CL_Year[1, 2:13]))
TA <- as.numeric(as.vector(Cell.CL_Year[3, 2:13]))

df$PREC <- PREC
df$TA <- TA

#Set X and Y Limits for Graph 20% higher than the max.
XPR <- max(PREC) + (max(PREC) * 0.2)
XTA <- max(TA) +   (max(TA) * 0.2)
NTA <- min(TA) - min(TA) * 0.2

#Make Cliamo Graph and Save to ouptput folder

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Climograph2.png"),
  width = 5 * dpi,
  height = 3 * dpi,
  res = dpi
)

par(mar = c(4, 4, 4, 4))
my_bar <- barplot(
  df$PREC,
  border = F ,
  names.arg = df$month ,
  las = 2 ,
  col = "darkblue" ,
  ylim = c(0, XPR) ,
  ylab = paste0("Rainfall [", RFUnit2, "]")
)#,

XPR
title(paste0("Monthly Climate: ", UNIT_Ns[u]),
      line = 0.8,
      cex.main = 1.5)
# text(my_bar, df$PREC+13 ,df$PREC ,cex=1)
par(new = TRUE)
plot(
  df$month,
  df$TA ,
  pch = 15 ,
  lty = 0,
  axes = FALSE,
  xlab = "",
  ylab = "",
  col = "red",
  ylim = c(NTA, XTA)
)
lines(
  x = df$month,
  y = df$TA,
  lwd = 2,
  col = "red"
)
points(
  x = df$month,
  y = df$TA,
  col = "red",
  pch = 16,
  cex = 2
)
mtext(
  paste0("Temperature [", TUnit, "]"),
  side = 4,
  line = 2.5,
  col = "red"
)
axis(
  4,
  ylim = c(NTA, XTA),
  col = "red",
  col.axis = "red",
  las = 1
)

# # Legend
# legend("topleft", legend = c("Rainfall", "Temperature"),
#        col = c("darkblue" , "red") ,
#        bty = "n", pch=c(15,20) , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.05, 0))

dev.off()

### Monthly Rainfall
print("PDKE: 53, Monthly Rainfall")

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Climograph_RF.png"),
  width = 5 * dpi,
  height = 3 * dpi,
  res = dpi
)

barplot(
  df$PREC,
  border = F ,
  names.arg = df$month ,
  las = 2 ,
  col = "darkblue" ,
  ylim = c(0, XPR) ,
  ylab = paste0("Rainfall (", RFUnit2, ".)")
)#,
title(paste0("Monthly Rainfall: ", UNIT_Ns[u]),
      line = 0.8,
      cex.main = 1.5)
par(new = TRUE)
plot(
  df$month,
  df$TA ,
  pch = 15 ,
  lty = 0,
  axes = FALSE,
  xlab = "",
  ylab = "",
  col = "red",
  ylim = c(NTA, XTA)
)

dev.off()

### Monthly Air Temp
col <- rgb(0.2, 0.2, 1, alpha = 0.15)

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Climograph_AT.png"),
  width = 5.2 * dpi,
  height = 3 * dpi,
  res = dpi
)

par(mar = c(4, 4, 4, 4))
barplot(
  df$PREC,
  border = F ,
  names.arg = df$month ,
  las = 2 ,
  col = col ,
  alpha = 0.5 ,
  ylim = c(0, XPR) ,
  ylab = paste0("Rainfall (", RFUnit2, ".)"),
  col.lab = "grey50"
)#,
title(paste0("Monthly Air Temperature"),
      line = 0,
      cex.main = 1.3)
par(new = TRUE)
plot(
  df$month,
  df$TA ,
  pch = 15 ,
  lty = 0,
  axes = FALSE,
  xlab = "",
  ylab = "",
  col = "red",
  ylim = c(NTA, XTA)
)
lines(
  x = df$month,
  y = df$TA,
  lwd = 2,
  col = "red"
)
points(
  x = df$month,
  y = df$TA,
  col = "red",
  pch = 16,
  cex = 2
)
mtext(
  paste0("Temperature (", TUnit, ")"),
  side = 4,
  line = 2.5,
  col = "red"
)
axis(
  4,
  ylim = c(NTA, XTA),
  col = "red",
  col.axis = "red",
  las = 1
)

dev.off()

########## Figure for Other Cliamte Variables
# Decide on a Break for Rainfall

#Mean Rainfall
print("PDKE: 53, Figure for Other Cliamte Variables ")

# Decide on a Break for Rainfall
RNGERF <- RFUP - RFLO
if (RNGERF > 8.99) {
  BI_brksRF <- round(seq(RFLO, RFUP, length = 10), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(10, "YlGnBu"))(50)
}
if (RNGERF < 9 && RNGERF > 7.99) {
  BI_brksRF <- round(seq(RFLO, RFUP, length = 9), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(9, "YlGnBu"))(50)
}
if (RNGERF < 8 && RNGERF > 6.99) {
  BI_brksRF <- round(seq(RFLO, RFUP, length = 8), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(8, "YlGnBu"))(50)
}
if (RNGERF < 7 && RNGERF > 5.99) {
  BI_brksRF <- round(seq(RFLO, RFUP, length = 7), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(7, "YlGnBu"))(50)
}
if (RNGERF < 6 && RNGERF > 4.99) {
  BI_brksRF <- round(seq(RFLO, RFUP, length = 5), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(5, "YlGnBu"))(50)
}
if (RNGERF < 5 && RNGERF > 3.99) {
  BI_brksRF <- round(seq(RFLO, RFUP, length = 5), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(5, "YlGnBu"))(50)
}
if (RNGERF < 4 && RNGERF > 2.99) {
  BI_brksRF <- round(seq(0, RFUP, length = 4), 1)
  
  colfuncRF <- colorRampPalette(brewer.pal(4, "YlGnBu"))(50)
}
if (RNGERF < 3 && RNGERF > 1.99) {
  BI_brksRF <- round(seq(0, RFUP, length = 3), 1)
  
  colfuncRF <- colorRampPalette(brewer.pal(3, "YlGnBu"))(50)
}
if (RNGERF < 2 && RNGERF > 0.99) {
  BI_brksRF <- round(seq(0, RFUP, length = 3), 1)
  
  colfuncRF <- colorRampPalette(brewer.pal(3, "YlGnBu"))(50)
}
if (RNGERF < 2 && RNGERF > 0) {
  BI_brksRF <- round(seq(0, 3, length = 3), 1)
  
  colfuncRF <- colorRampPalette(brewer.pal(3, "YlGnBu"))(50)
}

RNGERF

#Mean Rainfall
print("PDKE: 53, Mean Rainfall")

# Decide on a Break for Rainfall
RNGERFA <- ANNUP - ANNLO
if (RNGERFA >= 10) {
  BI_brksRFA <- round(seq(ANNLO, ANNUP, length = 11), 0)
  
  colfuncRFA <- colorRampPalette(brewer.pal(9, "YlGnBu"))(50)
}
if (RNGERFA < 10 && RNGERFA > 8.99) {
  BI_brksRFA <- round(seq(ANNLO, ANNUP, length = 10), 0)
  
  colfuncRFA <- colorRampPalette(brewer.pal(9, "YlGnBu"))(50)
}
if (RNGERFA < 9 && RNGERFA > 7.99) {
  BI_brksRFA <- round(seq(ANNLO, ANNUP, length = 9), 0)
  
  colfuncRFA <- colorRampPalette(brewer.pal(9, "YlGnBu"))(50)
}
if (RNGERFA < 8 && RNGERFA > 6.99) {
  BI_brksRFA <- round(seq(ANNLO, ANNUP, length = 8), 0)
  
  colfuncRFA <- colorRampPalette(brewer.pal(8, "YlGnBu"))(50)
}
if (RNGERFA < 7 && RNGERFA > 5.99) {
  BI_brksRFA <- round(seq(ANNLO, ANNUP, length = 7), 0)
  
  colfuncRFA <- colorRampPalette(brewer.pal(7, "YlGnBu"))(50)
}
if (RNGERFA < 6 && RNGERFA > 4.99) {
  BI_brksRFA <- round(seq(0, ANNUP, length = 6), 0)
  
  colfuncRFA <- colorRampPalette(brewer.pal(6, "YlGnBu"))(50)
}
if (RNGERFA < 5 && RNGERFA > 3.99) {
  BI_brksRFA <- round(seq(0, ANNUP, length = 5), 0)
  
  colfuncRFA <- colorRampPalette(brewer.pal(5, "YlGnBu"))(50)
}
if (RNGERFA < 4 && RNGERFA > 0.99) {
  BI_brksRFA <- round(seq(0, ANNUP, length = 3), 0)
  
  colfuncRFA <- colorRampPalette(brewer.pal(4, "YlGnBu"))(50)
}
if (RNGERFA < 4 && RNGERFA < 0.99) {
  BI_brksRFA <- round(seq(0, ANNUP, length = 3), 0)
  
  colfuncRFA <- colorRampPalette(brewer.pal(3, "YlGnBu"))(50)
}

#if (RNGERFA < 3){
#  BI_brksRFA<-round(seq(RFLO, RFUP, length = 3),0);
#  colfuncRFA<-colorRampPalette(brewer.pal(3,"YlGnBu"))(50)
#}

#Use Same Scale for TA
BI_brksTA <- round(seq(TnLO, TxUP, length = 9), 0)
colfuncTA <- colorRampPalette(brewer.pal(9, "YlOrRd"))(50)

#For RH
RNGERH <- RHUP - RHLO
if (RNGERH > 8) {
  BI_brksRH <- round(seq(RHLO, RHUP, length = 9), 0)
  colfuncRH <- colorRampPalette(brewer.pal(9, "BuPu"))(50)
}
# RH
if (RNGERH < 9 && RNGERH > 7.99) {
  BI_brksRH <- round(seq(RHLO, RHUP, length = 8), 0)
  colfuncRH <- colorRampPalette(brewer.pal(8, "BuPu"))(50)
}
if (RNGERH < 8 && RNGERH > 6.99) {
  BI_brksRH <- round(seq(RHLO, RHUP, length = 7), 0)
  colfuncRH <- colorRampPalette(brewer.pal(7, "BuPu"))(50)
}
if (RNGERH < 7 && RNGERH > 5.99) {
  BI_brksRH <- round(seq(RHLO, RHUP, length = 6), 0)
  colfuncRH <- colorRampPalette(brewer.pal(6, "BuPu"))(50)
}
if (RNGERH < 6 && RNGERH > 4.99) {
  BI_brksRH <- round(seq(RHLO, RHUP, length = 5), 0)
  colfuncRH <- colorRampPalette(brewer.pal(5, "BuPu"))(50)
}
if (RNGERH < 5 && RNGERH > 3.99) {
  BI_brksRH <- round(seq(RHLO, RHUP, length = 4), 0)
  colfuncRH <- colorRampPalette(brewer.pal(4, "BuPu"))(50)
}
if (RNGERH < 4 && RNGERH > 2.99) {
  BI_brksRH <- round(seq(RHLO, RHUP, length = 4), 0)
  colfuncRH <- colorRampPalette(brewer.pal(4, "BuPu"))(50)
}

#For SM
BI_brksSM <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
colfuncSM <- colorRampPalette(brewer.pal(9, "YlGn"))(50)
#For KD (Solar)
if ((KDUP - KDLO) >= 9) {
  BI_brksKD <- round(seq(KDLO, KDUP, length = 9), 0)
  colfuncKD <- colorRampPalette(brewer.pal(9, "OrRd"))(50)
}
if ((KDUP - KDLO) < 9) {
  BI_brksKD <- round(seq(KDLO, KDUP, length = (KDUP - KDLO)), 0)
  colfuncKD <- colorRampPalette(brewer.pal(1, "OrRd"))(50)
}
# if((KDUP - KDLO)<1){
#   BI_brksKD<-
#   colfuncKD <-colorRampPalette(brewer.pal(1,"OrRd"))(50)
# }

#For CF
BI_brksCF <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
colfuncCF <- colorRampPalette(brewer.pal(9, "PuRd"))(50)
#For ET
BI_brksET <- round(seq(ETLO, ETUP, length = 9), 0)
if (ETUP - ETLO < 9) {
  BI_brksET <- round(seq(ETLO, ETUP, length = 9), 1)
}
colfuncET <- colorRampPalette(brewer.pal(9, "PuBu"))(50)
#For WS
BI_brksWS <- round(seq(WSLO, WSUP, length = 9), 2)
colfuncWS <- colorRampPalette(brewer.pal(9, "Purples"))(50)

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Climate_less_RF.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
spplot(
  ANN_CropRF,
  col.regions = colfuncRF,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksRFA,
  main = list(
    label = paste0("Rainfall (", ANNMRFM , RFUnit, ")"),
    cex = 2
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 2)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
dev.off()

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Climate_less_TA.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
spplot(
  Tair_P_Crop,
  col.regions = colfuncTA,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksTA,
  main = list(
    label = paste0("Air Temperature (", Tair_P_M , TUnit2, ")"),
    cex = 2
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 2)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
dev.off()

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Climate_less_RH.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
spplot(
  RH_P_Crop,
  col.regions = colfuncRH,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksRH,
  main = list(
    label = paste0("Relative Humidity (", RH_P_M , " %)"),
    cex = 2
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 2)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
dev.off()


png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Climate_less_SR.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
spplot(
  KD_P_Crop,
  col.regions = colfuncKD,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksKD,
  main = list(
    label = paste0("Solar Radiation (", KD_P_M , " W/m2)"),
    cex = 2
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 2)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
dev.off()

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Climate_less_SM.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
spplot(
  SM_P_Crop,
  col.regions = colfuncSM,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksSM,
  main = list(
    label = paste0("Soil Moisture (", SM_P_M , ")"),
    cex = 2
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 2)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
dev.off()

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Climate_less_ET.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
spplot(
  ET_P_Crop,
  col.regions = colfuncET,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksET,
  main = list(
    label = paste0("Evapotranspiration (", ET_P_M , RFUnit, ")"),
    cex = 2
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 2)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
dev.off()

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Climate_less_WS.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)
spplot(
  WS_P_Crop,
  col.regions = colfuncWS,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksWS,
  main = list(
    label = paste0("Windspeed (", WS_P_M, " mph)"),
    cex = 2
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 2)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
dev.off()

##########  Temperature Maps Figure
print("PDKE: 54, Temperature Maps Figure ")

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " TA12.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

title1 = textGrob(
  paste("Monthly Temperature:", UNIT_Ns[u]),
  gp = gpar(
    col = "darkred",
    fontface = "bold",
    fontsize = 15
  )
)
grid.arrange(
  top = title1,
  spplot(
    Jan_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("JAN (", JanMTA , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Feb_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("FEB (", FebMTA , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Mar_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("MAR (", MarMTA  , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Apr_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("APR (", AprMTA  , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    May_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("MAY (", MayMTA , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Jun_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("JUN (", JunMTA  , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Jul_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("JUL (", JulMTA  , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Aug_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("AUG (", AugMTA , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Sep_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("SEP (", SepMTA  , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Oct_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("OCT (", OctMTA  , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Nov_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("NOV (", NovMTA  , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Dec_CropTA,
    col.regions = colfuncTA,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksTA,
    main = list(label = paste0("DEC (", DecMTA , TUnit2, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1))
)

dev.off()

##### Hottest and Coldest month maps

# find hottest and coldest months
print("PDKE: 55, find hottest and coldest months")

Cell.CL_Year
t <- Cell.CL_Year[3, 2:13]
min.col <- function(m, ...)
  max.col(-m, ...)
hm <- colnames(t)[min.col(t, ties.method = "first")]
cm <- colnames(t)[max.col(t, ties.method = "first")]

# Select correct month map based on dry and wet months
if (hm == "JAN") {
  hmap <- Jan_CropTA
  mnh <- "January"
}
if (hm == "FEB") {
  hmap <- Feb_CropTA
  mnh <- "February"
}
if (hm == "MAR") {
  hmap <- Mar_CropTA
  mnh <- "March"
}
if (hm == "APR") {
  hmap <- Apr_CropTA
  mnh <- "April"
}
if (hm == "MAY") {
  hmap <- May_CropTA
  mnh <- "May"
}
if (hm == "JUN") {
  hmap <- Jun_CropTA
  mnh <- "June"
}
if (hm == "JUL") {
  hmap <- Jul_CropTA
  mnh <- "July"
}
if (hm == "AUG") {
  hmap <- Aug_CropTA
  mnh <- "August"
}
if (hm == "SEP") {
  hmap <- Sep_CropTA
  mnh <- "September"
}
if (hm == "OCT") {
  hmap <- Oct_CropTA
  mnh <- "October"
}
if (hm == "NOV") {
  hmap <- Nov_CropTA
  mnh <- "November"
}
if (hm == "DEC") {
  hmap <- Dec_CropTA
  mnh <- "December"
}

if (cm == "JAN") {
  cmap <- Jan_CropTA
  mnc <- "January"
}
if (cm == "FEB") {
  cmap <- Feb_CropTA
  mnc <- "February"
}
if (cm == "MAR") {
  cmap <- Mar_CropTA
  mnc <- "March"
}
if (cm == "APR") {
  cmap <- Apr_CropTA
  mnc <- "April"
}
if (cm == "MAY") {
  cmap <- May_CropTA
  mnc <- "May"
}
if (cm == "JUN") {
  cmap <- Jun_CropTA
  mnc <- "June"
}
if (cm == "JUL") {
  cmap <- Jul_CropTA
  mnc <- "July"
}
if (cm == "AUG") {
  cmap <- Aug_CropTA
  mnc <- "August"
}
if (cm == "SEP") {
  cmap <- Sep_CropTA
  mnc <- "September"
}
if (cm == "OCT") {
  cmap <- Oct_CropTA
  mnc <- "October"
}
if (cm == "NOV") {
  cmap <- Nov_CropTA
  mnc <- "November"
}
if (cm == "DEC") {
  cmap <- Dec_CropTA
  mnc <- "December"
}

### make plots

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " TA_hm.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

spplot(
  hmap,
  col.regions = colfuncTA,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksTA,
  main = list(label = paste0(mnh, " Temp (°F)"), cex = 1.8),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 1.5)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))

dev.off()

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " TA_cm.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

spplot(
  cmap,
  col.regions = colfuncTA,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksTA,
  main = list(label = paste0(mnc, " Temp (°F)"), cex = 1.8),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 1.5)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))

dev.off()

##########   MEAN ANNUAL Relative Humidity 12-maps
print("PDKE: 56, MEAN ANNUAL Relative Humidity 12-maps")

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " RH12.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

title1 = textGrob(
  paste("Monthly Relative Humidity:", UNIT_Ns[u]),
  gp = gpar(
    col = "darkred",
    fontface = "bold",
    fontsize = 15
  )
)
grid.arrange(
  top = title1,
  spplot(
    Jan_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("JAN (", JanMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Feb_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("FEB (", FebMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Mar_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("MAR (", MarMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Apr_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("APR (", AprMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    May_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("MAY (", MayMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Jun_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("JUN (", JunMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Jul_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("JUL (", JulMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Aug_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("AUG (", AugMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Sep_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("SEP (", SepMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Oct_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("OCT (", OctMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Nov_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("NOV (", NovMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Dec_CropRH,
    col.regions = colfuncRH,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRH,
    main = list(label = paste0("DEC (", DecMRH , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1))
)
dev.off()

dpi = 300
##########   MEAN ANNUAL Rainfall 12-maps
#BI_brksRF<-round(seq(0, RFUP+0.3, length = 4),0) #For dry areas.
print("PDKE: 57, MEAN ANNUAL Rainfall 12-maps ")

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " RF12.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

title1 = textGrob(
  paste("Monthly Rainfall:", UNIT_Ns[u]),
  gp = gpar(
    col = "darkred",
    fontface = "bold",
    fontsize = 15
  )
)
grid.arrange(
  top = title1,
  spplot(
    Jan_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("JAN (", JanMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Feb_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("FEB (", FebMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Mar_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("MAR (", MarMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Apr_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("APR (", AprMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    May_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("MAY (", MayMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Jun_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("JUN (", JunMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Jul_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("JUL (", JulMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Aug_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("AUG (", AugMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Sep_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("SEP(", SepMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Oct_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("OCT (", OctMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Nov_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("NOV(", NovMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Dec_CropRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksRF,
    main = list(label = paste0("DEC (", DecMRF , RFUnit, ")"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1))
)

dev.off()

##### Driest and Wettest month maps
print("PDKE: 58, Driest and Wettest month maps")

# find driest and wettest months
Cell.CL_Year
d <- Cell.CL_Year[1, 2:13]
min.col <- function(m, ...)
  max.col(-m, ...)
dm <- colnames(d)[min.col(d, ties.method = "first")]
wm <- colnames(d)[max.col(d, ties.method = "first")]

# Select correct month map based on dry and wet months
if (dm == "JAN") {
  dmap <- Jan_CropRF
  mn <- "January"
}
if (dm == "FEB") {
  dmap <- Feb_CropRF
  mn <- "February"
}
if (dm == "MAR") {
  dmap <- Mar_CropRF
  mn <- "March"
}
if (dm == "APR") {
  dmap <- Apr_CropRF
  mn <- "April"
}
if (dm == "MAY") {
  dmap <- May_CropRF
  mn <- "May"
}
if (dm == "JUN") {
  dmap <- Jun_CropRF
  mn <- "June"
}
if (dm == "JUL") {
  dmap <- Jul_CropRF
  mn <- "July"
}
if (dm == "AUG") {
  dmap <- Aug_CropRF
  mn <- "August"
}
if (dm == "SEP") {
  dmap <- Sep_CropRF
  mn <- "September"
}
if (dm == "OCT") {
  dmap <- Oct_CropRF
  mn <- "October"
}
if (dm == "NOV") {
  dmap <- Nov_CropRF
  mn <- "November"
}
if (dm == "DEC") {
  dmap <- Dec_CropRF
  mn <- "December"
}

if (wm == "JAN") {
  wmap <- Jan_CropRF
  mnw <- "January"
}
if (wm == "FEB") {
  wmap <- Feb_CropRF
  mnw <- "February"
}
if (wm == "MAR") {
  wmap <- Mar_CropRF
  mnw <- "March"
}
if (wm == "APR") {
  wmap <- Apr_CropRF
  mnw <- "April"
}
if (wm == "MAY") {
  wmap <- May_CropRF
  mnw <- "May"
}
if (wm == "JUN") {
  wmap <- Jun_CropRF
  mnw <- "June"
}
if (wm == "JUL") {
  wmap <- Jul_CropRF
  mnw <- "July"
}
if (wm == "AUG") {
  wmap <- Aug_CropRF
  mnw <- "August"
}
if (wm == "SEP") {
  wmap <- Sep_CropRF
  mnw <- "September"
}
if (wm == "OCT") {
  wmap <- Oct_CropRF
  mnw <- "October"
}
if (wm == "NOV") {
  wmap <- Nov_CropRF
  mnw <- "November"
}
if (wm == "DEC") {
  wmap <- Dec_CropRF
  mnw <- "December"
}

### make plots

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " RF_dm.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

spplot(
  dmap,
  col.regions = colfuncRF,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksRF,
  main = list(label = paste0(mn, " Rainfall (in.)"), cex = 2),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 1.5)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))

dev.off()

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " RF_wm.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

spplot(
  wmap,
  col.regions = colfuncRF,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksRF,
  main = list(label = paste0(mnw, " Rainfall (in.)"), cex = 2),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 1.5)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))

dev.off()


#########   Seasonal Rainfall Maps
print("PDKE: 59, Seasonal Rainfall Maps ")

DrySeasonRF <-
  (May_CropRF + Jun_CropRF + Jul_CropRF + Aug_CropRF + Sep_CropRF + Oct_CropRF)
WetSeasonRF <-
  (Nov_CropRF + Dec_CropRF + Jan_CropRF + Feb_CropRF + Mar_CropRF + Apr_CropRF)

### Make wet and dry season average monthly rainfall maps and values
WetSM <- (WetSeasonRF / 6)
WetSMV <- round(cellStats(WetSM, 'mean'), 1)

DrySM <- (DrySeasonRF / 6)
DrySMV <- round(cellStats(DrySM, 'mean'), 1)

# get min and max values for figure scale
DryUPm   <- round(cellStats(DrySM, 'max'), 1)
DryLOm   <- round(cellStats(DrySM, 'min'), 1)
WetUPm   <- round(cellStats(WetSM, 'max'), 1)
WetLOm  <- round(cellStats(WetSM, 'min'), 1)

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

BI_brksSEA <- round(seq(SEALO, SEAUP, length = 10), 0)

RNGESEA <- SEAUP - SEALO
if (RNGESEA > 8.99) {
  BI_brksSEA <- round(seq(SEALO, SEAUP, length = 11), 1)
  
  colfuncRF <- colorRampPalette(brewer.pal(9, "YlGnBu"))(50)
}
if (RNGESEA < 9 && RNGESEA > 7.99) {
  BI_brksSEA <- round(seq(SEALO, SEAUP, length = 9), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(9, "YlGnBu"))(50)
}
if (RNGESEA < 8 && RNGESEA > 6.99) {
  BI_brksSEA <- round(seq(SEALO, SEAUP, length = 8), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(8, "YlGnBu"))(50)
}
if (RNGESEA < 7 && RNGESEA > 5.99) {
  BI_brksSEA <- round(seq(SEALO, SEAUP, length = 7), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(7, "YlGnBu"))(50)
}
if (RNGESEA < 6 && RNGESEA > 4.99) {
  BI_brksSEA <- round(seq(SEALO, SEAUP, length = 6), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(6, "YlGnBu"))(50)
}
if (RNGESEA < 5 && RNGESEA > 3.99) {
  BI_brksSEA <- round(seq(SEALO, SEAUP, length = 5), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(5, "YlGnBu"))(50)
}
if (RNGESEA < 4) {
  BI_brksSEA <- round(seq(SEALO, SEAUP, length = 4), 0)
  
  colfuncRF <- colorRampPalette(brewer.pal(4, "YlGnBu"))(50)
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

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " SeaRF.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

title1 = textGrob(
  paste("Seasonal Rainfall:", UNIT_Ns[u]),
  gp = gpar(
    col = "darkred",
    fontface = "bold",
    fontsize = 15
  )
)
grid.arrange(
  top = title1,
  spplot(
    WetSeasonRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksSEA,
    main = list(
      label = paste0("Wet Season (NOV-APR) ", WSeaMRF  , RFUnit),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ),
  spplot(
    DrySeasonRF,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksSEA,
    main = list(
      label = paste0("Dry Season (MAY-OCT) ", DSeaMRF  , RFUnit),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  )
)
dev.off()

# export seasonal monthly rainfall figure from above
png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " SeaMRF.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

title1 = textGrob(
  paste("Monthly Rainfall:", UNIT_Ns[u]),
  gp = gpar(
    col = "darkred",
    fontface = "bold",
    fontsize = 15
  )
)
grid.arrange(
  top = title1,
  spplot(
    WetSM,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksSEAm,
    main = list(
      label = paste0("Wet Season (NOV-APR) ", WetSMV  , RFUnit),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ),
  spplot(
    DrySM,
    col.regions = colfuncRF,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksSEAm,
    main = list(
      label = paste0("Dry Season (MAY-OCT) ", DrySMV  , RFUnit),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  )
)
dev.off()

### calculate the seasonal avg. monthly rainfall percentile vs. whole state
print("PDKE: 59A")
#stateRF
plot(stateRFM)
summary(stateRFMd)
#print(paste("PDKE: 59B", stateRFMd))
#dput(stateRFMd,file = "") # print out stateRFMd
print("PDKE: 59C")
per <-
  ecdf(stateRFMd$X) # was stateRFMd$layer, but there's no 'layer' in stateRFMd, so switched to 'X'
#per<-ecdf(stateRFMd$layer)
print("PDKE: 59D")

WetSMP <- round(per(WetSMV) * 100)
DrySMP <- round(per(DrySMV) * 100)
print("PDKE: 59E")

P <- data.frame(WetSMP, DrySMP)
print("PDKE: 59F")

# write seasonal monthly rainfall percentiles to csv
write.csv(P,
          paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "RF percentiles.csv"),
          row.names = F)

########## Add Data To tABLE

Cell.DataCL[u, 5:13] <-
  c(ANNMRFM ,
    AnnMTA,
    AnnMTX,
    AnnMTN,
    ANNMRH,
    SM_P_M,
    KD_P_M,
    ET_P_M,
    CF_P_M)

Cell.DataCLR[1, 3:12] <-
  c(ANNMRFM ,
    AnnMTA,
    AnnMTX,
    AnnMTN,
    ANNMRH,
    SM_P_M,
    KD_P_M,
    ET_P_M,
    CF_P_M,
    WS_P_M)
Cell.DataCLR[2, 3:12] <-
  c(ANNUP,
    AnnMTAx,
    AnnMTXx,
    AnnMTNx,
    ANNMRHx,
    SMUP,
    KDUP,
    ETUP,
    CFUP,
    WSUP)
Cell.DataCLR[3, 3:12] <-
  c(ANNLO,
    AnnMTAn,
    AnnMTXn,
    AnnMTNn,
    ANNMRHn,
    SMLO,
    KDLO,
    ETLO,
    CFLO,
    WSLO)

Cell.RF_Year[u, 3:15] <-  MEANRF2

write.csv(
  Cell.DataCLR,
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "Mean Climate.csv"),
  row.names = F
)
Cell.DataCLR


#######################################################################################################################





##########  Downscaling
print("PDKE: 60, Downscaling")

##########  Dynamical and Statistical Downscaling

print("Downscaling")

PWW_T <- UNIT_X[[u]] #Name Variable UNIT

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

DyRF100UP <-
  max(c(Dy4.5A_x, Dy4.5D_x, Dy4.5W_x, Dy8.5A_x, Dy8.5D_x, Dy8.5W_x))
DyRF100LO <-
  min(c(Dy4.5A_n, Dy4.5D_n, Dy4.5W_n, Dy8.5A_n, Dy8.5D_n, Dy8.5W_n))

########## DyDs Temperature 2100
#Percent Change DyDs RCP 4.5 & 8.5
print("PDKE: 61, Percent Change DyDs RCP 4.5 & 8.5")

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
print("PDKE: 62, Calculate Mean ")

Dy4.5A_MMT <- round(cellStats(Dy4.5A_CT , 'mean'), 1)
Cell.Data_DS[7, 1:2] <- c("Dy4.5A_MMT", Dy4.5A_MMT)
Dy8.5A_MMT <- round(cellStats(Dy8.5A_CT , 'mean'), 1)
Cell.Data_DS[8, 1:2] <- c("Dy8.5A_MMT", Dy8.5A_MMT)

Dy4.5A_x <- round(cellStats(Dy4.5A_CT , 'max'), 1)
Dy8.5A_x <- round(cellStats(Dy8.5A_CT , 'max'), 1)
Dy4.5A_n <- round(cellStats(Dy4.5A_CT , 'min'), 1)
Dy8.5A_n <- round(cellStats(Dy8.5A_CT , 'min'), 1)

##########   StDs Rainfall 2100
print("PDKE: 63, StDs Rainfall 2100")

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

StRF100UP <-
  max(c(St4.5A_x, St4.5D_x, St4.5W_x, St8.5A_x, St8.5D_x, St8.5W_x))
StRF100LO <-
  min(c(St4.5A_n, St4.5D_n, St4.5W_n, St8.5A_n, St8.5D_n, St8.5W_n))


##########   StDs Rainfall 2040-2070
print("PDKE: 64, StDs Rainfall 2040-2070 ")

##########   Percent Change StDs RCP 4.5 & 8.5

StDs_P_4.5_ANN40 <- raster(DyDs_files[49])
StDs_P_4.5_ANN40
St4.5A_M40 <- mask(x = StDs_P_4.5_ANN40 , mask = PWW_T)
St4.5A_C40 <- crop(x = St4.5A_M40, y = extent(PWW_T))

StDs_P_4.5_Dry40 <- raster(DyDs_files[51])
StDs_P_4.5_Dry40
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
print("PDKE: 65, StDs Temperature 2100")

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

DsTAUP <-
  max(Dy4.5A_x,
      Dy8.5A_x,
      Dy4.5A_n,
      Dy8.5A_n,
      St4.5A_x,
      St8.5A_x,
      St4.5A_n,
      St8.5A_n)
DsTALO <-
  min(Dy4.5A_x,
      Dy8.5A_x,
      Dy4.5A_n,
      Dy8.5A_n,
      St4.5A_x,
      St8.5A_x,
      St4.5A_n,
      St8.5A_n)
DsTaup <- max(abs(DsTAUP), abs(DsTALO))
DsTAlo <-  DsTaup * -1

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




########## LULIN ##########################

#Don't do anything with this data yet
########## Percent Change Lulin DyDs 8.5
print("PDKE: 66, Percent Change Lulin DyDs 8.5")

D2_P_8.5_ANN <- raster(D2_files[13])
D28.5A_M <- mask(x = D2_P_8.5_ANN, mask = PWW_T)
D28.5A_C <- crop(x = D28.5A_M, y = extent(PWW_T))

D2_P_8.5_Dry <- raster(D2_files[14])
D2.5D_M <- mask(x = D2_P_8.5_Dry, mask = PWW_T)
D2.5D_C <- crop(x = D2.5D_M, y = extent(PWW_T))

D2_P_8.5_Wet <- raster(D2_files[15])
D28.5W_M <- mask(x = D2_P_8.5_Wet, mask = PWW_T)
D28.5W_C <- crop(x = D28.5W_M, y = extent(PWW_T))

D28.5A_MM <- round(cellStats(Dy8.5A_C , 'mean'), 0)
Cell.Data_DS[25, 1:2] <- c("D28.5A_MM", D28.5A_MM)

D28.5D_MM <- round(cellStats(Dy8.5D_C , 'mean'), 0)
Cell.Data_DS[26, 1:2] <- c("D28.5D_MM", D28.5D_MM)

D28.5W_MM <- round(cellStats(Dy8.5W_C , 'mean'), 0)
Cell.Data_DS[27, 1:2] <- c("D28.5W_MM", D28.5W_MM)

D28.5A_x <- round(cellStats(Dy8.5A_C , 'max'), 0)
D28.5D_x <- round(cellStats(Dy8.5D_C , 'max'), 0)
D28.5W_x <- round(cellStats(Dy8.5W_C , 'max'), 0)

D28.5A_n <- round(cellStats(Dy8.5A_C , 'min'), 0)
D28.5D_n <- round(cellStats(Dy8.5D_C , 'min'), 0)
D28.5W_n <- round(cellStats(Dy8.5W_C , 'min'), 0)

#####################################################

D2RF100UP <- max(c(D28.5A_x, D28.5D_x, D28.5W_x))
D2RF100LO <- min(c(D28.5A_n, D28.5D_n, D28.5W_n))

########## D2Ds Temperature 2100 HERE
print("PDKE: 67, D2Ds Temperature 2100")

#Percent Change DyDs RCP 4.5 & 8.5
#HERE
D2Ds_P_8.5_T_ChgF  <- raster(D2_files[28])
if (TUnit == "\u00B0F") {
  D2Ds_P_8.5_T_ChgF  <- raster(D2_files[31])
}
#crs(D2Ds_P_8.5_T_ChgF)<-crs(EXAMP)
#D2Ds_P_8.5_T_ChgF <- spTransform(D2Ds_P_8.5_T_ChgF, crs(EXAMP))
D28.5A_MT <- mask(x = D2Ds_P_8.5_T_ChgF, mask = PWW_T)
D28.5A_CT <- crop(x = D28.5A_MT, y = extent(PWW_T))
plot(D28.5A_MT)
###########   Calculate Mean

D28.5A_MMT <- round(cellStats(D28.5A_CT , 'mean'), 1)
Cell.Data_DS[28, 1:2] <- c("D28.5A_MMT", D28.5A_MMT)

D28.5A_x <- round(cellStats(D28.5A_CT , 'max'), 1)
D28.5A_n <- round(cellStats(D28.5A_CT , 'min'), 1)


StRF100UP <-
  max(c(St4.5A_x, St4.5D_x, St4.5W_x, St8.5A_x, St8.5D_x, St8.5W_x))
StRF100LO <-
  min(c(St4.5A_n, St4.5D_n, St4.5W_n, St8.5A_n, St8.5D_n, St8.5W_n))
RFdsup <-
  max(c(
    abs(DyRF100UP),
    abs(DyRF100LO),
    abs(StRF100UP),
    abs(StRF100LO),
    abs(D2RF100UP),
    abs(D2RF100LO)
  ))
RFdslo <- RFdsup * -1


##########   Write all DS results

#DY - Dynamical
#ST - Statistical
#A - Annual
#W - Wet season
#D - Dry Season
#8.5 is RCP 8.5
#4.5 is RCP 4.5

write.csv(
  Cell.Data_DS,
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Downscaling.csv"),
  row.names = F
)

##########   FIGURES
print("PDKE: 68, FIGURES")

##########   Dynamical and Statistical Downscaling


############################################################################

#These graphs make use of North arrow and scale bar from RF Maps
colfunc2 <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
BI_brks2 <- round(seq(RFdslo, RFdsup , length = 9), 0)

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " DS_RF_8.5_v2.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

title1 = textGrob("Changes in Annual Rainfall by 2100",
                  gp = gpar(
                    col = "darkred",
                    fontface = "bold",
                    fontsize = 15
                  ))

grid.arrange(
  top = title1,
  spplot(
    Dy4.5A_C,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(
      label = paste0("Dynamical 4.5 (", Dy4.5A_MM , "% change)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1))
  ,
  
  spplot(
    St4.5A_C,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(
      label = paste0("Statistical 4.5 (", St4.5A_MM , "% change)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1))
  ,
  
  spplot(
    Dy8.5A_C,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(
      label = paste0("Dynamical 8.5 (", Dy8.5A_MM , "% change)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1))
  ,
  
  spplot(
    St8.5A_C,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(
      label = paste0("Statistical 8.5 (", St8.5A_MM , "% change)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1))
)



# spplot(Dy8.5D_C, col.regions = colfunc2, equal=FALSE,
#        axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
#        main=list(label=paste0("DyDs Dry (",Dy8.5D_MM ,"%)"),cex=0.9),
#        colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
#        sp.layout = list(UNIT_X[u])) ,

# spplot(St8.5D_C, col.regions = colfunc2, equal=FALSE,
#        axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
#        main=list(label=paste0("StDs Dry (",St8.5D_MM ,"%)"),cex=0.9),
#        colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
#        sp.layout = list(UNIT_X[u])) ,
#
# spplot(Dy8.5W_C, col.regions = colfunc2, equal=FALSE,
#        axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
#        main=list(label=paste0("DyDs Wet (",Dy8.5W_MM ,"%)"),cex=0.9),
#        colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
#        sp.layout = list(UNIT_X[u])) ,
#
# spplot(St8.5W_C, col.regions = colfunc2, equal=FALSE,
#        axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
#        main=list(label=paste0("StDs Wet (",St8.5W_MM ,"%)"),cex=0.9),
#        colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
#        sp.layout = list(UNIT_X[u])))

dev.off()

########## Dynamical and Statistical Downscaling
print("PDKE: 69, Dynamical and Statistical Downscaling ")

########## Downscaling Compare RF RCP 4.5

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " DS_RF_2100_4.5.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

title1 = textGrob(
  "Dynamical and Statistical DS RCP 4.5, 2100",
  gp = gpar(
    col = "darkred",
    fontface = "bold",
    fontsize = 15
  )
)


grid.arrange(
  top = title1,
  spplot(
    Dy4.5A_C,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(label = paste0("DyDs ANN (", Dy4.5A_MM , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)),
  
  spplot(
    St4.5A_C,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(label = paste0("StDs ANN (", Dy4.5A_MM , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Dy4.5D_C,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(label = paste0("DyDs DRY (", Dy4.5D_MM , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    St4.5D_C,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(label = paste0("StDs Dry (", St4.5D_MM , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Dy4.5W_C,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(label = paste0("DyDs WET (", Dy4.5W_MM , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    St4.5W_C,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(label = paste0("StDs WET (", St4.5W_MM , "%)"), cex =
                  0.9),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1))
)

dev.off()

########### Statistical 2040-2070
print("PDKE: 70, Statistical 2040-2070")

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " StDsRF2040.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

title1 = textGrob(
  "Changes in Annual Rainfall by Mid-Century",
  gp = gpar(
    col = "darkred",
    fontface = "bold",
    fontsize = 15
  )
)

grid.arrange(
  top = title1,
  spplot(
    St4.5A_C40,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(
      label = paste0("RCP 4.5 (", St4.5A_MM40 , "% change)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    St8.5A_C40,
    col.regions = colfunc2,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brks2,
    main = list(
      label = paste0("RCP 8.5 (", St8.5A_MM40 , "% change)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)),
  
  ncol = 2
)

# spplot(St4.5D_C40, col.regions = colfunc2, equal=FALSE,
#        axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
#        main=list(label=paste0("StDs DRY (",St4.5D_MM40 ,"%)"),cex=0.9),
#        colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
#        sp.layout = list(UNIT_X[u])) ,
#
# spplot(St8.5D_C40, col.regions = colfunc2, equal=FALSE,
#        axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
#        main=list(label=paste0("StDs DRY (",St8.5D_MM40 ,"%)"),cex=0.9),
#        colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
#        sp.layout = list(UNIT_X[u])) ,
#
# spplot(St4.5W_C40, col.regions = colfunc2, equal=FALSE,
#        axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
#        main=list(label=paste0("StDs WET (",St4.5W_MM40 ,"%)"),cex=0.9),
#        colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
#        sp.layout = list(UNIT_X[u])) ,
#
# spplot(St8.5W_C, col.regions = colfunc2, equal=FALSE,
#        axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
#        main=list(label=paste0("StDs WET (",St8.5W_MM40 ,"%)"),cex=0.9),
#        colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
#        sp.layout = list(UNIT_X[u])))
dev.off()

##########   Downscaling Compare TEMP RCP 8.5 & 4.5 2100
print("PDKE: 71, Downscaling Compare TEMP RCP 8.5 & 4.5 2100")

colfunc3 <- colorRampPalette(brewer.pal(9, "Reds"))(100)
print("PDKE: 71A1")
BI_brksF = 0  # added declaration as next line was throwing "object 'BI_brksF' not found"
if (TUnit == "°C") {
  BI_brksF <- round(seq(0, 7 , length = 7), 0)
}
if (TUnit == "°F") {
  BI_brksF <- round(seq(0, 9 , length = 9), 0)
}
print("PDKE: 71A")
png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " DS_Temp2100.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)
print("PDKE: 71B")

title4 = textGrob("Changes in Air Temperature by 2100",
                  gp = gpar(
                    col = "darkred",
                    fontface = "bold",
                    fontsize = 15
                  ))
print("PDKE: 71C")

# pulled these out of the following loop so I can find the bug, remove these 4 chunks once found

mainList <-
  list(label = paste0("Dynamical 4.5 (", Dy4.5A_MMT, TUnit2, " average)"),
       cex = 0.9)
print("PDKE: 71C1")

colorKeyList <-
  list(space = "right",
       height = 1,
       labels = list(cex = 0.6))
print("PDKE: 71C2")

spLayoutList <- list(UNIT_X[u])
print("PDKE: 71C3")

polylayer <- layer(sp.polygons(SHAPE, lwd = 1))
print("PDKE: 71C4")

print(paste("PDKE: Dy4.5A_CT", Dy4.5A_CT))
print(paste("PDKE: colfunc3", colfunc3))
print(paste("PDKE: BI_brksF", BI_brksF))
print(paste("PDKE: mainList", mainList))
print(paste("PDKE: colorKeyList", colorKeyList))
print(paste("PDKE: spLayoutList", spLayoutList))


spplot(
  Dy4.5A_CT,
  col.regions = colfunc3,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksF,
  main = mainList,
  colorkey = colorKeyList,
  sp.layout = spLayoutList
)
# removed + polylayer

print("PDKE: 71C5")


spplot(
  Dy4.5A_CT,
  col.regions = colfunc3,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksF,
  main = list(
    label = paste0("Dynamical 4.5 (", Dy4.5A_MMT , TUnit2, " average)"),
    cex = 0.9
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 0.6)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
print("PDKE: 71D")

spplot(
  St4.5A_CT,
  col.regions = colfunc3,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksF,
  main = list(
    label = paste0("Statistical 4.5 (", St4.5A_MMT , TUnit2, " average)"),
    cex = 0.9
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 0.6)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
print("PDKE: 71E")


spplot(
  Dy8.5A_CT,
  col.regions = colfunc3,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksF,
  main = list(
    label = paste0("Dynamical 8.5 (", Dy8.5A_MMT , TUnit2, " average)"),
    cex = 0.9
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 0.6)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
print("PDKE: 71F")

spplot(
  St8.5A_CT,
  col.regions = colfunc3,
  equal = FALSE,
  axes = TRUE,
  las = 1,
  cex.axis = 0.7,
  at = BI_brksF,
  main = list(
    label = paste0("Statistical 8.5 (", St8.5A_MMT , TUnit2, " average)"),
    cex = 0.9
  ),
  colorkey = list(
    space = "right",
    height = 1,
    labels = list(cex = 0.6)
  ),
  sp.layout = list(UNIT_X[u])
) +
  layer(sp.polygons(SHAPE, lwd = 1))
print("PDKE: 71G")



grid.arrange(
  top = title4,
  spplot(
    Dy4.5A_CT,
    col.regions = colfunc3,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksF,
    main = list(
      label = paste0("Dynamical 4.5 (", Dy4.5A_MMT , TUnit2, " average)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    St4.5A_CT,
    col.regions = colfunc3,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksF,
    main = list(
      label = paste0("Statistical 4.5 (", St4.5A_MMT , TUnit2, " average)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    Dy8.5A_CT,
    col.regions = colfunc3,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksF,
    main = list(
      label = paste0("Dynamical 8.5 (", Dy8.5A_MMT , TUnit2, " average)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    St8.5A_CT,
    col.regions = colfunc3,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksF,
    main = list(
      label = paste0("Statistical 8.5 (", St8.5A_MMT , TUnit2, " average)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1))
)
print("PDKE: 71D")

dev.off()


######################## AIR TEMP 2040-2070
print("PDKE: 72, AIR TEMP 2040-2070 ")

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " StDs_Temp2040.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

title4 = textGrob(
  "Changes in Air Temperature by Mid-Century",
  gp = gpar(
    col = "darkred",
    fontface = "bold",
    fontsize = 15
  )
)

grid.arrange(
  top = title4,
  spplot(
    St4.5A_CT40,
    col.regions = colfunc3,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksF,
    main = list(
      label = paste0("RCP 4.5 (", St4.5A_MMT40 , TUnit2, " change)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)) ,
  
  spplot(
    St8.5A_CT40,
    col.regions = colfunc3,
    equal = FALSE,
    axes = TRUE,
    las = 1,
    cex.axis = 0.7,
    at = BI_brksF,
    main = list(
      label = paste0("RCP 8.5 (", St8.5A_MMT40 , TUnit2, " change)"),
      cex = 0.9
    ),
    colorkey = list(
      space = "right",
      height = 1,
      labels = list(cex = 0.6)
    ),
    sp.layout = list(UNIT_X[u])
  ) +
    layer(sp.polygons(SHAPE, lwd = 1)),
  ncol = 2
)

dev.off()

##########################################################################################################################




##########  Create A Monthly Rainfall Time Series

print("Rainfall Extract")
print("Frazier et al 2016")

#########   Extract RF Data From Monthly Maps (2 DATASETS)
#########   Call in the Shapefile
PWW_T <- UNIT_X[[u]]

#########   Load RF MAPS
RF_Map_Path_A
RF_Tif_files = dir(
  RF_Map_Path_A,
  pattern = "*.tif",
  recursive = T,
  full.names = T
)  #Monthly RF
nfiles <- length(RF_Tif_files)

#Create a matrix for each cell
Cell.AF_Maps <- data.frame(matrix(nrow = 1116, ncol = 4))
colnames(Cell.AF_Maps) <- c("Date", "Year", "Month", "RF")

for (i in 1:nfiles)        {
  if (i == 1) {
    print("AF")
  }
  if (i == 250) {
    print(i)
  } #Just to see where Im at in the process
  if (i == 500) {
    print(i)
  }
  if (i == 750) {
    print(i)
  }
  if (i == 1000) {
    print(i)
  }
  
  #Open Map as raster and change projection
  RF_Map <- RF_Tif_files[i]
  RF_Map2 <- raster(RF_Map)
  #crs(RF_Map2) <- " +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  plot(RF_Map2)
  
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

print("Lucas")


#Load Daily RF MAPS
RF_Map_Path
RF_Tif_files = dir(
  RF_Map_Path,
  pattern = "*.tif",
  recursive = T,
  full.names = T
)  #Monthly RF
nfiles <- length(RF_Tif_files)

#Create a matrix for each cell
Cell.ML_Maps <- data.frame(matrix(nrow = 396, ncol = 4))
colnames(Cell.ML_Maps) <- c("Date", "Year", "Month", "RF")

for (i in 1:nfiles) {
  if (i == 1) {
    print("ML")
  }
  if (i == 100) {
    print(i)
  } #Just to see where Im at in the process
  if (i == 200) {
    print(i)
  }
  if (i == 300) {
    print(i)
  }
  
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
if (RFUnit == " in") {
  MRF_A3 <-  MRF_A3 * 0.0393701
}
summary(MRF_A3)

### Matty's data
# subet for 1990 - 2012
MRF_N2 =  Cell.ML_Maps[c(1:276), ]
head(MRF_N2)

# just keep RF values
MRF_N3 =  as.numeric(MRF_N2[, 4])

# convert mm to inches
if (RFUnit == " in") {
  MRF_N3  <-  MRF_N3  * 0.0393701
}
summary(MRF_N3)



MBE <- round(mean(MRF_A3 - MRF_N3), 1)
MAE <- round(mean(abs(MRF_A3 - MRF_N3)), 1)

D_Comp <- cbind(MRF_A3, MRF_N3)
Mx <- max(D_Comp, na.rm = T)

FNAME <- paste0("RF_Compare_23_", UNIT_Ns[u], ".csv")

##########   Comparison Figures For the two datasets

LM1 <- lm(MRF_A3 ~ MRF_N3)
LM1P <- round(coefficients(summary(LM1))[2, 4], 4)
LM1R <- round(summary(LM1)$r.squared, 2)
TITLE = paste(UNIT_Ns[u], " 23-yr RF Compare (", RFUnit2, ")")

dpi = 300

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " 23yr_RF_Compare.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

plot(
  MRF_A3 ~ MRF_N3,
  ylim = c(0, Mx),
  xlim = c(0, Mx),
  main = TITLE,
  ylab = "Frazier et al. (2016)",
  xlab = "Lucas et al. (In Review)"
)
abline(0, 1)
legend("topleft", c(paste("R2 = ", LM1R), paste("MBE = ", MBE), paste("MAE = ", MAE)))

dev.off()

##########   Merge Datasets to create full time period (1920 - current)
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
if (RFUnit == " in") {
  MRF100$RF <-  as.numeric(MRF100$RF) * 0.0393701
}

write.csv(
  MRF100,
  paste0(
    OUTPUTS_FOLDER,
    UNIT_N[u],
    "/",
    UNIT_N[u],
    " Monthly Rainfall_",
    RFUnit2,
    ".csv"
  ),
  row.names = F
)

# ## Can read in the csv from above to start script here
# MRF100<-read.csv(paste0(OUTPUTS_FOLDER,UNIT_N[u],"/",UNIT_N[u]," Monthly Rainfall_in.csv"))

# fix month values
MRF100$Month <- sub(".*/", "", MRF100$Date)

# get end date
ey <- as.numeric(MRF100[nrow(MRF100), ]$Year)
em <- as.numeric(MRF100[nrow(MRF100), ]$Month)
ed <- as.Date(MRF100[nrow(MRF100), ]$Ndate)

head(MRF100)
tail(MRF100)

ggplot(data = MRF100, aes(x = Year, y = RF)) +
  geom_line()

RF_IN <- as.numeric(MRF100[, 4])
summary(RF_IN)

MRF100$Day <- 1
MRF100$Ndate <-
  as.Date(with(MRF100, paste(Year, Month, Day, sep = "-")), "%Y-%m-%d")

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
if (coef(LM1)[[2]] > 0) {
  LM1s <- c("Increase")
}
if (coef(LM4)[[2]] > 0) {
  LM4s <- c("Increase")
}
if (coef(LM6)[[2]] > 0) {
  LM6s <- c("Increase")
}
if (coef(LM1)[[2]] < 0) {
  LM1s <- c("Decrease")
}
if (coef(LM4)[[2]] < 0) {
  LM4s <- c("Decrease")
}
if (coef(LM6)[[2]] < 0) {
  LM6s <- c("Decrease")
}

l <- paste0("1920 - ", ey)
m <- paste0("1980 - ", ey)
s <- paste0("2010 - ", ey)

RFT <- data.frame(Period = c(l, m, s))
RFT
RFT$Trend <- c(LM1s, LM4s, LM6s)
RFT
write.csv(
  RFT,
  paste0(
    OUTPUTS_FOLDER,
    UNIT_N[u],
    "/",
    UNIT_N[u],
    " RF_trend_dirctions.csv"
  ),
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
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Annual_RF_in.csv"),
  row.names = F
)


##########   Plot Annual RF

YrRF.ts <- ts(Mean_Y_RF$RF,
              c(1920),
              end = c(ey),
              frequency = 1)
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
WET_RF5 <-
  as.vector(tapply(WET_RF4, gl(length(WET_RF4) / 6, 6), mean, na.rm = T))

# WET_RF<-MRF100[MRF100$Month !="5" & MRF100$Month !="6" & MRF100$Month != "7" &
#                  MRF100$Month != "8" & MRF100$Month != "9" & MRF100$Month != "10",]
# WET_RF5<-as.numeric(WET_RF$RF)

# Dry Season
DRY_RF <- MRF100[MRF100a$Month != 1 & MRF100a$Month  != 02 &
                   MRF100a$Month  != 3 & MRF100a$Month  != 04 &
                   MRF100a$Month  != 11 & MRF100a$Month  != 12, ]
head(DRY_RF, 12)
DRY_RF2 <- as.numeric(DRY_RF$RF)

#Get seasonal average
DRY_RF3 <-
  as.vector(tapply(DRY_RF2, gl(length(DRY_RF2) / 6, 6), mean, na.rm = T))

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
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " RF_Trend.png"),
  width = 6.3 * dpi,
  height = 7 * dpi,
  res = dpi
)

##########   Annual and Seasonal Plot

par(mfrow = c(3, 1))
par(mar = c(4, 4, 4, 2))

MLIM1 <- max(c(myts1Y, myts1YW, myts1YD), na.rm = T)
YLIM <-  min(myts1Y)
MLIM <-  (MLIM1 + (MLIM1 * 0.45))

par(mai = c(0.3, 0.6, 0.2, 0.2))
plot(
  myts1Y ~ YDateT1,
  ylab = paste0("Average Rainfall (", RFUnit2, "/month)"),
  type = "l",
  col = "blue",
  xlab = "",
  xaxt = "n",
  ylim = c(YLIM, MLIM),
  cex.axis = 1.3,
  las = 1
)
title(
  paste("Rainfall Trend 1920-", ey, ":", UNIT_Ns[u]),
  line = 0.5,
  cex.main = 1.5
)

legend(
  "topright",
  c(
    paste0("1920-", ey, " R2 = ", LM1RY, " p = ", LM1PY),
    #paste0("1940-",ey," Trend =", T2Y,", R2 =",LM2RY, " p = ",LM2PY),
    #paste0("1960-",ey," R2 = ",LM3RY, " p = ",LM3PY),
    paste0("1980-", ey, " R2 = ", LM4RY, " p = ", LM4PY),
    #paste0("2000-",ey," R2 = ",LM5RY, " p = ",LM5PY)),
    paste0("2010-", ey, " R2 = ", LM6RY, " p = ", LM6PY)
  ),
  #lty = 1, col = c("darkred","darkorange","darkgreen","darkblue","purple","darkcyan"), lwd = 3)
  lty = 1,
  col = c("grey70", "grey30", "grey1"),
  lwd = 3
)

legend(
  "topleft",
  c("Annual"),
  cex = 1.5,
  text.font = 2,
  bty = "n"
)

ablineclip(
  lm(myts1Y ~ YDateT1),
  x1 = -19000,
  x2 = 20000,
  col = "grey70",
  lwd = 3
)
#ablineclip(lm(myts2Y~YDateT2),x1=-11000,x2=19000,col="darkorange",lwd=3)
#ablineclip(lm(myts3Y~YDateT3),x1=-4500,x2=19000,col=alpha("grey30",0.7),lwd=3)
ablineclip(
  lm(myts4Y ~ YDateT4),
  x1 = 3500,
  x2 = 20000,
  col = "grey30",
  lwd = 3
)
#ablineclip(lm(myts5Y~YDateT5),x1=11000,x2=19000,col=alpha("grey1",0.7),lwd=3)
ablineclip(
  lm(myts6Y ~ YDateT6),
  x1 = 14000,
  x2 = 20000,
  col = "grey1",
  lwd = 3
)

####### Wet Season ###########
par(mai = c(0.3, 0.6, 0.2, 0.2))
YLIM <-  min(myts1YW, na.rm = T)
plot(
  myts1YW ~ YDateT1,
  ylab = paste0("Average Rainfall (", RFUnit2, "/month)"),
  type = "l",
  col = "blue",
  xlab = "",
  xaxt = "n",
  ylim = c(YLIM, MLIM),
  cex.axis = 1.3,
  las = 1
)

legend(
  "topright",
  c(
    paste0("1920-", ey, " R2 = ", LM1RYW, " p = ", LM1PYW),
    #paste0("1940-",ey," Trend =", T2YW,", R2 =",LM2RYW, " p = ",LM2PYW),
    #paste0("1960-",ey," R2 = ",LM3RYW, " p = ",LM3PYW),
    paste0("1980-", ey, " R2 = ", LM4RYW, " p = ", LM4PYW),
    #paste0("2000-",ey," R2 = ",LM5RYW, " p = ",LM5PYW)),
    paste0("2010-", ey, " R2 = ", LM6RYW, " p = ", LM6PYW)
  ),
  #lty = 1, col = c("darkred","darkorange","darkgreen","darkblue","purple","darkcyan"), lwd = 3,cex=1)
  lty = 1,
  col = c("grey70", "grey30", "grey1"),
  lwd = 3
)

legend(
  "topleft",
  c("Wet Season"),
  cex = 1.5,
  text.font = 2,
  bty = "n"
)

ablineclip(
  lm(myts1YW ~ YDateT1),
  x1 = -19000,
  x2 = 20000,
  col = "grey70",
  lwd = 3
)
#ablineclip(lm(myts2YW~YDateT2),x1=-11000,x2=20000,col="darkorange",lwd=3)
#ablineclip(lm(myts3YW~YDateT3),x1=-4500,x2=20000,col="grey30",lwd=3)
ablineclip(
  lm(myts4YW ~ YDateT4),
  x1 = 3500,
  x2 = 20000,
  col = "grey30",
  lwd = 3
)
#ablineclip(lm(myts5YW~YDateT5),x1=11000,x2=20000,col="grey1",lwd=3)
ablineclip(
  lm(myts6YW ~ YDateT6),
  x1 = 14000,
  x2 = 20000,
  col = "grey1",
  lwd = 3
)

########## Dry Season

par(mai = c(0.3, 0.6, 0.2, 0.2))
YLIM <-  min(myts1YD, na.rm = T)
plot(
  myts1YD ~ YDateT1,
  ylab = paste0("Average Rainfall (", RFUnit2, "/month)"),
  type = "l",
  col = "blue",
  xlab = "",
  ylim = c(YLIM, MLIM),
  cex.axis = 1.3,
  las = 1
)
# axis(1, labels = T)
# title(main = "Average Wet Season Rainfall Pu'u Wa'awa'a (1920-2021)", line = 1)
legend(
  "topright",
  c(
    paste0("1920-", ey, " R2 = ", LM1RYD, " p = ", LM1PYD),
    #paste0("1940-",ey," Trend =", T2YD,", R2 =",LM2RYD, " p = ",LM2PYD),
    #paste0("1960-",ey," R2 = ",LM3RYD, " p = ",LM3PYD),
    paste0("1980-", ey, " R2 = ", LM4RYD, " p = ", LM4PYD),
    #paste0("2000-",ey," R2 = ",LM5RYD, " p = ",LM5PYD)),
    paste0("2010-", ey, " R2 = ", LM6RYD, " p = ", LM6PYD)
  ),
  #lty = 1, col = c("darkred","darkorange","darkgreen","darkblue","purple","darkcyan"), lwd = 3)
  lty = 1,
  col = c("grey70", "grey30", "grey1"),
  lwd = 3
)

legend(
  "topleft",
  c("Dry Season"),
  cex = 1.5,
  text.font = 2,
  bty = "n"
)

ablineclip(
  lm(myts1YD ~ YDateT1),
  x1 = -20000,
  x2 = 20000,
  col = "grey70",
  lwd = 3
)
#ablineclip(lm(myts2YD~YDateT2),x1=-11000,x2=20000,col="darkorange",lwd=3)
#ablineclip(lm(myts3YD~YDateT3),x1=-4500,x2=20000,col="grey30",lwd=3)
ablineclip(
  lm(myts4YD ~ YDateT4),
  x1 = 3500,
  x2 = 20000,
  col = "grey30",
  lwd = 3
)
#ablineclip(lm(myts5YD~YDateT5),x1=11000,x2=20000,col="grey1",lwd=3)
ablineclip(
  lm(myts6YD ~ YDateT6),
  x1 = 14000,
  x2 = 20000,
  col = "grey1",
  lwd = 3
)

dev.off()


### SPI/DROUGHT ###

#############################################################################
# Count All Drought Events Long-term 12 and Short term 12 and 3.
Cell.SPICNT <- data.frame(matrix(ncol = 4, nrow = 4))
colnames(Cell.SPICNT) <-
  c("Total", "SPI_12_Long", "SPI_3_Short", "SPI_12_Short")
Cell.SPICNT[1:4, 1] <-
  c("Drought Events", "Moderate", "Severe", "Extreme")

### Ryan's code
print("SPI FIG 12")
RF <- as.numeric(MRF100$RF)
RFDATA <- cbind(MRF100[, c(2, 3)])

#SPI-12 Calculation
SPI12 <- spi(RF, scale = 12, distribution = 'Gamma')
#SPI-3 Calculation
SPI3 <- spi(RF, scale = 3, distribution = 'Gamma')

#PlOT ALL SPI
png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " SPI.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

plot(spi(
  ts(RF, freq = 12, start = c(1920, 1)),
  scale = 12,
  distribution = 'Gamma'
),
main = paste0("SPI-12 1920-", ey, ": ", UNIT_Ns[u]))

dev.off()
##########   Table Metrics SPI 12

Cell.DataSPI <- data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI) <-
  c("Start",
    "End",
    "Duration",
    "A Intensity",
    "P Intensity",
    "Magnitude")
Cell.DataSPI

##### Derek's Drought Code from Guam
# load rainfall dataset created around line 2828 above
wd <- paste0(OUTPUTS_FOLDER, UNIT_N, "/")
setwd(wd)
csv <- paste0(UNIT_N[u], " Monthly Rainfall_", RFUnit2, ".csv")
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

spi3_val <- spi3$fitted
spi3_val <-
  data.frame(spi = as.matrix(spi3_val), date = time(spi3_val))
spi3_val$m.scale <- 3
spi3_val$date <- sub("\\..*", "", spi3_val$date)

head(spi3_val, 20)
tail(spi3_val)

# combine dataframes
SPI_ALL <- rbind(spi_val, spi3_val)
head(SPI_ALL, 20)
tail(SPI_ALL, 20)

# determine if last year is complete, or how many months are missing
# complete will be 24 months in 1 year because two time scales
ly <- nrow(SPI_ALL[which(SPI_ALL$date == max(SPI_ALL$date)), ])
m <- 24 - ly

# make extra rows so the last year is complete
extra <- data.frame(matrix(ncol = 3, nrow = m))
colnames(extra) <- c("spi", "date", "m.scale")
extra$date <- max(SPI_ALL$date)
extra[1:(m / 2), ]$m.scale <- 3
extra[((m / 2) + 1):nrow(extra), ]$m.scale <- 12

extra
head(SPI_ALL)
# add extra rows to dataset
SPI_ALL <- rbind(SPI_ALL, extra)

# sort by year
SPI_ALL <- SPI_ALL[order(SPI_ALL$date, SPI_ALL$m.scale), ]

# fix date column
SPI_ALL$date2 <- NA
SPI_ALL$date2 <- rep(c(1:12))
SPI_ALL$date <-
  as.Date(paste0(SPI_ALL$date, "/", SPI_ALL$date2, "/01"), format = "%Y/%m/%d")
head(SPI_ALL, 20)
tail(SPI_ALL, 20)

SPI_ALL <- SPI_ALL[1:3]
colnames(SPI_ALL) <- c("SPI", "date", "m.scale")

# sort by date
SPI_ALL <- SPI_ALL[order(as.Date(SPI_ALL$date)), ]

# Convert drought events (negative values) to positive
SPI_ALL$spi_negs <- ifelse(SPI_ALL$SPI > 0, 0, SPI_ALL$SPI)
SPI_ALL$spi_negs <- abs(SPI_ALL$spi_negs)

summary(SPI_ALL$spi_negs)
head(SPI_ALL, 50)

# Save SPI_ALL drought intensity (inverted SPI) dataset
write.csv(
  SPI_ALL,
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " SPI_NEGS_ALL.csv"),
  row.names = F
)

### Drought Event Count for SPI-12 ###

spi2 <- SPI_ALL[which(SPI_ALL$m.scale == 12), ]

# create binary column of drought (SPI>0) yes or no
spi2$drought <- ifelse(spi2$spi_negs > 0, 1, 0)
head(spi2, 50)

# make column of consecutive drought months
spi2$DRGHT_ct <- with(spi2, (drought == 1) *
                        ave(drought, rleid(drought == 1),
                            FUN = seq_along))

summary(spi2$DRGHT_ct)

# subset data for drought event months only
spi3 <- spi2[which(spi2$DRGHT_ct >= 1), ]
head(spi3, 20)

# create empty event count column
spi3$event_ct <- 0

### fill in event_ct column
for (r in 1:nrow(spi3)) {
  SPI_M <- spi3[r, ]
  SPI_M
  
  if (r == 1) {
    spi3[r, ]$event_ct <- r
  }
  if (SPI_M$DRGHT_ct > 1) {
    spi3[r, ]$event_ct <- spi3[r - 1, ]$event_ct
  }
  if (SPI_M$DRGHT_ct == 1 &&
      r > 1) {
    spi3[r, ]$event_ct <- spi3[r - 1, ]$event_ct + 1
  }
}

head(spi3, 100)
summary(spi3$event_ct)

### keep only events with peak spi_negs >=1
event <- unique(spi3$event_ct)

for (x in event) {
  sub <- spi3[which(spi3$event_ct == x), ]
  
  if (max(sub$spi_negs) < 1) {
    spi3[which(spi3$event_ct == x), ]$event_ct <- 0
  }
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
spi5 <- merge(spi2, spi4, all = T)

head(spi5, 10)
summary(spi5$event_ct2)

# get rid of wrong event column
spi5 <- subset(spi5, select = -c(event_ct))

colnames(spi5)[which(names(spi5) == "event_ct2")] <-
  "event_ct"

# sort by date
spi5 <- spi5[order(as.Date(spi5$date)), ]
head(spi5, 20)

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
  
  if (max(SPI_I$spi_negs) <= 1.5) {
    SPI_I$intensity <- 1
  }
  if (max(SPI_I$spi_negs) > 1.5 &&
      max(SPI_I$spi_negs) <= 2) {
    SPI_I$intensity <- 2
  }
  if (max(SPI_I$spi_negs) > 2) {
    SPI_I$intensity <- 3
  }
  
  SPI_I$peak <- max(SPI_I$spi_negs)
  SPI_I$mean <- mean(SPI_I$spi_negs)
  SPI_I$mag <- sum(SPI_I$spi_negs)
  
  # keep only these columns
  SPI_I <-
    SPI_I[c("date", "months", "intensity", "peak", "mean", "mag")]
  SPI_I
  
  spi6 <- merge(spi6, SPI_I, by = "date", all.x = T)
  
  # merge intensity columns together
  if (x > 1) {
    spi6$months <-
      ifelse(is.na(spi6$months.x), spi6$months.y, spi6$months.x)
  }
  if (x > 1) {
    spi6$intensity <-
      ifelse(is.na(spi6$intensity.x),
             spi6$intensity.y,
             spi6$intensity.x)
  }
  if (x > 1) {
    spi6$peak <- ifelse(is.na(spi6$peak.x), spi6$peak.y, spi6$peak.x)
  }
  if (x > 1) {
    spi6$mean <- ifelse(is.na(spi6$mean.x), spi6$mean.y, spi6$mean.x)
  }
  if (x > 1) {
    spi6$mag <- ifelse(is.na(spi6$mag.x), spi6$mag.y, spi6$mag.x)
  }
  
  if (x > 1) {
    spi6 <- subset(
      spi6,
      select = -c(
        months.x,
        months.y,
        intensity.x,
        intensity.y,
        peak.x,
        peak.y,
        mean.x,
        mean.y,
        mag.x,
        mag.y
      )
    )
  }
}

head(spi6)
head(spi6, 100)
summary(spi6$intensity)
summary(spi6$peak)
summary(spi6$months)

### make table of start and end dates for each drought intensity event

Cell.DataSPI <- data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI) <-
  c("Start",
    "End",
    "Duration",
    "A Intensity",
    "P Intensity",
    "Magnitude")
Cell.DataSPI$Start <- as.Date(Cell.DataSPI$Start)
Cell.DataSPI$End <- as.Date(Cell.DataSPI$End)
Cell.DataSPI
colnames(Cell.DataSPI)

### loop through rows and fill table

for (y in 1:nrow(spi6)) {
  # get row
  row <- spi6[y, ]
  
  # get previous and next row
  prev <- spi6[y - 1, ]
  nex <- spi6[y + 1, ]
  
  # get event number
  number <- row$event_ct
  
  if (!is.na(row$intensity) &&
      is.na(prev$intensity)) {
    Cell.DataSPI[number, ]$Start <- as.Date(nex$date)
  }
  #if(!is.na(row$intensity) && is.na(prev$intensity)) {Cell.DataSPI[number,]$intensity<-row$intensity}
  if (!is.na(row$intensity) &&
      is.na(nex$intensity)) {
    Cell.DataSPI[number, ]$End <- as.Date(nex$date)
  }
  if (!is.na(row$intensity)) {
    Cell.DataSPI[number, ]$Duration <- row$months
  }
  #if(!is.na(row$intensity)) {Cell.DataSPI[number,]$event_ct<-row$event_ct}
  if (!is.na(row$intensity)) {
    Cell.DataSPI[number, ]$`P Intensity` <- row$peak
  }
  if (!is.na(row$intensity)) {
    Cell.DataSPI[number, ]$`A Intensity` <- row$mean
  }
  if (!is.na(row$intensity)) {
    Cell.DataSPI[number, ]$Magnitude <- row$mag
  }
}

# fix date formats
Cell.DataSPI$Start2 <-
  as.yearmon(sub(" .*", "", Cell.DataSPI$Start), "%Y-%m-%d")
Cell.DataSPI$End2 <-
  as.yearmon(sub(" .*", "", Cell.DataSPI$End), "%Y-%m-%d")

Cell.DataSPI

###### Back to Ryan's code #####

write.csv(
  Cell.DataSPI,
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Drought History.csv"),
  row.names = F
)

##########   Count Droughts For Figure - Derek's edits
Cell.DataSPI$`P Intensity` <- as.numeric(Cell.DataSPI$`P Intensity`)

EX_Cnt <- sum(Cell.DataSPI[5] > 2)
SV_Cnt <- sum(Cell.DataSPI[5] > 1.5)
MO_Cnt <- sum(Cell.DataSPI[5] > 1)

D_Cnt <- MO_Cnt
SV_Cnt2 <- SV_Cnt - EX_Cnt
MO_Cnt2 <- D_Cnt -  SV_Cnt2 - EX_Cnt

# create SPI count dataframe
Cell.SPICNT <- data.frame(matrix(ncol = 4, nrow = 4))
colnames(Cell.SPICNT) <-
  c("Total", "SPI_12_Long", "SPI_3_Short", "SPI_12_Short")
Cell.SPICNT[1:4, 1] <-
  c("Drought Events", "Moderate", "Severe", "Extreme")

Cell.SPICNT[1:4, 2] <- c(D_Cnt, MO_Cnt2, SV_Cnt2, EX_Cnt)
Cell.SPICNT

##########   Remove Postive SPI and make absloute values
### Derek's code
SPIVEC <- SPI_ALL[which(SPI_ALL$m.scale == 12), ]$SPI
# SPIVEC[SPIVEC > 0] <- 0
# SPIVEC_Abs <- as.vector(abs(SPIVEC))

M_SPI.ts <- ts(SPIVEC,
               c(1920, 1),
               end = c(ey, em),
               frequency = 12)
myts66 <-
  as.vector(window(M_SPI.ts, start = c(1920, 1), end = c(ey, em)))
DateT1 <-
  as.Date(seq(as.Date("1920-01-01"), as.Date(ed), by = "months"))
#short.date_M = strftime(MonthlyRF$Ndate, "%Y-%m")

xx <- data.frame(DateT1, myts66)
colnames(xx) <- c("DT", "SP")
xx$DT <- as.Date(xx$DT)

write.csv(xx,
          paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "SPI_12.csv"),
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

write.csv(spi12am,
          paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "SPI_12 monthly.csv"))

# change positive SPI values to 0 and make drought figure
SPIVEC <- SPI_ALL[which(SPI_ALL$m.scale == 12), ]$SPI
SPIVEC[SPIVEC > 0] <- 0
SPIVEC_Abs <- as.vector(abs(SPIVEC))

M_SPI.ts <-
  ts(SPIVEC_Abs,
     c(1920, 1),
     end = c(ey, em),
     frequency = 12)
myts66 <-
  as.vector(window(M_SPI.ts, start = c(1920, 1), end = c(ey, em)))
DateT1 <-
  as.Date(seq(as.Date("1920-01-01"), as.Date(ed), by = "months"))
#short.date_M = strftime(MonthlyRF$Ndate, "%Y-%m")

xx <- data.frame(DateT1, myts66)
colnames(xx) <- c("DT", "SP")
xx$DT <- as.Date(xx$DT)
xx
tail(xx)

dpi = 300
png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "Drought_History.png"),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

print

(
  ggplot(xx, aes(x = DT, y = SP)) +
    geom_area(fill = "darkorange", color = "black") +
    xlab("") +
    
    labs(
      title = paste0("SPI-12 Drought Events 1920 -", ey, ": ", UNIT_Ns[u]),
      x = "",
      y = "Drought Intensity"
    ) +
    geom_hline(
      yintercept = 2,
      linetype = "dashed",
      color = "darkred",
      size = 1
    ) +
    geom_hline(
      yintercept = 1.5,
      linetype = "dashed",
      color = "red",
      size = 1
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "orange",
      size = 1
    )  +
    annotate(
      "text",
      x = xx$DT[350],
      y = 3.4,
      label = paste0("Drought Events = ",  D_Cnt),
      size = 5,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = xx$DT[350],
      y = 3.1,
      label = paste0("Moderate Droughts = ", MO_Cnt2),
      size = 5,
      fontface = "bold",
      colour = "orange"
    ) +
    annotate(
      "text",
      x = xx$DT[350],
      y = 2.8,
      label = paste0("Severe Droughts = ", SV_Cnt2),
      size = 5,
      fontface = "bold",
      colour = "red"
    ) +
    annotate(
      "text",
      x = xx$DT[350],
      y = 2.5,
      label = paste0("Extreme Droughts = ", EX_Cnt),
      size = 5,
      fontface = "bold",
      colour = "darkred"
    )
)


dev.off()





##########################################################################################################################################

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
spi2$DRGHT_ct <- with(spi2, (drought == 1) *
                        ave(drought, rleid(drought == 1),
                            FUN = seq_along))

# subset data for drought event months only
spi3 <- spi2[which(spi2$drought >= 1), ]
head(spi3, 20)

# create empty event count column
spi3$event_ct <- 0

### fill in event_ct column
for (r in 1:nrow(spi3)) {
  SPI_M <- spi3[r, ]
  SPI_M
  
  if (r == 1) {
    spi3[r, ]$event_ct <- r
  }
  if (SPI_M$DRGHT_ct > 1) {
    spi3[r, ]$event_ct <- spi3[r - 1, ]$event_ct
  }
  if (SPI_M$DRGHT_ct == 1 &&
      r > 1) {
    spi3[r, ]$event_ct <- spi3[r - 1, ]$event_ct + 1
  }
}

head(spi3, 50)
summary(spi3$event_ct)

### keep only events with peak spi_negs >=1
event <- unique(spi3$event_ct)

for (x in event) {
  sub <- spi3[which(spi3$event_ct == x), ]
  
  if (max(sub$spi_negs) < 1) {
    spi3[which(spi3$event_ct == x), ]$event_ct <- 0
  }
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
spi5 <- merge(spi2, spi4, all = T)

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
  
  if (max(SPI_I$spi_negs) <= 1.5) {
    SPI_I$intensity <- 1
  }
  if (max(SPI_I$spi_negs) > 1.5 &&
      max(SPI_I$spi_negs) <= 2) {
    SPI_I$intensity <- 2
  }
  if (max(SPI_I$spi_negs) > 2) {
    SPI_I$intensity <- 3
  }
  
  SPI_I$peak <- max(SPI_I$spi_negs)
  SPI_I$mean <- mean(SPI_I$spi_negs)
  SPI_I$mag <- sum(SPI_I$spi_negs)
  
  # keep only these columns
  SPI_I <- SPI_I[c("date", "months", "intensity", "peak", "mean", "mag")]
  SPI_I
  
  spi6 <- merge(spi6, SPI_I, by = "date", all.x = T)
  #
  # # sort by date
  # spi5<-spi5[order(as.Date(spi5$date)),]
  # head(spi5, 100)
  
  # merge intensity columns together
  if (x > 1) {
    spi6$months <-
      ifelse(is.na(spi6$months.x), spi6$months.y, spi6$months.x)
  }
  if (x > 1) {
    spi6$intensity <-
      ifelse(is.na(spi6$intensity.x),
             spi6$intensity.y,
             spi6$intensity.x)
  }
  if (x > 1) {
    spi6$peak <- ifelse(is.na(spi6$peak.x), spi6$peak.y, spi6$peak.x)
  }
  if (x > 1) {
    spi6$mean <- ifelse(is.na(spi6$mean.x), spi6$mean.y, spi6$mean.x)
  }
  if (x > 1) {
    spi6$mag <- ifelse(is.na(spi6$mag.x), spi6$mag.y, spi6$mag.x)
  }
  
  if (x > 1) {
    spi6 <- subset(
      spi6,
      select = -c(
        months.x,
        months.y,
        intensity.x,
        intensity.y,
        peak.x,
        peak.y,
        mean.x,
        mean.y,
        mag.x,
        mag.y
      )
    )
  }
}

head(spi6, 100)
summary(spi6$intensity)
summary(spi6$peak)
summary(spi6$months)

### make table of start and end dates for each drought intensity event

Cell.DataSPI3_S <- data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI3_S) <-
  c("Start",
    "End",
    "Duration",
    "A Intensity",
    "P Intensity",
    "Magnitude")

Cell.DataSPI3_S$Start <- as.POSIXct(Cell.DataSPI3_S$Start)
Cell.DataSPI3_S$End <- as.POSIXct(Cell.DataSPI3_S$End)
Cell.DataSPI3_S
colnames(Cell.DataSPI3_S)

### loop through rows and fill table

for (y in 1:nrow(spi6)) {
  # get row
  row <- spi6[y, ]
  row
  
  # get previous and next row
  prev <- spi6[y - 1, ]
  
  nex <- spi6[y + 1, ]
  
  # get event number
  number <- row$event_ct
  
  if (!is.na(row$intensity) &&
      is.na(prev$intensity)) {
    Cell.DataSPI3_S[number, ]$Start <- as.Date(nex$date)
  }
  # if(!is.na(row$intensity) && is.na(prev$intensity)) {event.t[number,]$intensity<-row$intensity}
  if (!is.na(row$intensity) &&
      is.na(nex$intensity)) {
    Cell.DataSPI3_S[number, ]$End <- as.Date(nex$date)
  }
  if (!is.na(row$intensity)) {
    Cell.DataSPI3_S[number, ]$Duration <- row$months
  }
  # if(!is.na(row$intensity)) {event.t[number,]$event_ct<-row$event_ct}
  if (!is.na(row$intensity)) {
    Cell.DataSPI3_S[number, ]$`P Intensity` <- row$peak
  }
  if (!is.na(row$intensity)) {
    Cell.DataSPI3_S[number, ]$`A Intensity` <- row$mean
  }
  if (!is.na(row$intensity)) {
    Cell.DataSPI3_S[number, ]$Magnitude <- row$mag
  }
}

# fix date formats
Cell.DataSPI3_S$Start <-
  as.yearmon(sub(" .*", "", Cell.DataSPI3_S$Start), "%Y-%m-%d")
Cell.DataSPI3_S$End <-
  as.yearmon(sub(" .*", "", Cell.DataSPI3_S$End), "%Y-%m-%d")

Cell.DataSPI3_S

write.csv(
  Cell.DataSPI3_S,
  paste0(
    OUTPUTS_FOLDER,
    UNIT_N[u],
    "/",
    UNIT_N[u],
    "Drought History SPI_3.csv"
  ),
  row.names = F
)

##########   Count Droughts For Figure - Derek's edits
Cell.DataSPI3_S$`P Intensity` <-
  as.numeric(Cell.DataSPI3_S$`P Intensity`)

EX_Cnt <- sum(Cell.DataSPI3_S[5] > 2)
SV_Cnt <- sum(Cell.DataSPI3_S[5] > 1.5)
MO_Cnt <- sum(Cell.DataSPI3_S[5] > 1)

D_Cnt <- MO_Cnt
SV_Cnt2 <- SV_Cnt - EX_Cnt
MO_Cnt2 <- D_Cnt -  SV_Cnt2 - EX_Cnt

Cell.SPICNT[1:4, 3] <- c(D_Cnt, MO_Cnt2, SV_Cnt2, EX_Cnt)
Cell.SPICNT

##########   Remove Postive SPI and make absloute values
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
M_SPI.ts <-
  ts(SPIVEC_Abs,
     c(1990, 1),
     end = c(ey, em),
     frequency = 12)
myts66 <-
  as.vector(window(M_SPI.ts, start = c(1990, 1), end = c(ey, em)))
DateT1 <-
  as.Date(seq(as.Date("1990-01-01"), as.Date(ed), by = "months"))
xx <- data.frame(DateT1, myts66)
colnames(xx) <- c("DT", "SP")
xx$DT <- as.Date(xx$DT)
xx

write.csv(xx,
          paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "SPI_3.csv"),
          row.names = F)

# plot and write figure
dpi = 300
png(
  paste0(
    OUTPUTS_FOLDER,
    UNIT_N[u],
    "/",
    UNIT_N[u],
    "Drought_HistoryS_3.png"
  ),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

print(
  ggplot(xx, aes(x = DT, y = SP)) +
    geom_area(fill = "darkorange", color = "black") +
    xlab("") +
    
    labs(
      title = paste0("SPI-3 Drought Events 1990-", ey, ": ", UNIT_Ns[u]),
      x = "",
      y = "Drought Intensity"
    ) +
    geom_hline(
      yintercept = 2,
      linetype = "dashed",
      color = "darkred",
      size = 1
    ) +
    geom_hline(
      yintercept = 1.5,
      linetype = "dashed",
      color = "red",
      size = 1
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "orange",
      size = 1
    )  +
    annotate(
      "text",
      x = xx$DT[70],
      y = 3.4,
      label = paste0("Drought Events = ",  D_Cnt),
      size = 5,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = xx$DT[70],
      y = 3.1,
      label = paste0("Moderate Droughts = ", MO_Cnt2),
      size = 5,
      fontface = "bold",
      colour = "orange"
    ) +
    annotate(
      "text",
      x = xx$DT[70],
      y = 2.8,
      label = paste0("Severe Droughts = ", SV_Cnt2),
      size = 5,
      fontface = "bold",
      colour = "red"
    ) +
    annotate(
      "text",
      x = xx$DT[70],
      y = 2.5,
      label = paste0("Extreme Droughts = ", EX_Cnt),
      size = 5,
      fontface = "bold",
      colour = "darkred"
    ) +
    annotate(
      "text",
      x = xx$DT[300],
      y = 3.1,
      label = "SPI-3",
      size = 12,
      fontface = "bold",
      colour = "black"
    )
)

dev.off()


######################################################################################################################################




#### SPI-12 short

head(SPI_ALL)
spi2 <- SPI_ALL[which(SPI_ALL$m.scale == 12), ]
spi2 <- spi2[which(spi2$date >= as.Date("1990-01-01")), ]
head(spi2)

# create binary column of drought (SPI>=1) yes or no
spi2$drought <- ifelse(spi2$spi_negs > 0, 1, 0)
head(spi2, 30)

# make column of consecutive drought months
spi2$DRGHT_ct <- with(spi2, (drought == 1) *
                        ave(drought, rleid(drought == 1),
                            FUN = seq_along))

# subset data for drought event months only
spi3 <- spi2[which(spi2$DRGHT_ct >= 1), ]
head(spi3, 20)

# create empty event count column
spi3$event_ct <- 0

### fill in event_ct column
for (r in 1:nrow(spi3)) {
  SPI_M <- spi3[r, ]
  SPI_M
  
  if (r == 1) {
    spi3[r, ]$event_ct <- r
  }
  if (SPI_M$DRGHT_ct > 1) {
    spi3[r, ]$event_ct <- spi3[r - 1, ]$event_ct
  }
  if (SPI_M$DRGHT_ct == 1 &&
      r > 1) {
    spi3[r, ]$event_ct <- spi3[r - 1, ]$event_ct + 1
  }
}

head(spi3, 100)
summary(spi3$event_ct)

### keep only events with peak spi_negs >=1
event <- unique(spi3$event_ct)

for (x in event) {
  sub <- spi3[which(spi3$event_ct == x), ]
  
  if (max(sub$spi_negs) < 1) {
    spi3[which(spi3$event_ct == x), ]$event_ct <- 0
  }
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
spi5 <- merge(spi2, spi4, all = T)

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
  
  if (max(SPI_I$spi_negs) <= 1.5) {
    SPI_I$intensity <- 1
  }
  if (max(SPI_I$spi_negs) > 1.5 &&
      max(SPI_I$spi_negs) <= 2) {
    SPI_I$intensity <- 2
  }
  if (max(SPI_I$spi_negs) > 2) {
    SPI_I$intensity <- 3
  }
  
  SPI_I$peak <- max(SPI_I$spi_negs)
  SPI_I$mean <- mean(SPI_I$spi_negs)
  SPI_I$mag <- sum(SPI_I$spi_negs)
  
  # keep only these columns
  SPI_I <- SPI_I[c("date", "months", "intensity", "peak", "mean", "mag")]
  SPI_I
  
  spi6 <- merge(spi6, SPI_I, by = "date", all.x = T)
  
  # merge intensity columns together
  if (x > 1) {
    spi6$months <-
      ifelse(is.na(spi6$months.x), spi6$months.y, spi6$months.x)
  }
  if (x > 1) {
    spi6$intensity <-
      ifelse(is.na(spi6$intensity.x),
             spi6$intensity.y,
             spi6$intensity.x)
  }
  if (x > 1) {
    spi6$peak <- ifelse(is.na(spi6$peak.x), spi6$peak.y, spi6$peak.x)
  }
  if (x > 1) {
    spi6$mean <- ifelse(is.na(spi6$mean.x), spi6$mean.y, spi6$mean.x)
  }
  if (x > 1) {
    spi6$mag <- ifelse(is.na(spi6$mag.x), spi6$mag.y, spi6$mag.x)
  }
  
  if (x > 1) {
    spi6 <- subset(
      spi6,
      select = -c(
        months.x,
        months.y,
        intensity.x,
        intensity.y,
        peak.x,
        peak.y,
        mean.x,
        mean.y,
        mag.x,
        mag.y
      )
    )
  }
}

head(spi6, 100)
summary(spi6$intensity)
summary(spi6$peak)
summary(spi6$months)

### make table of start and end dates for each drought intensity event

Cell.DataSPI12_S <- data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI12_S) <-
  c("Start",
    "End",
    "Duration",
    "A Intensity",
    "P Intensity",
    "Magnitude")

Cell.DataSPI12_S$Start <- as.POSIXct(Cell.DataSPI12_S$Start)
Cell.DataSPI12_S$End <- as.POSIXct(Cell.DataSPI12_S$End)
Cell.DataSPI12_S
colnames(Cell.DataSPI12_S)

### loop through rows and fill table

for (y in 1:nrow(spi6)) {
  # get row
  row <- spi6[y, ]
  row
  
  # get previous and next row
  prev <- spi6[y - 1, ]
  
  nex <- spi6[y + 1, ]
  
  # get event number
  number <- row$event_ct
  
  if (!is.na(row$intensity) &&
      is.na(prev$intensity)) {
    Cell.DataSPI12_S[number, ]$Start <- as.Date(nex$date)
  }
  # if(!is.na(row$intensity) && is.na(prev$intensity)) {event.t[number,]$intensity<-row$intensity}
  if (!is.na(row$intensity) &&
      is.na(nex$intensity)) {
    Cell.DataSPI12_S[number, ]$End <- as.Date(nex$date)
  }
  if (!is.na(row$intensity)) {
    Cell.DataSPI12_S[number, ]$Duration <- row$months
  }
  # if(!is.na(row$intensity)) {event.t[number,]$event_ct<-row$event_ct}
  if (!is.na(row$intensity)) {
    Cell.DataSPI12_S[number, ]$`P Intensity` <- row$peak
  }
  if (!is.na(row$intensity)) {
    Cell.DataSPI12_S[number, ]$`A Intensity` <- row$mean
  }
  if (!is.na(row$intensity)) {
    Cell.DataSPI12_S[number, ]$Magnitude <- row$mag
  }
}

# fix date formats
Cell.DataSPI12_S$Start <-
  as.yearmon(sub(" .*", "", Cell.DataSPI12_S$Start), "%Y-%m-%d")
Cell.DataSPI12_S$End <-
  as.yearmon(sub(" .*", "", Cell.DataSPI12_S$End), "%Y-%m-%d")

Cell.DataSPI12_S
Cell.DataSPI12_S[which(Cell.DataSPI12_S$`P Intensity` > 1.5), ]


##########   Count Droughts For Figure - Derek's edits
Cell.DataSPI12_S$`P Intensity` <-
  as.numeric(Cell.DataSPI12_S$`P Intensity`)

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
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " Drought Count.csv"),
  row.names = F
)

##########   Remove Postive SPI and make absloute values
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
M_SPI.ts <-
  ts(SPIVEC_Abs,
     c(1990, 1),
     end = c(ey, em),
     frequency = 12)
myts66 <-
  as.vector(window(M_SPI.ts, start = c(1990, 1), end = c(ey, em)))
DateT1 <-
  as.Date(seq(as.Date("1990-01-01"), as.Date(ed), by = "months"))
xx <- data.frame(DateT1, myts66)
colnames(xx) <- c("DT", "SP")
xx$DT <- as.Date(xx$DT)
xx

# plot and write figure
dpi = 300
png(
  paste0(
    OUTPUTS_FOLDER,
    UNIT_N[u],
    "/",
    UNIT_N[u],
    "Drought_HistoryS_12.png"
  ),
  width = 6.5 * dpi,
  height = 4 * dpi,
  res = dpi
)

print(
  ggplot(xx, aes(x = DT, y = SP)) +
    geom_area(fill = "darkorange", color = "black") +
    xlab("") +
    
    labs(
      title = paste0("SPI-12 Drought Events 1990-", ey, ": ", UNIT_Ns[u]),
      x = "",
      y = "Drought Intensity"
    ) +
    geom_hline(
      yintercept = 2,
      linetype = "dashed",
      color = "darkred",
      size = 1
    ) +
    geom_hline(
      yintercept = 1.5,
      linetype = "dashed",
      color = "red",
      size = 1
    ) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      color = "orange",
      size = 1
    )  +
    annotate(
      "text",
      x = xx$DT[70],
      y = 3.4,
      label = paste0("Drought Events = ",  D_Cnt),
      size = 5,
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = xx$DT[70],
      y = 3.1,
      label = paste0("Moderate Droughts = ", MO_Cnt2),
      size = 5,
      fontface = "bold",
      colour = "orange"
    ) +
    annotate(
      "text",
      x = xx$DT[70],
      y = 2.8,
      label = paste0("Severe Droughts = ", SV_Cnt2),
      size = 5,
      fontface = "bold",
      colour = "red"
    ) +
    annotate(
      "text",
      x = xx$DT[70],
      y = 2.5,
      label = paste0("Extreme Droughts = ", EX_Cnt),
      size = 5,
      fontface = "bold",
      colour = "darkred"
    ) +
    annotate(
      "text",
      x = xx$DT[300],
      y = 3.1,
      label = "SPI-12",
      size = 12,
      fontface = "bold",
      colour = "black"
    )
)

dev.off()


################################################################################################################################################

##########  MEI

##########   UNIT

print("MEI")

Cell.MEI <- data.frame(matrix(nrow = 5, ncol = 9))
colnames(Cell.MEI) <-
  c("Phase",
    "W-Mean",
    "D-Mean",
    "W-Max",
    "D-Max",
    "W-Min",
    "D-Min",
    "W-Count",
    "D-Count")
Cell.MEI[1:5, 1] <-
  c("Strong EL", "Weak EL", "Neutral", "Weak LA", "Strong LA")

##########   All MEI
head(MEI)
tail(MEI)

# get last year
ly <- max(MRF100$Year)
MEI$Year <- seq(1950, ly, 1)

#Cell.MEI_All[u:1] <- UNIT_N[u]
#Cell.MEI_All[u,8] <- W_MRF_MEAN
#Cell.MEI_All[UNIT_C+u,8] <- D_MRF_MEAN
# remove first row from wet season (1950 wet season not available because it starts in 1949)
MEI_W <- subset(MEI, select = -c(MEI_D))
MEI_W <- MEI_W[2:nrow(MEI_W), ]
MEI_D <- subset(MEI, select = -c(MEI_W))

# ## Can read in the csv from above to start script here
# MRF100<-read.csv(paste0(OUTPUTS_FOLDER,UNIT_N[u],"/",UNIT_N[u]," Monthly Rainfall_in.csv"))

head(MRF100)
tail(MRF100)

# # If starting this section from monthly rainfall csv and only want up to 2012
# Cell.AF_Maps<-MRF100[which(MRF100$Year<2013),]

# # Using MEI and Abby's rainfall dataset
# MRF  =  Cell.AF_Maps[c(1:1116),]

# Using ONI and the combo dataset 1950 (row 361) - current
# MRF = MRF100[c(361:(nrow(MRF100)-12)),]
MRF = MRF100[361:nrow(MRF100), ]
head(MRF)
tail(MRF)

##########  Need to remove Rows to get seasons correct
# first row always starts on May 1950
MRF2 =  MRF[-c(1:4), ]
head(MRF2)

# last row changes depending on dataset end point
# get last row in dataset
lr <- MRF2[nrow(MRF2), ]

# get count of how many months in the last year
yc <- max(MRF2$Year)
yc <- nrow(MRF2[which(MRF2$Year == yc), ])

# get row number of April or October of last year
mrn <- data.frame(which(MRF2$Year == lr$Year))

if (nrow(mrn > 4)) {
  arn <- mrn[4, ]
}
if (nrow(mrn > 10)) {
  orn <- mrn[10, ]
}

# make dataset end at the end of the last complete season
if (yc < 4) {
  MRF2 <- MRF2[which(MRF2$Year < lr), ]
}
if (yc > 4 && yc < 10) {
  MRF2 <- MRF2[1:arn, ]
}
if (yc > 10) {
  MRF2 <- MRF2[1:orn, ]
}

tail(MRF2)

MRF3 <- MRF2
MRF3$Month <- as.numeric(MRF3$Month)

### for each consecutive 6 months, aggregate ANOM by season and keep maximum
# make list of values 6 values apart
rows <- seq(from = 1,
            to = nrow(MRF3),
            by = 6)
rows

# create empty dataframe
seasons <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                    c("Year", "RF", "season"))
head(seasons)
head(MRF3)

# loop through consecutive months and aggregate season oni2 values

n <- 1

for (y in rows) {
  b <- MRF3[c(y:(y + 5)), ]
  b
  # calculate average seasonal rainfall value and set season
  seasons[n, ]$RF <- mean(b$RF)
  seasons[n, ]$Year <- max(b$Year)
  if (as.numeric(min(b$Month) == 5)) {
    seasons[n, ]$season <- "dry"
  }
  if (as.numeric(min(b$Month) == 1)) {
    seasons[n, ]$season <- "wet"
  }
  
  n <- n + 1
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
DryRF <- seasons[which(seasons$season == "dry"), ]

DryRF
L0_D <- cbind(DryRF, MEI_D[which(!is.na(MEI_D$MEI_D)), ])
L0_D

WetRF <- seasons[which(seasons$season == "wet"), ]
L0_W <- cbind(WetRF, MEI_W[which(!is.na(MEI_W$MEI_W)), ])
L0_W

##########   Wet Season
##########   Separate by ENSO Phase

ELN_W  <- subset(L0_W, MEI_W > 0.5)
LAN_W  <- subset(L0_W, MEI_W < -0.5)
NUT_Wx <- subset(L0_W, MEI_W > -0.5)
NUT_W  <- subset(NUT_Wx, MEI_W < 0.5)

##########   Strong and Weak EL Wet Season

ELN_W_Weak <- subset(ELN_W , MEI_W  <= 1.5)
ELN_W_Strong <- subset(ELN_W , MEI_W  > 1.5)

##########   Strong and Weak La Wet Season

LAN_W_Weak <- subset(LAN_W , MEI_W  >= -1.5)
LAN_W_Strong <- subset(LAN_W , MEI_W  < -1.5)

##########   DRY Season
ELN_D <- subset(L0_D, MEI_D > 0.5)
LAN_D <- subset(L0_D, MEI_D < -0.5)
NUT_Dx <- subset(L0_D, MEI_D > -0.5)
NUT_D  <- subset(NUT_Dx, MEI_D < 0.5)

##########   Strong and Weak EL Dry Season

ELN_D_Weak <- subset(ELN_D , MEI_D  <= 1.5)
ELN_D_Weak
ELN_D_Strong <- subset(ELN_D , MEI_D  > 1.5)
ELN_D_Strong

##########   Strong and Weak La Dry Season

LAN_D_Weak <- subset(LAN_D , MEI_D  >= -1.5)
LAN_D_Strong <- subset(LAN_D , MEI_D  < -1.5)

##########   Extract RF values
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

Cell.MEI[1:5, 8] <-
  c(C_EL_W_S, C_EL_W_W, C_NU_W, C_LA_W_W, C_LA_W_S)
Cell.MEI[1:5, 9] <-
  c(C_EL_D_S, C_EL_D_W, C_NU_D, C_LA_D_W, C_LA_D_S)
Cell.MEI
##########   Mean RF values for each season-phase
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

Cell.MEI[1:5, 2] <-
  c(Me_EL_W_S, Me_EL_W_W, Me_NU_W, Me_LA_W_W, Me_LA_W_S)
Cell.MEI[1:5, 3] <-
  c(Me_EL_D_S, Me_EL_D_W, Me_NU_D, Me_LA_D_W, Me_LA_D_S)

##########   MAX

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

Cell.MEI[1:5, 4] <-
  c(Mx_EL_W_S, Mx_EL_W_W, Mx_NU_W, Me_LA_W_W, Mx_LA_W_S)
Cell.MEI[1:5, 5] <-
  c(Mx_EL_D_S, Mx_EL_D_W, Mx_NU_D, Me_LA_D_W, Mx_LA_D_S)

##########   MIN
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

Cell.MEI[1:5, 6] <-
  c(Mn_EL_W_S, Mn_EL_W_W, Mn_NU_W, Mn_LA_W_W, Mn_LA_W_S)
Cell.MEI[1:5, 7] <-
  c(Mn_EL_D_S, Mn_EL_D_W, Mn_NU_D, Mn_LA_D_W, Mn_LA_D_S)

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
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "MEI_DRY.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

DRY_ENSO <- c(EL_D_S, EL_D_W, NU_D, LA_D_S, LA_D_W)

boxplot(
  c(EL_D_S),
  (EL_D_W),
  (NU_D),
  (LA_D_W),
  (LA_D_S),
  col = c("darkred", "red", "grey", "blue", "darkblue"),
  names = c("SEL" , "WEL", "NUT", "WLA", "SLA"),
  ylab = paste0("Avg. Monthly Rainfall (", RFUnit2, ")"),
  ylim = c(0, MAXD)
)


title(
  paste0("Dry Sea. RF: ", UNIT_Ns[u]),
  line = 1,
  cex.main = 1.1,
  col.main = "gold4"
)

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

MAXA <- max(c(Mx_EL_W_S, Mx_EL_W_W, Mx_LA_W_S, Mx_LA_W_W, Mx_NU_W))
MAXB <- MAXA * 0.05
MAXB2 <- MAXA * 0.10
MAXC <- MAXA * 0.15
MAXD <- (MAXC * 1.2) + MAXA
MAXD2 <- MAXD - MAXB2
MAXD22 <- MAXD - MAXB
MAXD3 <- MAXD - MAXC

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "MEI_WET.png"),
  width = 5 * dpi,
  height = 5 * dpi,
  res = dpi
)

boxplot(
  c(EL_W_S),
  (EL_W_W),
  (NU_W),
  (LA_W_W),
  (LA_W_S),
  col = c("darkred", "red", "grey", "blue", "darkblue"),
  names = c("SEL" , "WEL", "NUT", "WLA", "SLA"),
  ylab = paste0("Avg. Monthly Rainfall (", RFUnit2, ")"),
  ylim = c(0, MAXD)
)
title(
  paste0("Wet Sea RF: ", UNIT_Ns[u]),
  line = 1,
  cex.main = 1.1,
  col.main = "darkgreen"
)

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
          paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " MEI_A.csv"),
          row.names = F)



#####################################################
##### ENSO rainfall barplots (Derek's addition)

# Average Monthly Rainfall by ENSO Phase and Season barplot
# spell out the ENSO phases
Cell.MEI2 <- Cell.MEI
Cell.MEI2$Phase2 <- NA
Cell.MEI2

p <- as.list(unique(Cell.MEI2$Phase))
p

for (i in 1:nrow(Cell.MEI2)) {
  y <- Cell.MEI2[i, ]
  if (y$Phase == p[[1]]) {
    Cell.MEI2[i, ]$Phase2 <- "Strong El Nino"
  }
  if (y$Phase == p[[2]]) {
    Cell.MEI2[i, ]$Phase2 <- "Weak El Nino"
  }
  if (y$Phase == p[[3]]) {
    Cell.MEI2[i, ]$Phase2 <- "Neutral"
  }
  if (y$Phase == p[[4]]) {
    Cell.MEI2[i, ]$Phase2 <- "Weak La Nina"
  }
  if (y$Phase == p[[5]]) {
    Cell.MEI2[i, ]$Phase2 <- "Strong La Nina"
  }
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
cc$Phase <-
  factor(
    cc$Phase,
    levels = c(
      "Strong El Nino",
      "Weak El Nino",
      "Neutral",
      "Weak La Nina",
      "Strong La Nina"
    )
  )

cc

# set ylim
ylim <- max(cc$Mean) * 1.2

png(
  paste0(
    OUTPUTS_FOLDER,
    UNIT_N[u],
    "/",
    UNIT_N[u],
    "ENSO_season_barplot.png"
  ),
  width = 6 * dpi,
  height = 2.5 * dpi,
  res = dpi
)

ggplot(data = cc,
       aes(x = Season, y = Mean, group = Phase)) +
  geom_bar(
    aes(fill = Phase),
    position = position_dodge(width = 0.7),
    stat = "identity",
    color = "black",
    alpha = .7,
    width = .55
  ) +
  labs(title = "Average Monthly Rainfall by Season and ENSO Phase",
       y = "Rainfall (inches)", x = "Season") +
  scale_fill_manual(
    values = c("darkred", "red", "grey", "lightskyblue1", "darkblue"),
    limits = c(
      "Strong El Nino",
      "Weak El Nino",
      "Neutral",
      "Weak La Nina",
      "Strong La Nina"
    )
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  guides(fill = guide_legend(title = "ENSO Phase")) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.7),
            vjust = -0.8) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13)
  )

dev.off()

cc
### Export table of average monthly rainfall by season and ENSO phase
write.csv(cc, paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], " MEI_S.csv"))

#####################################################
##### Air temperature graph
########## Month Year Air Temperature Maps
#AT_Map_Path_A <- ("F:/PDKE/CCVD/CCVD INPUTS/air_temp/data_map_newer/")
AT_Map_Path_A <- paste0(INPUTS_FOLDER, "/air_temp/data_map_newer/")
maps <- AT_Map_Path_A

# make list of years
years <- as.list(list.files(maps))

# create empty table for raster stats
table <- data.frame(matrix(
  ncol = 5,
  nrow = 0,
  dimnames = list(NULL, c("year", "month", "min", "max", "mean"))
))

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

y <- as.list(unique(dat$year))

for (i in y) {
  if (length(which(dat$year == i)) == 12) {
    dat[which(dat$year == i), ]$full <- "Y"
  }
}

ic <- unique(dat[which(dat$full == "N"), ]$year)

# dat<-dat[which(dat$full == "Y"),]
head(dat)
tail(dat)
#############################
### Annual air temp trend

# aggregate monthly to annual air temp
rm(min)
rm(max)
rm(mean)

dat.y <- aggregate(mean ~ year, dat, mean)
dat.y$mean <- round(dat.y$mean, digit = 2)
dat.y$meanmin <-
  round((aggregate(mean ~ year, dat, min))$mean, digits = 2)
dat.y$meanmax <-
  round((aggregate(mean ~ year, dat, max))$mean, digits = 2)
dat.y$min <- round((aggregate(min ~ year, dat, min))$min, digits = 2)
dat.y$max <- round((aggregate(max ~ year, dat, max))$max, digits = 2)

dat.y

# remove stats from incomplete year
dat.y <- dat.y[which(dat.y$year < ic), ]
dat.y[(nrow(dat.y) + 1), ]$year <- ic

# make date column
dat.y$date <- as.Date(paste0(dat.y$year, "-01-01"))

# write to csv table of annual air temp values from monthly means (min, max, mean)
write.csv(dat.y,
          paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "_monthly_airtemp.csv"))

slope <-
  formatC((coef(lm(
    dat.y$mean ~ dat.y$date
  ))[2]), format = "e", digits = 2)
slope

dpi = 300

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "_annual_airtemp.png"),
  width = 6 * dpi,
  height = 3 * dpi,
  res = dpi
)

ggplot(dat.y, aes(x = date, y = mean)) +
  geom_line(size = 1.2, color = "orange") +
  geom_smooth(
    method = lm,
    se = F,
    size = 1,
    color = "black"
  ) +
  stat_cor(
    method = "pearson",
    label.x = as.Date("1991-01-01"),
    label.y = 58
  ) +
  scale_x_date(date_breaks = "5 years", labels = date_format(format = "%Y")) +
  labs(title = "Annual Air Temperature and Extremes",
       y = "Temperature (F)", x = "Year") +
  geom_line(aes(x = date, y = min), size = 1.2, color = "blue") +
  geom_line(aes(x = date, y = max), size = 1.2, color = "red") +
  annotate(
    "text",
    x = as.Date("1994-07-01"),
    y = 54,
    label = paste0("Slope = ", slope)
  ) +
  theme(
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

dev.off()

##### add annual mean values
dat2 <-
  merge(dat, dat.y[, c("year", "mean", "min", "max")], by = "year", all.x =
          T)
head(dat2)
tail(dat2)

# write to csv table of monthly air temp values (min, max, mean)
write.csv(dat2,
          paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "_daily_airtemp.csv"))

# set y-axis limits
ylow <- min(dat.y$min, na.rm = T) * .5
yhi <- max(dat.y$max, na.rm = T) * 1.0005

# set slope
slope <-
  formatC((coef(lm(
    dat2$mean.x ~ dat2$date
  ))[2]), format = "e", digits = 2)
slope

# set location for linear trend values
yl <- min(dat.y$min, na.rm = T) * .75

dpi = 300

png(
  paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "_monthly_airtemp.png"),
  width = 6 * dpi,
  height = 3 * dpi,
  res = dpi
)

ggplot(dat2, aes(x = date, y = mean.x)) +
  geom_line(color = "grey60") +
  geom_smooth(
    span = 0.2,
    se = F,
    size = 1.3,
    color = "orange"
  ) +
  geom_smooth(method = lm,
              se = F,
              color = "black") +
  stat_cor(
    method = "pearson",
    label.x = as.Date("1993-01-01"),
    label.y = yl - 5
  ) +
  scale_x_date(date_breaks = "4 years", labels = date_format(format = "%Y")) +
  ylim(ylow, yhi) +
  labs(title = "Air Temperature Trends",
       y = "Temperature (F)", x = "") +
  geom_hline(yintercept = 0) +
  annotate(
    "text",
    x = as.Date("1996-07-01"),
    y = yl,
    label = paste0("Slope = ", slope)
  ) +
  geom_line(
    data = dat.y,
    aes(x = date, y = min),
    size = 1.2,
    color = "blue"
  ) +
  geom_line(
    data = dat.y,
    aes(x = date, y = max),
    size = 1.2,
    color = "red"
  ) +
  theme(
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

dev.off()

######## Air Temp Anomalies ########
#setwd("F:/PDKE/CCVD/MINI_Phase2/")
setwd(WORKING_DIR)

# make list of climatology files
list <- list.files(AT_CLIM_Path_A, pattern = "\\.tif$")

### Make dataframe of average monthly climatology values within study area
# create empty table for mean AT values
clim <- data.frame(matrix(
  ncol = 2,
  nrow = 0,
  dimnames = list(NULL, c("month", "mean"))
))
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
anom$date <-
  as.Date(paste0(anom$year, "/", anom$month, "/", "01"), format = "%Y/%m/%d")

# write to csv table of monthly air temp anomaly values
write.csv(anom,
          paste0(OUTPUTS_FOLDER, UNIT_N[u], "/", UNIT_N[u], "_anomaly_airtemp.csv"))

# set y-axis limits
ylow <- min(anom$anom_f) * .95
yhi <- max(anom$anom_f) * 1.0005

# set slope
slope <-
  formatC((coef(lm(
    anom$mean_f ~ anom$date
  ))[2]), format = "e", digits = 2)
slope

# set location for linear trend values
yl <- (((min(anom$anom_f) - max(anom$anom_f)) / 2) + max(anom$anom_f))

dpi = 300

png(
  paste0(
    OUTPUTS_FOLDER,
    UNIT_N[u],
    "/",
    UNIT_N[u],
    "_monthly_airtemp_anomaly.png"
  ),
  width = 6 * dpi,
  height = 3 * dpi,
  res = dpi
)

ggplot(anom, aes(x = date, y = anom_f)) +
  geom_line(color = "grey60") +
  geom_smooth(
    span = 0.2,
    se = F,
    size = 1.3,
    color = "orange"
  ) +
  geom_smooth(method = lm,
              se = F,
              color = "black") +
  # stat_cor(method="pearson", label.x=as.Date("1993-01-01"), label.y=yl-3) +
  scale_x_date(date_breaks = "4 years", labels = date_format(format = "%Y")) +
  ylim(ylow, yhi) +
  labs(title = "Air Temperature Anomaly Trends",
       y = "Temperature (F)", x = "") +
  # geom_hline(yintercept=0) +
  # annotate("text", x=as.Date("1996-07-01"), y=yl,
  #          label=paste0("Slope = ",slope)) +
  # geom_line(data = dat.y, aes(x=date, y=min), size=1.2, color="blue") +
  # geom_line(data = dat.y, aes(x=date,y=max), size=1.2, color="red") +
  theme(
    panel.background = element_rect(fill = NA, color = "black"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

dev.off()

### END! ###
