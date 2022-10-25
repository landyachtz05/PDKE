library(gstat)
library(raster)
library(sp) 
library(maptools)
library(rgdal)
library(RColorBrewer)
library(gridExtra)
library(stats.xbar)
library(ggplot2)
library(grid)
library(gridExtra)
library(devtools)
library(lubridate)
library(rgeos)
library(latticeExtra)
library(rasterVis)
library(plotrix)
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
library(plyr)
library(SPEI)
library(sf)
library(ggpubr)


##########   This code will the generate the inputs for a CCVD portfolio
##########   Reuired Inputs: 1. Folder destinations 2. Measurementuints 3. Shpaefiles 4. List of files and output names
##########   Searchable sections of the Code:
##########   Maps - Elevation - Mean Climate - Downscaling - Rainfall Extract - SPI 1 - SPI 2 - MEI

##########################################################################################################################

setwd("E:/PDKE/CCVD/MINI_Phase2/")               # WORKING DIRECTORY
IFOLDER <- "E:/PDKE/CCVD/CCVD INPUTS/"       # INPUT FOLDER
RFOLDER <- "E:/PDKE/CCVD/MINI_Phase2/"           # OUTPUT FOLDER 1

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

Mean_CLIM = dir(paste0(IFOLDER,"Mean Climate/"), pattern="*x.adf", recursive=T, full.names=T)
EXAMP <- raster(Mean_CLIM[71])

plot(EXAMP)

############################################################################################################################

##########   Coastal Shape Files 

Coast <- readOGR(paste0(IFOLDER,"COAST/coast_2/coast_geo_shp.dbf"))

# ALL ISLANDS
Coast_Crop <- crop(Coast , extent(-160.0, -154.8066, 18.91069, 22.23))
Coast_Crop_T <- spTransform(Coast_Crop, crs(EXAMP))  
# BIG ISLAND
Coast_Crop <- crop(Coast , extent(-156.1, -154.8066, 18.91069, 20.3))
Coast_BI <- spTransform(Coast_Crop, crs(EXAMP))  
# MAUI
Coast_Crop <- crop(Coast, extent(-156.75, -155.9257, 20.343, 21.1)) #Maui 
Coast_MN <- spTransform(Coast_Crop, crs(EXAMP))  
#Molokai
Coast_Crop <- crop(Coast, extent(-157.58, -156.7, 21.031, 21.997)) #Molokai 
Coast_MO <- spTransform(Coast_Crop, crs(EXAMP))   
#Lanai
Coast_Crop <- crop(Coast, extent(-157.07, -156.8, 20.05, 21.02)) #Lanai
Coast_LA <- spTransform(Coast_Crop, crs(EXAMP))
#Oahu
Coast_Crop <- crop(Coast, extent(-158.29, -157.35, 21.05, 22.02)) #Oahu
Coast_OA <- spTransform(Coast_Crop, crs(EXAMP))
#KAUI
Coast_Crop <- crop(Coast, extent(-160.0, -158.35, 21.8, 22.23)) #Kaui
Coast_KA <- spTransform(Coast_Crop, crs(EXAMP))
extent(Coast)
#Kahoolawe 
Coast_Crop <- crop(Coast, extent(-156.75, -156.5, 20.343, 20.72)) #Maui 
Coast_KO <- spTransform(Coast_Crop, crs(EXAMP))  


# #National Park 
# NP_ALL <- readOGR("E:/PDKE/CCVD/NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service.shp")
# NP_ALL <- spTransform(NP_ALL, crs(EXAMP))  # Set spatial extent to Example 
# NPALL <- NP_ALL[NP_ALL@data$STATE == "HI",]
# NPALL@data
# NPALL@data$UNIT_CODE
# HALE <- NPALL[NPALL@data$UNIT_CODE == "HALE",]  #"Haleakala National Park"

# #Waimea Valley
# NP_ALL <- readOGR("E:/PDKE/CCVD/waimeavalley_watershed.shp")
# HALE <- NP_ALL
# HALE@data

#Kahikinui HHL, Maui
NP_ALL <- readOGR("E:/PDKE/CCVD/HHL_kahikinui_maui_WGS.shp")
HALE <- NP_ALL
HALE@data

# Check map - NEED TO CHANGE MANUALLY BASED ON ISLAND
plot(Coast_MN, main = "Maui")
plot(HALE ,add = T, col="red")


##################################################

##########   Define Units   


##########   SHAPE FILES FOR ANALYSIS
UNIT_X <-   c(HALE,HALE)
 
##########   ISLAND - NEED TO CHANGE MANUALLY BASED ON ISLAND
UNIT_I <-  c("MN","MN")

##########  UNIT NAME
##########  UNIT NAME (For CCVD Narrative)
UNIT_N <- as.vector(as.character(c("Hawaiian Homelands - Kahikinui","Hawaiian Homelands - Kahikinui")))


##########  Short Name (For Figure Titles)
UNIT_Ns <- as.vector(as.character(c("Kahikinui","Kahikinui")))
#20,24,-26 34


##########   COUNT INPUTS

UNIT_Shape <- sum(!is.na(UNIT_X))
UNIT_C <- sum(!is.na(UNIT_N))  # This will mark the end of the Loop
UNIT_Island <- sum(!is.na(UNIT_I))
ShortNm <- sum(!is.na(UNIT_Ns))

##########   CONFIRM INPUTS ARE THE SAME

UNIT_C
ShortNm
UNIT_Shape
UNIT_Island

##########  Mean Rainfall DATA   

MeanRF_ALL = dir(paste0(IFOLDER,"Mean_RF_Data/StateMaps/"), pattern="*x.adf", recursive=T, full.names=T)


##########   FIRE OCCURRENCE

#Fire Occurrence Shape 2019 (From Clay)
u<-1

FIRE_Shape <- readOGR(paste0(IFOLDER,"StateFire_1999/2020_1999_Hawaii_Fire_Perimeters.shp"))
FIRE_Shape_T <- spTransform(FIRE_Shape, crs(EXAMP))  
crs(FIRE_Shape_T) <- CRS(proj4string(EXAMP))
Fire_Mask <- mask(x = FIRE_Shape_T, mask = UNIT_X[[u]])
Fire_Crop <- crop(x = FIRE_Shape_T, y = extent(UNIT_X[[u]]))

##########   Forest Roads

F_Roads <- readOGR(paste0(IFOLDER,"Forestry Roads/Forestry_Roads.shp")) 
F_Roads <- spTransform(F_Roads, crs(EXAMP))  
plot(F_Roads)

##########   Trail Inventory

I_Trails <- readOGR(paste0(IFOLDER,"Inventory Trails/InventoryTrails.shp")) 
I_Trails <- spTransform(I_Trails, crs(EXAMP))   
plot(I_Trails)

##########   Digital Elevation Model 

ELEV2 = dir(paste0(IFOLDER,"ned_dem/"), pattern="*x.adf", recursive=T, full.names=T)
ELEV <- raster(ELEV2[1])
if(ELUnit == " ft") {ELEV =  ELEV * 3.28084 }
crs(ELEV) <- "+proj=longlat +datum=WGS84 +no_defs" 
plot(ELEV)

##########  Landcover

LC2 = dir(paste0(IFOLDER, "Landcover/"), recursive=T, full.names=T)
LC <- raster(LC2[3])
crs(LC) <- "+proj=longlat +datum=WGS84 +no_defs"
plot(LC)

##########  Downscaling

DyDs_files = dir(paste0(IFOLDER,"HI_Downscaling_Final/"), recursive=T, full.names=T)  
D2_files = dir(paste0(IFOLDER,"Lulin_Downscaling/"), recursive=T, full.names=T)  

##########   Month Year Rainfall Maps (Frazier et al., 2016)

RF_Map_Path_A <- ("E:/PDKE/CCVD/data/production/rainfall/legacy/month/statewide/data_map/")

##########   Month-Year Rainfall Maps New Lucas et al., In Review

# Need to download maps from HCDP
# Legacy maps 1920-1989
# Matty maps 1990-NRT
# Put in INput folder 
#For now Read this path in seperate as not in input folder 
RF_Map_Path <- ("E:/PDKE/CCVD/NEW RF MAPS/statewide_rf_mm/rf_mm/")

##########  Multi-Variate ENSO Index

#Replace this with the ONI Dataset at some point 
MEI <- read.csv(paste0(IFOLDER,"MEI_Season.csv"),sep=",")

########## Month Year Air Temperature Maps
AT_Map_Path_A <- ("E:/PDKE/CCVD/CCVD INPUTS/air_temp/data_map/")

#################################################################################################################################

##########   DataFrames For Output 

#Set up Group matrix for all climate Variables 
#Climate For all Ranches
Cell.DataCL <-data.frame(matrix(ncol=15, nrow=UNIT_C))
colnames(Cell.DataCL) <- c("Unit","Elev Min","Elev Mean", "Elev Max","RF","Tavg","Tmax","Tmin","RH","SM","KD","ET","CF","Island","Short Name")
Cell.DataCL[1] <- UNIT_N
Cell.DataCL[15] <- UNIT_Ns

#Cliamte For Individual Ranches
Cell.DataCLR <-data.frame(matrix(ncol = 15, nrow = 3))
colnames(Cell.DataCLR) <- c("Unit","Elev","RF","Tavg","Tmax","Tmin","RH","SM","KD","ET","CF","DryRF","WetRF","Island","Short Name")

#Create a matrix For a site Comparision of rainfall 
Cell.RF_Year <-data.frame(matrix(ncol = 15, nrow = UNIT_C))
colnames(Cell.RF_Year) <- c("Unit","Mean ELEV","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC","ANN")
Cell.RF_Year[1] <- UNIT_N




u<-1
# for (u in 1:1) {

print("Analysis UNIT ")  # Unit Analyzed
print(u)                 # Number In Loop
print(UNIT_N[u])         # Name Of Unit 
print("Maps") 

#Create A Directory for Output 
path <- paste0(RFOLDER,UNIT_N[u])
dir.create(path, showWarnings = TRUE, recursive = FALSE)

IS <- UNIT_I[u]       # Island 
SH <- UNIT_X[[u]]     # Shape File

##########   Get Coast for Island 

if(IS == "BI") {CoastM <- Coast_BI; Iname <- "Hawaii"}
if(IS == "MN") {CoastM <- Coast_MN; Iname <- "Maui"}
if(IS == "MO") {CoastM <- Coast_MO; Iname <- "Molokai"}
if(IS == "LA") {CoastM <- Coast_LA; Iname <- "Lanai"}
if(IS == "OA") {CoastM <- Coast_OA; Iname <- "Oahu"}
if(IS == "KA") {CoastM <- Coast_KA; Iname <- "Kauai"}
if(IS == "KO") {CoastM <- Coast_KO; Iname <- "Kahoolawe"}
if(IS == "ALL") {CoastM <- Coast_KO; Iname <- "Hawaiian Islands"}

##########   GET Central Lattitude and Longitude for shapefile

Xmin <- extent(SH)[1]
Xmax <- extent(SH)[2]
ymin <- extent(SH)[3]
ymax <- extent(SH)[4]

LONG <- round(mean(Xmin,Xmax),3)
LAT <- round(mean(ymin,ymax),3)


##########   Make a Plot of Shapefile with Island Coast and LAT and LON

TIT <- paste0(UNIT_Ns[u]," (",Iname,")")  # Title For Figure 
dpi=300

# find width of map
width<-abs(extent(CoastM)[1])-abs(extent(CoastM)[2])

# convert width units (DD to miles)
library(DescTools)

# set latitude in radians
  lat<-mean(c(extent(CoastM)[3], extent(CoastM)[4]))
  lat.r<-DegToRad(lat)

  # calculate width (miles) of 1 degree at latitude
  width.d<-abs(cos(lat)*69.172)
  
  # calculate width (miles) of map
  width.m<-round(width*width.d,1)
  
  # set scalebar width and convert to kilometers (to match coordinate system)
  sb<-width.m*.25
  
  if(sb<=3) {sb2<-1}
  if(sb>3 && sb<=5) {sb2<-5}
  if(sb>5 && sb<=10) {sb2<-10}
  if(sb>10 && sb>=20) {sb2<-20}
  
  sb2<-sb2*1.60934

# plot
png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Map.png"),width=5*dpi,height=5*dpi,res=dpi) 

par(mfrow=c(1,1))

plot(CoastM,lwd=1,lty = "solid",axes=TRUE,las=2)
title(TIT , line = 0.5,cex.main = 1.5)
plot(SH,add=TRUE,col='red')
legend("topright", legend = c(paste0("Latitude = ", LAT), paste0("Longitude = ", LONG)), 
bty = "n", horiz = FALSE, inset = c(0.05, 0.05),cex=0.8)
scalebar(sb2, type="bar", divs=4, label=c(0,5,10), below="miles")



dev.off()

########## Map with all of the Islands 

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," ISMap.png"),width=5*dpi,height=5*dpi,res=dpi) 
par(mfrow=c(1,1))
plot(Coast_Crop_T,lwd=1,lty = "solid",axes=TRUE,las=2)
title(TIT , line = 0.5,cex.main = 1.5)
plot(SH,add=TRUE,col='red')
#Legend
legend("topright", legend = c(paste0("Latitude = ", LAT), paste0("Longitude = ", LONG)), 
       bty = "n", horiz = FALSE, inset = c(0.05, 0.05))

dev.off()

##########  Fire Map 

FIRE_CropI <- FIRE_Shape_T
TITF <- "Fire Occurrence 1999-2020"

dpi=300
png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Fire.png"),width=5*dpi,height=5*dpi,res=dpi) 

par(mfrow=c(1,1))
plot(CoastM,lwd=1,lty = "solid",axes=TRUE,las=2)
title(TITF , line = 0.5,cex.main = 1.5)
plot(SH,add=TRUE,lwd=2,col="cyan")
plot(FIRE_CropI,add=TRUE,col="darkred")

legend("topright", legend = c("Fire Occurrence",UNIT_Ns[[u]]), pch=c(15,15),col=c("darkred","cyan"),
       bty = "n", horiz = FALSE, inset = c(0.05, 0.05),cex=0.8)

dev.off()

########## Road Map 

Roads_CropI <- crop(x = F_Roads, y = extent(CoastM))

TITF <- paste0("Forest Roads (",Iname,")")

dpi=300
png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Forest_Roads.png"),width=5*dpi,height=5*dpi,res=dpi) 

par(mfrow=c(1,1))
plot(CoastM,lwd=1,lty = "solid",axes=TRUE,las=2)
title(TITF , line = 0.5,cex.main = 1.5)
plot(SH,add=TRUE,lwd=2,col="cyan")
plot(Roads_CropI,add=TRUE,col="darkgreen")

legend("topright", legend = c("Forest Roads",UNIT_Ns[[u]]), pch=c(15,15),col=c("darkgreen","cyan"),
       bty = "n", horiz = FALSE, inset = c(0.05, 0.05),cex=0.8)

dev.off()

########## Trail Map

Trails_CropI <- crop(x = I_Trails, y = extent(CoastM))

TITF <- paste0("Trail Inventory (",Iname,")")

dpi=300
png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Trails.png"),width=5*dpi,height=5*dpi,res=dpi) 

par(mfrow=c(1,1))
plot(CoastM,lwd=1,lty = "solid",axes=TRUE,las=2)
title(TITF , line = 0.5,cex.main = 1.5)
plot(SH,add=TRUE,lwd=2,col="cyan")
plot(Trails_CropI ,add=TRUE,col="brown")

legend("topright", legend = c("Trails",UNIT_Ns[[u]]), pch=c(15,15),col=c("brown","cyan"),
       bty = "n", horiz = FALSE, inset = c(0.05, 0.05),cex=0.8)

dev.off()

########### Elevation Extract and Map

print("Elevation")  

ELEV_MaskI <- mask(x = ELEV, mask = CoastM)
ELEV_CropI <- crop(x = ELEV_MaskI, y = extent(CoastM))

#Mask For Geography
ELEV_Mask <- mask(x = ELEV, mask = UNIT_X[[u]])
ELEV_Crop <- crop(x = ELEV_Mask, y = extent(UNIT_X[[u]]))
ELEV_Mean <- round(cellStats(ELEV_Crop, 'mean'),1)
ELEV_Max <- round(cellStats(ELEV_Crop, 'max'),1)
ELEV_Min <- round(cellStats(ELEV_Crop, 'min'),1)

Cell.DataCL[u,2:4] <- c(ELEV_Min,ELEV_Mean,ELEV_Max)
Cell.DataCL[u,17] <- Iname
Cell.DataCLR[1:3,1] <- c("Mean","Max","Min")
Cell.DataCLR[1:3,2] <- c(ELEV_Mean,ELEV_Max,ELEV_Min)
Cell.DataCLR[1,14] <- Iname
Cell.DataCLR[1,15] <- UNIT_Ns[u]

Cell.DataCL
Cell.DataCLR

#BI_brksRH<-round(seq(RHLO, RHUP, length = 9),0)
colfuncEL <-colorRampPalette(brewer.pal(9,"PuBuGn"))(100)

SHAPE <- UNIT_X[[u]]

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," ELMap.png"),width=5*dpi,height=5*dpi,res=dpi) 

print(spplot(ELEV_CropI, col.regions = colfuncEL, equal=FALSE,
             axes = TRUE,las = 1, cex.axis=0.7,
             main=list(label=paste0("Elevation"," ",Iname," (",ELUnit2,")"),cex=0.8),
             colorkey = list(space = "right", height = 1, labels=list(cex=0.6)))+
        layer(sp.polygons(SHAPE,lwd=2)) +
        layer(sp.polygons(CoastM,lwd=1)))

dev.off()   

########## Landcover Extraction and Map

print("Landcover")  

LC_MaskI <- mask(x = LC, mask = CoastM)
LC_CropI <- crop(x = LC_MaskI, y = extent(CoastM))

#Mask For Geography
LC_Mask <- mask(x = LC, mask = UNIT_X[[u]])
LC_Crop <- crop(x = LC_Mask, y = extent(UNIT_X[[u]]))

plot(LC_Crop)

#Convert raster to dataframe for stats
LC_Crop2 <- as.data.frame(LC_Crop, xy=T)
head(LC_Crop2)

# remove LC "0" (no data)
LC_Crop3<-LC_Crop2[which(LC_Crop2$LCMAP_HI_2020_V10_LCPRI_crs != 0),]

#Count frequencies
LC_ct<-as.data.frame(table(LC_Crop3$LCMAP_HI_2020_V10_LCPRI_crs))
colnames(LC_ct)<-c("LC","count")
LC_ct2<-LC_ct[order(-LC_ct$count),]

#convert into area (1 pixel = 900 sqare meters = 0.222395 acres)
LC_ct2$acres<-round(LC_ct2$count*0.222395, 1)

#Calculate percentages
all<-sum(LC_ct2$count)

LC_ct2$pct<-round((LC_ct2$count/all)*100, 0)

for (x in 1:nrow(LC_ct2)) {
  if(LC_ct2[x,]$pct == 0) {LC_ct2[x,]$pct<-"<1"}}

LC_ct2$pct<-paste0(LC_ct2$pct,"%")
LC_ct2

#define class names and colors
LC_ct2$class_name<-as.character(NA)
LC_ct2$color<-as.character(NA)

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

# set class names as leveled factors
LC_ct2$class_name<-factor(LC_ct2$class_name, 
                          levels=c("Developed","Cropland","Grass/Shrub","Tree Cover",
                                   "Water","Wetland","Barren"))

# set first thru third landcover classes
LC1<-paste0(LC_ct2[1,]$class_name,"(",LC_ct2[1,]$pct,")")
LC2<-paste0(LC_ct2[2,]$class_name,"(",LC_ct2[2,]$pct,")")
LC3<-paste0(LC_ct2[3,]$class_name,"(",LC_ct2[3,]$pct,")")

### Plot Bar Graph
# set bargraph ylim
ylim<-max(LC_ct2$acres) + (max(LC_ct2$acres)*0.1)

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," LC_barchart.png"),width=5*dpi,height=5*dpi,res=dpi) 

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
geom_text(aes(label = pct), vjust = -0.8) +
theme(axis.text.x = element_text(angle=35, vjust = .65, size=12)) +
labs(title=paste0("Landcover"," ",UNIT_N[u]), y="Area (acres)", x="") +
theme(legend.position="none"))

dev.off()   


# 
# Cell.DataCL[u,14:16] <- c(LC1,LC2,LC3)
# Cell.DataCLR[1,12:14] <- c(LC1,LC2,LC3)

SHAPE <- UNIT_X[[u]]

### remove "0" values from raster
LC_CropI2<-reclassify(LC_CropI, cbind(0, NA))
summary(LC_CropI2)

# export map
TITF<-paste0("Landcover"," ",Iname)

dpi=300
png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," LCMap.png"),width=9*dpi,height=5*dpi,res=dpi) 
unique(LC_CropI2)

plot(LC_CropI2, legend=F, col=c("black","darkgoldenrod1","darkolivegreen1","darkgreen",
                     "aquamarine","blue","chocolate4"),
     xaxt='n', yaxt='n')
plot(SHAPE, lwd=1, border="black", lwd=2, add=T)
plot(CoastM, add=T)
title(TITF , line = 0.5,cex.main = 1.5)

legend("topright",legend=c("Developed","Cropland","Grass/Shrub","Tree Cover","Water",
                "Wetland","Barren"),
       fill=c("black","darkgoldenrod1","darkolivegreen1","darkgreen",
               "blue","aquamarine","chocolate4"))


dev.off()  

########## Mean Climate Extract 

print("Mean Climate")  

#Create a climate matrix 
Cell.matrix5 <- matrix(ncol = 14, nrow = 9)
Cell.CL_Year <-data.frame(Cell.matrix5)
colnames(Cell.CL_Year) <- c("Variable","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC","ANN")
Cell.CL_Year[1,1] <-  paste0("RF [",RFUnit2,"]")
Cell.CL_Year[2,1] <- paste0("Min TA [",TUnit,"]")
Cell.CL_Year[3,1] <-  paste0("Mean TA [",TUnit,"]")
Cell.CL_Year[4,1] <-  paste0("Max TA [",TUnit,"]")
Cell.CL_Year[5,1] <- "RH [%]"
Cell.CL_Year[6,1] <- "CF [%]"
Cell.CL_Year[7,1] <- paste0("ET [",RFUnit2,"]")
Cell.CL_Year[8,1] <- "SM [%]"
Cell.CL_Year[9,1] <- "S [W m/2]"

Cell.CL_Year

##########   Solar Radiation

Jan <- raster(Mean_CLIM[58])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropKD <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
JanMKD   <- round(cellStats(Jan_CropKD, 'mean'),0)
JanMKDx   <- round(cellStats(Jan_CropKD, 'max'),0)
JanMKDn   <- round(cellStats(Jan_CropKD, 'min'),0)

Feb <- raster(Mean_CLIM[57])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropKD <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
FebMKD   <- round(cellStats(Feb_CropKD, 'mean'),0)
FebMKDx   <- round(cellStats(Feb_CropKD, 'max'),0)
FebMKDn   <- round(cellStats(Feb_CropKD, 'min'),0)

Mar <- raster(Mean_CLIM[61])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropKD <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
MarMKD   <- round(cellStats(Mar_CropKD, 'mean'),0)
MarMKDx   <- round(cellStats(Mar_CropKD, 'max'),0)
MarMKDn   <- round(cellStats(Mar_CropKD, 'min'),0)

Apr <- raster(Mean_CLIM[54])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropKD <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
AprMKD   <- round(cellStats(Apr_CropKD, 'mean'),0)
AprMKDx   <- round(cellStats(Apr_CropKD, 'max'),0)
AprMKDn   <- round(cellStats(Apr_CropKD, 'min'),0)

May <- raster(Mean_CLIM[62])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropKD <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
MayMKD   <- round(cellStats(May_CropKD, 'mean'),0)
MayMKDx   <- round(cellStats(May_CropKD, 'max'),0)
MayMKDn   <- round(cellStats(May_CropKD, 'min'),0)

Jun <- raster(Mean_CLIM[60])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropKD <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
JunMKD   <- round(cellStats(Jun_CropKD, 'mean'),0)
JunMKDx   <- round(cellStats(Jun_CropKD, 'max'),0)
JunMKDn   <- round(cellStats(Jun_CropKD, 'min'),0)

Jul <- raster(Mean_CLIM[59])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropKD <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
JulMKD   <- round(cellStats(Jul_CropKD, 'mean'),0)
JulMKDx   <- round(cellStats(Jul_CropKD, 'max'),0)
JulMKDn   <- round(cellStats(Jul_CropKD, 'min'),0)

Aug <- raster(Mean_CLIM[55])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropKD <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
AugMKD   <- round(cellStats(Aug_CropKD, 'mean'),0)
AugMKDx   <- round(cellStats(Aug_CropKD, 'max'),0)
AugMKDn   <- round(cellStats(Aug_CropKD, 'min'),0)

Sep <- raster(Mean_CLIM[65])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropKD <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
SepMKD   <- round(cellStats(Sep_CropKD, 'mean'),0)
SepMKDx   <- round(cellStats(Sep_CropKD, 'max'),0)
SepMKDn   <- round(cellStats(Sep_CropKD, 'min'),0)

Oct <- raster(Mean_CLIM[64]) 
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropKD <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
OctMKD   <- round(cellStats(Oct_CropKD, 'mean'),0)
OctMKDx   <- round(cellStats(Oct_CropKD, 'max'),0)
OctMKDn   <- round(cellStats(Oct_CropKD, 'min'),0)

Nov <- raster(Mean_CLIM[63]) 
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropKD <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
NovMKD   <- round(cellStats(Nov_CropKD, 'mean'),0)
NovMKDx   <- round(cellStats(Nov_CropKD, 'max'),0)
NovMKDn   <- round(cellStats(Nov_CropKD, 'min'),0)

Dec <- raster(Mean_CLIM[56]) 
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropKD <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
DecMKD   <- round(cellStats(Dec_CropKD, 'mean'),0)
DecMKDx   <- round(cellStats(Dec_CropKD, 'max'),0)
DecMKDn   <- round(cellStats(Dec_CropKD, 'min'),0)

Ann <- raster(Mean_CLIM[53]) 
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropKD <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
AnnMKD   <- round(cellStats(Ann_CropKD, 'mean'),0)
AnnMKDx   <- round(cellStats(Ann_CropKD, 'max'),0)
AnnMKDn   <- round(cellStats(Ann_CropKD, 'min'),0)

MEANKD2 <- c(JanMKD,FebMKD,MarMKD,AprMKD,MayMKD,JunMKD,JulMKD,AugMKD,SepMKD,OctMKD,NovMKD,DecMKD,AnnMKD)
MEANKD3 <- MEANKD2 * 100

#GKD thresholds for scale Bar
KDUP <- max(JanMKDx,FebMKDx,MarMKDx,AprMKDx,MayMKDx,JunMKDx,JulMKDx,AugMKDx,SepMKDx,OctMKDx,NovMKDx,DecMKDx,AnnMKDx)
KDLO <- min(JanMKDn,FebMKDn,MarMKDn,AprMKDn,MayMKDn,JunMKDn,JulMKDn,AugMKDn,SepMKDn,OctMKDn,NovMKDn,DecMKDn,AnnMKDn)

Cell.CL_Year[9,2:14] <- MEANKD2 


##########   Soil Moisture

Jan <- raster(Mean_CLIM[45])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropSM <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
JanMSM   <- round(cellStats(Jan_CropSM, 'mean'),2)
JanMSMx   <- round(cellStats(Jan_CropSM, 'max'),2)
JanMSMn   <- round(cellStats(Jan_CropSM, 'min'),2)

Feb <- raster(Mean_CLIM[44])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropSM <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
FebMSM   <- round(cellStats(Feb_CropSM, 'mean'),2)
FebMSMx   <- round(cellStats(Feb_CropSM, 'max'),2)
FebMSMn   <- round(cellStats(Feb_CropSM, 'min'),2)

Mar <- raster(Mean_CLIM[48])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropSM <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
MarMSM   <- round(cellStats(Mar_CropSM, 'mean'),2)
MarMSMx   <- round(cellStats(Mar_CropSM, 'max'),2)
MarMSMn   <- round(cellStats(Mar_CropSM, 'min'),2)

Apr <- raster(Mean_CLIM[41])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropSM <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
AprMSM   <- round(cellStats(Apr_CropSM, 'mean'),2)
AprMSMx   <- round(cellStats(Apr_CropSM, 'max'),2)
AprMSMn   <- round(cellStats(Apr_CropSM, 'min'),2)

May <- raster(Mean_CLIM[49])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropSM <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
MayMSM   <- round(cellStats(May_CropSM, 'mean'),2)
MayMSMx   <- round(cellStats(May_CropSM, 'max'),2)
MayMSMn   <- round(cellStats(May_CropSM, 'min'),2)

Jun <- raster(Mean_CLIM[47])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropSM <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
JunMSM   <- round(cellStats(Jun_CropSM, 'mean'),2)
JunMSMx   <- round(cellStats(Jun_CropSM, 'max'),2)
JunMSMn   <- round(cellStats(Jun_CropSM, 'min'),2)

Jul <- raster(Mean_CLIM[46])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropSM <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
JulMSM   <- round(cellStats(Jul_CropSM, 'mean'),2)
JulMSMx   <- round(cellStats(Jul_CropSM, 'max'),2)
JulMSMn   <- round(cellStats(Jul_CropSM, 'min'),2)

Aug <- raster(Mean_CLIM[42])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropSM <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
AugMSM   <- round(cellStats(Aug_CropSM, 'mean'),2)
AugMSMx   <- round(cellStats(Aug_CropSM, 'max'),2)
AugMSMn   <- round(cellStats(Aug_CropSM, 'min'),2)

Sep <- raster(Mean_CLIM[52])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropSM <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
SepMSM   <- round(cellStats(Sep_CropSM, 'mean'),2)
SepMSMx   <- round(cellStats(Sep_CropSM, 'max'),2)
SepMSMn   <- round(cellStats(Sep_CropSM, 'min'),2)

Oct <- raster(Mean_CLIM[51]) 
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropSM <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
OctMSM   <- round(cellStats(Oct_CropSM, 'mean'),2)
OctMSMx   <- round(cellStats(Oct_CropSM, 'max'),2)
OctMSMn   <- round(cellStats(Oct_CropSM, 'min'),2)

Nov <- raster(Mean_CLIM[50]) 
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropSM <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
NovMSM   <- round(cellStats(Nov_CropSM, 'mean'),2)
NovMSMx   <- round(cellStats(Nov_CropSM, 'max'),2)
NovMSMn   <- round(cellStats(Nov_CropSM, 'min'),2)

Dec <- raster(Mean_CLIM[43]) 
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropSM <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
DecMSM   <- round(cellStats(Dec_CropSM, 'mean'),2)
DecMSMx   <- round(cellStats(Dec_CropSM, 'max'),2)
DecMSMn   <- round(cellStats(Dec_CropSM, 'min'),2)

Ann <- raster(Mean_CLIM[40]) 
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropSM <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
AnnMSM   <- round(cellStats(Ann_CropSM, 'mean'),2)
AnnMSMx   <- round(cellStats(Ann_CropSM, 'max'),2)
AnnMSMn   <- round(cellStats(Ann_CropSM, 'min'),2)

MEANSM2 <- c(JanMSM,FebMSM,MarMSM,AprMSM,MayMSM,JunMSM,JulMSM,AugMSM,SepMSM,OctMSM,NovMSM,DecMSM,AnnMSM)
MEANSM3 <- MEANSM2 * 100

#GSM thresholds for scale Bar
SMUP <- max(JanMSMx,FebMSMx,MarMSMx,AprMSMx,MayMSMx,JunMSMx,JulMSMx,AugMSMx,SepMSMx,OctMSMx,NovMSMx,DecMSMx,AnnMSMx)
SMLo <- min(JanMSMn,FebMSMn,MarMSMn,AprMSMn,MayMSMn,JunMSMn,JulMSMn,AugMSMn,SepMSMn,OctMSMn,NovMSMn,DecMSMn,AnnMSMn)

Cell.CL_Year[8,2:14] <- MEANSM3




##########   EVAPOTRANSIPRATION

Jan <- raster(Mean_CLIM[19])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropET <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Jan_CropET <- Jan_CropET  * 0.0393701}
JanMET   <- round(cellStats(Jan_CropET, 'mean'),0)
JanMETx   <- round(cellStats(Jan_CropET, 'max'),0)
JanMETn   <- round(cellStats(Jan_CropET, 'min'),0)

Feb <- raster(Mean_CLIM[18])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropET <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Feb_CropET <- Feb_CropET  * 0.0393701}
FebMET   <- round(cellStats(Feb_CropET, 'mean'),0)
FebMETx   <- round(cellStats(Feb_CropET, 'max'),0)
FebMETn   <- round(cellStats(Feb_CropET, 'min'),0)

Mar <- raster(Mean_CLIM[22])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropET <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Mar_CropET <- Mar_CropET  * 0.0393701}
MarMET   <- round(cellStats(Mar_CropET, 'mean'),0)
MarMETx   <- round(cellStats(Mar_CropET, 'max'),0)
MarMETn   <- round(cellStats(Mar_CropET, 'min'),0)

Apr <- raster(Mean_CLIM[15])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropET <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Apr_CropET <- Apr_CropET  * 0.0393701}
AprMET   <- round(cellStats(Apr_CropET, 'mean'),0)
AprMETx   <- round(cellStats(Apr_CropET, 'max'),0)
AprMETn   <- round(cellStats(Apr_CropET, 'min'),0)

May <- raster(Mean_CLIM[23])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropET <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {May_CropET <- May_CropET  * 0.0393701}
MayMET   <- round(cellStats(May_CropET, 'mean'),0)
MayMETx   <- round(cellStats(May_CropET, 'max'),0)
MayMETn   <- round(cellStats(May_CropET, 'min'),0)

Jun <- raster(Mean_CLIM[21])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropET <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Jun_CropET <- Jun_CropET  * 0.0393701}
JunMET   <- round(cellStats(Jun_CropET, 'mean'),0)
JunMETx   <- round(cellStats(Jun_CropET, 'max'),0)
JunMETn   <- round(cellStats(Jun_CropET, 'min'),0)

Jul <- raster(Mean_CLIM[20])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropET <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Jul_CropET <- Jul_CropET  * 0.0393701}
JulMET   <- round(cellStats(Jul_CropET, 'mean'),0)
JulMETx   <- round(cellStats(Jul_CropET, 'max'),0)
JulMETn   <- round(cellStats(Jul_CropET, 'min'),0)

Aug <- raster(Mean_CLIM[16])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropET <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Aug_CropET <- Aug_CropET  * 0.0393701}
AugMET   <- round(cellStats(Aug_CropET, 'mean'),0)
AugMETx   <- round(cellStats(Aug_CropET, 'max'),0)
AugMETn   <- round(cellStats(Aug_CropET, 'min'),0)

Sep <- raster(Mean_CLIM[26])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropET <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Sep_CropET <- Sep_CropET  * 0.0393701}
SepMET   <- round(cellStats(Sep_CropET, 'mean'),0)
SepMETx   <- round(cellStats(Sep_CropET, 'max'),0)
SepMETn   <- round(cellStats(Sep_CropET, 'min'),0)

Oct <- raster(Mean_CLIM[25]) 
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropET <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Oct_CropET <- Oct_CropET  * 0.0393701}
OctMET   <- round(cellStats(Oct_CropET, 'mean'),0)
OctMETx   <- round(cellStats(Oct_CropET, 'max'),0)
OctMETn   <- round(cellStats(Oct_CropET, 'min'),0)

Nov <- raster(Mean_CLIM[24]) 
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropET <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Nov_CropET <- Nov_CropET  * 0.0393701}
NovMET   <- round(cellStats(Nov_CropET, 'mean'),0)
NovMETx   <- round(cellStats(Nov_CropET, 'max'),0)
NovMETn   <- round(cellStats(Nov_CropET, 'min'),0)

Dec <- raster(Mean_CLIM[17]) 
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropET <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Dec_CropET <- Dec_CropET  * 0.0393701}
DecMET   <- round(cellStats(Dec_CropET, 'mean'),0)
DecMETx   <- round(cellStats(Dec_CropET, 'max'),0)
DecMETn   <- round(cellStats(Dec_CropET, 'min'),0)

Ann <- raster(Mean_CLIM[14]) 
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropET <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Ann_CropET <- Ann_CropET  * 0.0393701}
AnnMET   <- round(cellStats(Ann_CropET, 'mean'),0)
AnnMETx   <- round(cellStats(Ann_CropET, 'max'),0)
AnnMETn   <- round(cellStats(Ann_CropET, 'min'),0)

MEANET2 <- c(JanMET,FebMET,MarMET,AprMET,MayMET,JunMET,JulMET,AugMET,SepMET,OctMET,NovMET,DecMET,AnnMET)


#Get thresholds for scale Bar
ETUP <- max(JanMETx,FebMETx,MarMETx,AprMETx,MayMETx,JunMETx,JulMETx,AugMETx,SepMETx,OctMETx,NovMETx,DecMETx,AnnMETx)
ETLo <- min(JanMETn,FebMETn,MarMETn,AprMETn,MayMETn,JunMETn,JulMETn,AugMETn,SepMETn,OctMETn,NovMETn,DecMETn,AnnMETn)

Cell.CL_Year[7,2:14] <- MEANET2 

plot(Aug)
##########   Cloud Frequency

Jan <- raster(Mean_CLIM[6])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropCF <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
JanMCF   <- round(cellStats(Jan_CropCF, 'mean'),2)
JanMCFx   <- round(cellStats(Jan_CropCF, 'max'),2)
JanMCFn   <- round(cellStats(Jan_CropCF, 'min'),2)

Feb <- raster(Mean_CLIM[5])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropCF <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
FebMCF   <- round(cellStats(Feb_CropCF, 'mean'),2)
FebMCFx   <- round(cellStats(Feb_CropCF, 'max'),2)
FebMCFn   <- round(cellStats(Feb_CropCF, 'min'),2)

Mar <- raster(Mean_CLIM[9])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropCF <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
MarMCF   <- round(cellStats(Mar_CropCF, 'mean'),2)
MarMCFx   <- round(cellStats(Mar_CropCF, 'max'),2)
MarMCFn   <- round(cellStats(Mar_CropCF, 'min'),2)

Apr <- raster(Mean_CLIM[2])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropCF <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
AprMCF   <- round(cellStats(Apr_CropCF, 'mean'),2)
AprMCFx   <- round(cellStats(Apr_CropCF, 'max'),2)
AprMCFn   <- round(cellStats(Apr_CropCF, 'min'),2)

May <- raster(Mean_CLIM[10])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropCF <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
MayMCF   <- round(cellStats(May_CropCF, 'mean'),2)
MayMCFx   <- round(cellStats(May_CropCF, 'max'),2)
MayMCFn   <- round(cellStats(May_CropCF, 'min'),2)

Jun <- raster(Mean_CLIM[8])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropCF <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
JunMCF   <- round(cellStats(Jun_CropCF, 'mean'),2)
JunMCFx   <- round(cellStats(Jun_CropCF, 'max'),2)
JunMCFn   <- round(cellStats(Jun_CropCF, 'min'),2)

Jul <- raster(Mean_CLIM[7])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropCF <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
JulMCF   <- round(cellStats(Jul_CropCF, 'mean'),2)
JulMCFx   <- round(cellStats(Jul_CropCF, 'max'),2)
JulMCFn   <- round(cellStats(Jul_CropCF, 'min'),2)

Aug <- raster(Mean_CLIM[3])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropCF <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
AugMCF   <- round(cellStats(Aug_CropCF, 'mean'),2)
AugMCFx   <- round(cellStats(Aug_CropCF, 'max'),2)
AugMCFn   <- round(cellStats(Aug_CropCF, 'min'),2)

Sep <- raster(Mean_CLIM[13])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropCF <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
SepMCF   <- round(cellStats(Sep_CropCF, 'mean'),2)
SepMCFx   <- round(cellStats(Sep_CropCF, 'max'),2)
SepMCFn   <- round(cellStats(Sep_CropCF, 'min'),2)

Oct <- raster(Mean_CLIM[12]) 
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropCF <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
OctMCF   <- round(cellStats(Oct_CropCF, 'mean'),2)
OctMCFx   <- round(cellStats(Oct_CropCF, 'max'),2)
OctMCFn   <- round(cellStats(Oct_CropCF, 'min'),2)

Nov <- raster(Mean_CLIM[11]) 
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropCF <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
NovMCF   <- round(cellStats(Nov_CropCF, 'mean'),2)
NovMCFx   <- round(cellStats(Nov_CropCF, 'max'),2)
NovMCFn   <- round(cellStats(Nov_CropCF, 'min'),2)

Dec <- raster(Mean_CLIM[4]) 
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropCF <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
DecMCF   <- round(cellStats(Dec_CropCF, 'mean'),2)
DecMCFx   <- round(cellStats(Dec_CropCF, 'max'),2)
DecMCFn   <- round(cellStats(Dec_CropCF, 'min'),2)

Ann <- raster(Mean_CLIM[1]) 
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropCF <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
AnnMCF   <- round(cellStats(Ann_CropCF, 'mean'),2)
AnnMCFx   <- round(cellStats(Ann_CropCF, 'max'),2)
AnnMCFn   <- round(cellStats(Ann_CropCF, 'min'),2)

MEANCF2 <- c(JanMCF,FebMCF,MarMCF,AprMCF,MayMCF,JunMCF,JulMCF,AugMCF,SepMCF,OctMCF,NovMCF,DecMCF,AnnMCF)
MEANCF3 <- MEANCF2 * 100

#Get thresholds for scale Bar
CFUP <- max(JanMCFx,FebMCFx,MarMCFx,AprMCFx,MayMCFx,JunMCFx,JulMCFx,AugMCFx,SepMCFx,OctMCFx,NovMCFx,DecMCFx,AnnMCFx)
CFLO <- min(JanMCFn,FebMCFn,MarMCFn,AprMCFn,MayMCFn,JunMCFn,JulMCFn,AugMCFn,SepMCFn,OctMCFn,NovMCFn,DecMCFn,AnnMCFn)

Cell.CL_Year[6,2:14] <- MEANCF3

##########   

##########   MEAN TEMPERATURE 

Jan <- raster(Mean_CLIM[71])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropTA <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Jan_CropTA <- (Jan_CropTA * 1.8) + 32}
JanMTA   <- round(cellStats(Jan_CropTA, 'mean'),1)
JanMTAx   <- round(cellStats(Jan_CropTA, 'max'),1)
JanMTAn   <- round(cellStats(Jan_CropTA, 'min'),1)

Feb <- raster(Mean_CLIM[70])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropTA <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Feb_CropTA <- (Feb_CropTA * 1.8) + 32}
FebMTA   <- round(cellStats(Feb_CropTA, 'mean'),1)
FebMTAx   <- round(cellStats(Feb_CropTA, 'max'),1)
FebMTAn   <- round(cellStats(Feb_CropTA, 'min'),1)

Mar <- raster(Mean_CLIM[74])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropTA <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Mar_CropTA <- (Mar_CropTA * 1.8) + 32}
MarMTA   <- round(cellStats(Mar_CropTA, 'mean'),1)
MarMTAx   <- round(cellStats(Mar_CropTA, 'max'),1)
MarMTAn   <- round(cellStats(Mar_CropTA, 'min'),1)

Apr <- raster(Mean_CLIM[67])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropTA <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Apr_CropTA <- (Apr_CropTA * 1.8) + 32}
AprMTA   <- round(cellStats(Apr_CropTA, 'mean'),1)
AprMTAx   <- round(cellStats(Apr_CropTA, 'max'),1)
AprMTAn   <- round(cellStats(Apr_CropTA, 'min'),1)

May <- raster(Mean_CLIM[75])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropTA <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {May_CropTA <- (May_CropTA * 1.8) + 32}
MayMTA   <- round(cellStats(May_CropTA, 'mean'),1)
MayMTAx   <- round(cellStats(May_CropTA, 'max'),1)
MayMTAn   <- round(cellStats(May_CropTA, 'min'),1)

Jun <- raster(Mean_CLIM[73])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropTA <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Jun_CropTA <- (Jun_CropTA * 1.8) + 32}
JunMTA   <- round(cellStats(Jun_CropTA, 'mean'),1)
JunMTAx   <- round(cellStats(Jun_CropTA, 'max'),1)
JunMTAn   <- round(cellStats(Jun_CropTA, 'min'),1)

Jul <- raster(Mean_CLIM[72])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropTA <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Jul_CropTA <- (Jul_CropTA * 1.8) + 32}
JulMTA   <- round(cellStats(Jul_CropTA, 'mean'),1)
JulMTAx   <- round(cellStats(Jul_CropTA, 'max'),1)
JulMTAn   <- round(cellStats(Jul_CropTA, 'min'),1)

Aug <- raster(Mean_CLIM[68])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropTA <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Aug_CropTA <- (Aug_CropTA * 1.8) + 32}
AugMTA   <- round(cellStats(Aug_CropTA, 'mean'),1)
AugMTAx   <- round(cellStats(Aug_CropTA, 'max'),1)
AugMTAn   <- round(cellStats(Aug_CropTA, 'min'),1)

Sep <- raster(Mean_CLIM[78])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropTA <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Sep_CropTA <- (Sep_CropTA * 1.8) + 32}
SepMTA   <- round(cellStats(Sep_CropTA, 'mean'),1)
SepMTAx   <- round(cellStats(Sep_CropTA, 'max'),1)
SepMTAn   <- round(cellStats(Sep_CropTA, 'min'),1)

Oct <- raster(Mean_CLIM[77]) 
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropTA <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Oct_CropTA <- (Oct_CropTA * 1.8) + 32}
OctMTA   <- round(cellStats(Oct_CropTA, 'mean'),1)
OctMTAx   <- round(cellStats(Oct_CropTA, 'max'),1)
OctMTAn   <- round(cellStats(Oct_CropTA, 'min'),1)

Nov <- raster(Mean_CLIM[76]) 
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropTA <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Nov_CropTA <- (Nov_CropTA * 1.8) + 32}
NovMTA   <- round(cellStats(Nov_CropTA, 'mean'),1)
NovMTAx   <- round(cellStats(Nov_CropTA, 'max'),1)
NovMTAn   <- round(cellStats(Nov_CropTA, 'min'),1)

Dec <- raster(Mean_CLIM[69]) 
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropTA <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Dec_CropTA <- (Dec_CropTA * 1.8) + 32}
DecMTA   <- round(cellStats(Dec_CropTA, 'mean'),1)
DecMTAx   <- round(cellStats(Dec_CropTA, 'max'),1)
DecMTAn   <- round(cellStats(Dec_CropTA, 'min'),1)

Ann <- raster(Mean_CLIM[66]) 
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropTA <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Ann_CropTA <- (Ann_CropTA * 1.8) + 32}
AnnMTA   <- round(cellStats(Ann_CropTA, 'mean'),1)
AnnMTAx   <- round(cellStats(Ann_CropTA, 'max'),1)
AnnMTAn   <- round(cellStats(Ann_CropTA, 'min'),1)

MEANTA2 <- c(JanMTA,FebMTA,MarMTA,AprMTA,MayMTA,JunMTA,JulMTA,AugMTA,SepMTA,OctMTA,NovMTA,DecMTA,AnnMTA)

#Get thresholds for scale Bar
TaUP <- max(JanMTAx,FebMTAx,MarMTAx,AprMTAx,MayMTAx,JunMTAx,JulMTAx,AugMTAx,SepMTAx,OctMTAx,NovMTAx,DecMTAx,AnnMTAx)
TaLo <- min(JanMTAn,FebMTAn,MarMTAn,AprMTAn,MayMTAn,JunMTAn,JulMTAn,AugMTAn,SepMTAn,OctMTAn,NovMTAn,DecMTAn,AnnMTAn)

Cell.CL_Year[3,2:14] <- MEANTA2 

########## MAX TEMPERATURE 

Jan <- raster(Mean_CLIM[84])
crs(Jan) <- crs(EXAMP)
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropTX <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Jan_CropTX <- (Jan_CropTX * 1.8) + 32}
JanMTX  <- round(cellStats(Jan_CropTX , 'mean'),1)
JanMTXx   <- round(cellStats(Jan_CropTX, 'max'),1)
JanMTXn   <- round(cellStats(Jan_CropTX, 'min'),1)

Feb <- raster(Mean_CLIM[83])
crs(Feb) <- crs(EXAMP)
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropTX  <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Feb_CropTX <- (Feb_CropTX * 1.8) + 32}
FebMTX   <- round(cellStats(Feb_CropTX , 'mean'),1)
FebMTXx   <- round(cellStats(Feb_CropTX, 'max'),1)
FebMTXn   <- round(cellStats(Feb_CropTX, 'min'),1)

Mar <- raster(Mean_CLIM[87])
crs(Mar) <- crs(EXAMP)
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropTX  <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Mar_CropTX <- (Mar_CropTX * 1.8) + 32}
MarMTX   <- round(cellStats(Mar_CropTX , 'mean'),1)
MarMTXx   <- round(cellStats(Mar_CropTX, 'max'),1)
MarMTXn   <- round(cellStats(Mar_CropTX, 'min'),1)

Apr <- raster(Mean_CLIM[80])
crs(Apr) <- crs(EXAMP)
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropTX  <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Apr_CropTX <- (Apr_CropTX * 1.8) + 32}
AprMTX   <- round(cellStats(Apr_CropTX , 'mean'),1)
AprMTXx   <- round(cellStats(Apr_CropTX, 'max'),1)
AprMTXn   <- round(cellStats(Apr_CropTX, 'min'),1)

May <- raster(Mean_CLIM[88])
crs(May) <- crs(EXAMP)
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropTX  <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {May_CropTX <- (May_CropTX * 1.8) + 32}
MayMTX    <- round(cellStats(May_CropTX, 'mean'),1)
MayMTXx   <- round(cellStats(May_CropTX, 'max'),1)
MayMTXn   <- round(cellStats(May_CropTX, 'min'),1)

Jun <- raster(Mean_CLIM[86])
crs(Jun) <- crs(EXAMP)
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropTX  <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Jun_CropTX <- (Jun_CropTX * 1.8) + 32}
JunMTX   <- round(cellStats(Jun_CropTX , 'mean'),1)
JunMTXx   <- round(cellStats(Jun_CropTX, 'max'),1)
JunMTXn   <- round(cellStats(Jun_CropTX, 'min'),1)

Jul <- raster(Mean_CLIM[85])
crs(Jul) <- crs(EXAMP)
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropTX  <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Jul_CropTX <- (Jul_CropTX * 1.8) + 32}
JulMTX    <- round(cellStats(Jul_CropTX, 'mean'),1)
JulMTXx   <- round(cellStats(Jul_CropTX, 'max'),1)
JulMTXn   <- round(cellStats(Jul_CropTX, 'min'),1)

Aug <- raster(Mean_CLIM[81])
crs(Aug) <- crs(EXAMP)
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropTX  <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Aug_CropTX <- (Aug_CropTX * 1.8) + 32}
AugMTX   <- round(cellStats(Aug_CropTX, 'mean'),1)
AugMTXx   <- round(cellStats(Aug_CropTX, 'max'),1)
AugMTXn   <- round(cellStats(Aug_CropTX, 'min'),1)

Sep <- raster(Mean_CLIM[91])
crs(Sep) <- crs(EXAMP)
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropTX  <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Sep_CropTX <- (Sep_CropTX * 1.8) + 32}
SepMTX   <- round(cellStats(Sep_CropTX, 'mean'),1)
SepMTXx   <- round(cellStats(Sep_CropTX, 'max'),1)
SepMTXn   <- round(cellStats(Sep_CropTX, 'min'),1)

Oct <- raster(Mean_CLIM[90]) 
crs(Oct) <- crs(EXAMP)
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropTX  <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Oct_CropTX <- (Oct_CropTX * 1.8) + 32}
OctMTX   <- round(cellStats(Oct_CropTX, 'mean'),1)
OctMTXx   <- round(cellStats(Oct_CropTX, 'max'),1)
OctMTXn   <- round(cellStats(Oct_CropTX, 'min'),1)

Nov <- raster(Mean_CLIM[89]) 
crs(Nov) <- crs(EXAMP)
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropTX  <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Nov_CropTX <- (Nov_CropTX * 1.8) + 32}
NovMTX   <- round(cellStats(Nov_CropTX, 'mean'),1)
NovMTXx   <- round(cellStats(Nov_CropTX, 'max'),1)
NovMTXn   <- round(cellStats(Nov_CropTX, 'min'),1)

Dec <- raster(Mean_CLIM[82]) 
crs(Dec) <- crs(EXAMP)
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropTX  <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Dec_CropTX <- (Dec_CropTX * 1.8) + 32}
DecMTX   <- round(cellStats(Dec_CropTX, 'mean'),1)
DecMTXx   <- round(cellStats(Dec_CropTX, 'max'),1)
DecMTXn   <- round(cellStats(Dec_CropTX, 'min'),1)

Ann <- raster(Mean_CLIM[79]) 
crs(Ann) <- crs(EXAMP)
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropTX <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Ann_CropTX <- (Ann_CropTX * 1.8) + 32}
AnnMTX   <- round(cellStats(Ann_CropTX, 'mean'),1)
AnnMTXx   <- round(cellStats(Ann_CropTX, 'max'),1)
AnnMTXn   <- round(cellStats(Ann_CropTX, 'min'),1)

MEANTX2 <- c(JanMTX,FebMTX,MarMTX,AprMTX,MayMTX,JunMTX,JulMTX,AugMTX,SepMTX,OctMTX,NovMTX,DecMTX,AnnMTX)

##########   Add to dataframe

Cell.CL_Year[4,2:14] <- MEANTX2

##########   Get thresholds for figures 

TxUP <- max(JanMTXx,FebMTXx,MarMTXx,AprMTXx,MayMTXx,JunMTXx,JulMTXx,AugMTXx,SepMTXx,OctMTXx,NovMTXx,DecMTXx,AnnMTXx)
TxLO <- min(JanMTXn,FebMTXn,MarMTXn,AprMTXn,MayMTXn,JunMTXn,JulMTXn,AugMTXn,SepMTXn,OctMTXn,NovMTXn,DecMTXn,AnnMTXn)

##########    MIN TEMPERATURE 

Jan <- raster(Mean_CLIM[97])
crs(Jan) <- crs(EXAMP)
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropTN <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Jan_CropTN <- (Jan_CropTN * 1.8) + 32}
JanMTN   <- round(cellStats(Jan_CropTN, 'mean'),1)
JanMTNx   <- round(cellStats(Jan_CropTN , 'max'),1)
JanMTNn   <- round(cellStats(Jan_CropTN , 'min'),1)

Feb <- raster(Mean_CLIM[96])
crs(Feb) <- crs(EXAMP)
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropTN <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Feb_CropTN <- (Feb_CropTN * 1.8) + 32}
FebMTN  <- round(cellStats(Feb_CropTN, 'mean'),1)
FebMTNx   <- round(cellStats(Feb_CropTN , 'max'),1)
FebMTNn   <- round(cellStats(Feb_CropTN , 'min'),1)

Mar <- raster(Mean_CLIM[100])
crs(Mar) <- crs(EXAMP)
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropTN <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Mar_CropTN <- (Mar_CropTN * 1.8) + 32}
MarMTN   <- round(cellStats(Mar_CropTN, 'mean'),1)
MarMTNx   <- round(cellStats(Mar_CropTN , 'max'),1)
MarMTNn   <- round(cellStats(Mar_CropTN , 'min'),1)

Apr <- raster(Mean_CLIM[93])
crs(Apr) <- crs(EXAMP)
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropTN <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Apr_CropTN <- (Apr_CropTN * 1.8) + 32}
AprMTN   <- round(cellStats(Apr_CropTN, 'mean'),1)
AprMTNx   <- round(cellStats(Apr_CropTN , 'max'),1)
AprMTNn   <- round(cellStats(Apr_CropTN , 'min'),1)

May <- raster(Mean_CLIM[101])
crs(May) <- crs(EXAMP)
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropTN <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {May_CropTN <- (May_CropTN * 1.8) + 32}
MayMTN   <- round(cellStats(May_CropTN, 'mean'),1)
MayMTNx   <- round(cellStats(May_CropTN , 'max'),1)
MayMTNn   <- round(cellStats(May_CropTN , 'min'),1)

Jun <- raster(Mean_CLIM[99])
crs(Jun) <- crs(EXAMP)
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropTN <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Jun_CropTN <- (Jun_CropTN * 1.8) + 32}
JunMTN   <- round(cellStats(Jun_CropTN, 'mean'),1)
JunMTNx   <- round(cellStats(Jun_CropTN , 'max'),1)
JunMTNn   <- round(cellStats(Jun_CropTN , 'min'),1)

Jul <- raster(Mean_CLIM[98])
crs(Jul) <- crs(EXAMP)
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropTN <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Jul_CropTN <- (Jul_CropTN * 1.8) + 32}
JulMTN   <- round(cellStats(Jul_CropTN, 'mean'),1)
JulMTNx   <- round(cellStats(Jul_CropTN , 'max'),1)
JulMTNn   <- round(cellStats(Jul_CropTN , 'min'),1)

Aug <- raster(Mean_CLIM[94])
crs(Aug) <- crs(EXAMP)
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropTN <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Aug_CropTN <- (Aug_CropTN * 1.8) + 32}
AugMTN   <- round(cellStats(Aug_CropTN, 'mean'),1)
AugMTNx   <- round(cellStats(Aug_CropTN , 'max'),1)
AugMTNn   <- round(cellStats(Aug_CropTN , 'min'),1)

Sep <- raster(Mean_CLIM[104])
crs(Sep) <- crs(EXAMP)
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropTN <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Sep_CropTN <- (Sep_CropTN * 1.8) + 32}
SepMTN   <- round(cellStats(Sep_CropTN, 'mean'),1)
SepMTNx   <- round(cellStats(Sep_CropTN , 'max'),1)
SepMTNn   <- round(cellStats(Sep_CropTN , 'min'),1)

Oct <- raster(Mean_CLIM[103]) 
crs(Oct) <- crs(EXAMP)
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropTN <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Oct_CropTN <- (Oct_CropTN * 1.8) + 32}
OctMTN   <- round(cellStats(Oct_CropTN, 'mean'),1)
OctMTNx   <- round(cellStats(Oct_CropTN , 'max'),1)
OctMTNn   <- round(cellStats(Oct_CropTN , 'min'),1)

Nov <- raster(Mean_CLIM[102]) 
crs(Nov) <- crs(EXAMP)
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropTN <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Nov_CropTN <- (Nov_CropTN * 1.8) + 32}
NovMTN   <- round(cellStats(Nov_CropTN, 'mean'),1)
NovMTNx   <- round(cellStats(Nov_CropTN , 'max'),1)
NovMTNn   <- round(cellStats(Nov_CropTN , 'min'),1)

Dec <- raster(Mean_CLIM[95])
crs(Dec) <- crs(EXAMP)
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropTN <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Dec_CropTN <- (Dec_CropTN * 1.8) + 32}
DecMTN   <- round(cellStats(Dec_CropTN, 'mean'),1)
DecMTNx   <- round(cellStats(Dec_CropTN , 'max'),1)
DecMTNn   <- round(cellStats(Dec_CropTN , 'min'),1)

Ann <- raster(Mean_CLIM[92]) 
crs(Ann) <- crs(EXAMP)
Ann_Mask <- mask(x = Ann, mask = UNIT_X[[u]])
Ann_CropTN <- crop(x = Ann_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Ann_CropTN <- (Ann_CropTN * 1.8) + 32}
AnnMTN   <- round(cellStats(Ann_CropTN, 'mean'),1)
AnnMTNx   <- round(cellStats(Ann_CropTN , 'max'),1)
AnnMTNn   <- round(cellStats(Ann_CropTN , 'min'),1)

########## Add to Table 

MEANTN2 <- c(JanMTN,FebMTN,MarMTN,AprMTN,MayMTN,JunMTN,JulMTN,AugMTN,SepMTN,OctMTN,NovMTN,DecMTN,AnnMTN)
Cell.CL_Year[2,2:14] <- MEANTN2

########## Get thresholds for figures 

TnUP <- max(JanMTNx,FebMTNx,MarMTNx,AprMTNx,MayMTNx,JunMTNx,JulMTNx,AugMTNx,SepMTNx,OctMTNx,NovMTNx,DecMTNx,AnnMTNx)
TnLO <- min(JanMTNn,FebMTNn,MarMTNn,AprMTNn,MayMTNn,JunMTNn,JulMTNn,AugMTNn,SepMTNn,OctMTNn,NovMTNn,DecMTNn,AnnMTNn)

##########   RELATIVE HuMIDITY

Jan <- raster(Mean_CLIM[32])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropRH <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
JanMRH    <- round(cellStats(Jan_CropRH , 'mean'),0)
JanMRHx   <- round(cellStats(Jan_CropRH  , 'max'),0)
JanMRHn   <- round(cellStats(Jan_CropRH  , 'min'),0)

Feb <- raster(Mean_CLIM[31])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropRH  <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
FebMRH    <- round(cellStats(Feb_CropRH , 'mean'),0)
FebMRHx   <- round(cellStats(Feb_CropRH  , 'max'),0)
FebMRHn   <- round(cellStats(Feb_CropRH  , 'min'),0)

Mar <- raster(Mean_CLIM[35])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropRH  <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
MarMRH    <- round(cellStats(Mar_CropRH , 'mean'),0)
MarMRHx   <- round(cellStats(Mar_CropRH  , 'max'),0)
MarMRHn   <- round(cellStats(Mar_CropRH  , 'min'),0)

Apr <- raster(Mean_CLIM[28])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropRH  <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
AprMRH    <- round(cellStats(Apr_CropRH , 'mean'),0)
AprMRHx   <- round(cellStats(Apr_CropRH  , 'max'),0)
AprMRHn   <- round(cellStats(Apr_CropRH  , 'min'),0)

May <- raster(Mean_CLIM[36])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropRH  <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
MayMRH   <- round(cellStats(May_CropRH, 'mean'),0)
MayMRHx   <- round(cellStats(May_CropRH  , 'max'),0)
MayMRHn   <- round(cellStats(May_CropRH  , 'min'),0)

Jun <- raster(Mean_CLIM[34])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropRH  <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
JunMRH    <- round(cellStats(Jun_CropRH, 'mean'),0)
JunMRHx   <- round(cellStats(Jun_CropRH  , 'max'),0)
JunMRHn   <- round(cellStats(Jun_CropRH  , 'min'),0)

Jul <- raster(Mean_CLIM[33])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropRH  <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
JulMRH   <- round(cellStats(Jul_CropRH, 'mean'),0)
JulMRHx   <- round(cellStats(Jul_CropRH  , 'max'),0)
JulMRHn   <- round(cellStats(Jul_CropRH  , 'min'),0)

Aug <- raster(Mean_CLIM[29])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropRH  <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
AugMRH   <- round(cellStats(Aug_CropRH, 'mean'),0)
AugMRHx   <- round(cellStats(Aug_CropRH  , 'max'),0)
AugMRHn   <- round(cellStats(Aug_CropRH  , 'min'),0)

Sep <- raster(Mean_CLIM[39])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropRH  <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
SepMRH <- round(cellStats(Sep_CropRH, 'mean'),0)
SepMRHx   <- round(cellStats(Sep_CropRH  , 'max'),0)
SepMRHn   <- round(cellStats(Sep_CropRH  , 'min'),0)

Oct <- raster(Mean_CLIM[38]) 
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropRH  <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
OctMRH <- round(cellStats(Oct_CropRH, 'mean'),0)
OctMRHx   <- round(cellStats(Oct_CropRH  , 'max'),0)
OctMRHn   <- round(cellStats(Oct_CropRH  , 'min'),0)

Nov <- raster(Mean_CLIM[37]) 
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropRH  <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
NovMRH <- round(cellStats(Nov_CropRH , 'mean'),0)
NovMRHx   <- round(cellStats(Nov_CropRH  , 'max'),0)
NovMRHn   <- round(cellStats(Nov_CropRH  , 'min'),0)

Dec <- raster(Mean_CLIM[30]) 
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropRH  <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
DecMRH    <- round(cellStats(Dec_CropRH, 'mean'),0)
DecMRHx   <- round(cellStats(Dec_CropRH  , 'max'),0)
DecMRHn   <- round(cellStats(Dec_CropRH  , 'min'),0)

ANN <- raster(Mean_CLIM[27]) 
ANN_Mask <- mask(x = ANN, mask = UNIT_X[[u]])
ANN_CropRH  <- crop(x = ANN_Mask, y = extent(UNIT_X[[u]]))
ANNMRH    <- round(cellStats(ANN_CropRH, 'mean'),0)
ANNMRHx   <- round(cellStats(ANN_CropRH  , 'max'),0)
ANNMRHn   <- round(cellStats(ANN_CropRH  , 'min'),0)

##########   Add to Table 

MEANRH2 <- c(JanMRH,FebMRH ,MarMRH ,AprMRH ,MayMRH ,JunMRH ,JulMRH ,AugMRH ,SepMRH ,OctMRH ,NovMRH ,DecMRH, ANNMRH  )
Cell.CL_Year[5,2:14] <- MEANRH2 

##########   Get thresholds for figures 

RHUP <- max(JanMRHx,FebMRHx ,MarMRHx ,AprMRHx ,MayMRHx ,JunMRHx ,JulMRHx ,AugMRHx ,SepMRHx ,OctMRHx ,NovMRHx ,DecMRHx, ANNMRHx  )
RHLO <- min(JanMRHn,FebMRHn ,MarMRHn ,AprMRHn ,MayMRHn ,JunMRHn ,JulMRHn ,AugMRHn ,SepMRHn ,OctMRHn ,NovMRHn ,DecMRHn, ANNMRHn  )




#########   OTHER CLIMATE VARIABLES
### Read in the Annual rasters for each variable only and extract min and max to show spatial range. 
# Have extracted the min max from Annual layers already This next step can be streamlined 

#See Varible "ANN_CropRH" for example 
RH_P <- raster(Mean_CLIM[27])
RH_P_Mask <- mask(x = RH_P, mask = UNIT_X[[u]])
RH_P_Crop <- crop(x = RH_P_Mask, y = extent(UNIT_X[[u]]))
RH_P_M   <- round(cellStats(RH_P_Crop, 'mean'),0)
RHUPA   <- round(cellStats(RH_P_Crop  , 'max'),1)
RHLOA   <- round(cellStats(RH_P_Crop  , 'min'),1)

Tair_P <- raster(Mean_CLIM[66])
names(Tair_P) = "Mean Air Temp."
Tair_P_Mask <- mask(x = Tair_P, mask = UNIT_X[[u]])
Tair_P_Crop <- crop(x = Tair_P_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") {Tair_P_Crop <- (Tair_P_Crop * 1.8) + 32}
Tair_P_M   <- round(cellStats(Tair_P_Crop,  'mean'),1)
Tair_PName <- names(Tair_P)

SM_P <- raster(Mean_CLIM[40])
names(SM_P) = "Soil Moisture"
SM_P_Mask <- mask(x = SM_P, mask = UNIT_X[[u]])
SM_P_Crop <- crop(x = SM_P_Mask, y = extent(UNIT_X[[u]]))
SM_P_M   <- round(cellStats(SM_P_Crop, 'mean'),2)
SM_PName <- names(SM_P)
SMUP   <- round(cellStats(SM_P_Crop  , 'max'),2)
SMLO   <- round(cellStats(SM_P_Crop  , 'min'),2)

CF_P <- raster(Mean_CLIM[1])
names(CF_P) = "Cloud Freq"
CF_P_Mask <- mask(x = CF_P, mask = UNIT_X[[u]])
CF_P_Crop <- crop(x = CF_P_Mask, y = extent(UNIT_X[[u]]))
CF_P_M   <- round(cellStats(CF_P_Crop, 'mean'),2)
CF_PName <- names(CF_P)
CFUP   <- round(cellStats(CF_P_Crop  , 'max'),2)
CFLO   <- round(cellStats(CF_P_Crop  , 'min'),2)

Tmax_P <- raster(Mean_CLIM[79])
crs(Tmax_P)  <- crs(EXAMP)
names(Tmax_P) = "Max Air Temp."
Tmax_P_Mask <- mask(x = Tmax_P, mask = UNIT_X[[u]])
Tmax_P_Crop <- crop(x = Tmax_P_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") { Tmax_P_Crop <- (Tmax_P_Crop  * 1.8) + 32}
Tmax_P_M   <- round(cellStats(Tmax_P_Crop, 'mean'),1)
Tmax_PName <- names(Tmax_P)

Tmin_P <- raster(Mean_CLIM[92])
crs(Tmin_P)  <- crs(EXAMP)
names(Tmin_P) = "Mean Air Temp."
Tmin_P_Mask <- mask(x = Tmin_P, mask = UNIT_X[[u]])
Tmin_P_Crop <- crop(x = Tmin_P_Mask, y = extent(UNIT_X[[u]]))
if(TUnit == "°F") { Tmin_P_Crop <- (Tmin_P_Crop  * 1.8) + 32}
Tmin_P_M   <- round(cellStats(Tmin_P_Crop, 'mean'),1)
Tmin_PName <- names(Tmin_P)

KD_P <- raster(Mean_CLIM[53])
names(KD_P) = "Solar Radiation."
KD_P_Mask <- mask(x = KD_P, mask = UNIT_X[[u]])
KD_P_Crop <- crop(x = KD_P_Mask, y = extent(UNIT_X[[u]]))
KD_P_M   <- round(cellStats(KD_P_Crop, 'mean'),0)
KD_PName <- names(KD_P)
KDUP   <- round(cellStats(KD_P_Crop  , 'max'),0)
KDLO   <- round(cellStats(KD_P_Crop  , 'min'),0)

ET_P <- raster(Mean_CLIM[14])
names(ET_P) = "Evaporation"
ET_P_Mask <- mask(x = ET_P, mask = UNIT_X[[u]])
ET_P_Crop <- crop(x = ET_P_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") { ET_P_Crop <- ET_P_Crop * 0.0393701}
ET_P_M   <- round(cellStats(ET_P_Crop , 'mean'),0)
ET_PName <- names(ET_P)
ETUP   <- round(cellStats(ET_P_Crop  , 'max'),0)
ETLO   <- round(cellStats(ET_P_Crop  , 'min'),0)

VPD_P <- raster(Mean_CLIM[105])
names(VPD_P) = "VPD"
VPD_P_Mask <- mask(x = VPD_P, mask = UNIT_X[[u]])
VPD_P_Crop <- crop(x = VPD_P_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") { VPD_P_Crop <- VPD_P_Crop * 0.0393701}
VPD_P_M   <- round(cellStats(VPD_P_Crop, 'mean'),0)
VPD_PName <- names(VPD_P)
VPDUP   <- round(cellStats(VPD_P_Crop  , 'max'),0)
VPDLO   <- round(cellStats(VPD_P_Crop  , 'min'),0)

##########   MEAN RAINFALL Data

Jan <- raster(MeanRF_ALL[1])
Jan_Mask <- mask(x = Jan, mask = UNIT_X[[u]])
Jan_CropRF <- crop(x = Jan_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Jan_CropRF  <- Jan_CropRF  * 0.0393701}
JanMRF    <- round(cellStats(Jan_CropRF,  'mean'),1)
JanMRFx   <- round(cellStats(Jan_CropRF  , 'max'),1)
JanMRFn   <- round(cellStats(Jan_CropRF  , 'min'),1)

Feb <- raster(MeanRF_ALL[2])
Feb_Mask <- mask(x = Feb, mask = UNIT_X[[u]])
Feb_CropRF <- crop(x = Feb_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Feb_CropRF  <- Feb_CropRF  * 0.0393701}
FebMRF   <- round(cellStats(Feb_CropRF,  'mean'),1)
FebMRFx   <- round(cellStats(Feb_CropRF  , 'max'),1)
FebMRFn   <- round(cellStats(Feb_CropRF  , 'min'),1)

Mar <- raster(MeanRF_ALL[3])
Mar_Mask <- mask(x = Mar, mask = UNIT_X[[u]])
Mar_CropRF <- crop(x = Mar_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Mar_CropRF  <- Mar_CropRF  * 0.0393701}
MarMRF   <- round(cellStats(Mar_CropRF,  'mean'),1)
MarMRFx   <- round(cellStats(Mar_CropRF  , 'max'),1)
MarMRFn   <- round(cellStats(Mar_CropRF  , 'min'),1)

Apr <- raster(MeanRF_ALL[4])
Apr_Mask <- mask(x = Apr, mask = UNIT_X[[u]])
Apr_CropRF <- crop(x = Apr_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Apr_CropRF  <- Apr_CropRF  * 0.0393701}
AprMRF   <- round(cellStats(Apr_CropRF,  'mean'),1)
AprMRFx   <- round(cellStats(Apr_CropRF  , 'max'),1)
AprMRFn   <- round(cellStats(Apr_CropRF  , 'min'),1)

May <- raster(MeanRF_ALL[5])
May_Mask <- mask(x = May, mask = UNIT_X[[u]])
May_CropRF <- crop(x = May_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {May_CropRF  <- May_CropRF  * 0.0393701}
MayMRF   <- round(cellStats(May_CropRF,  'mean'),1)
MayMRFx   <- round(cellStats(May_CropRF  , 'max'),1)
MayMRFn   <- round(cellStats(May_CropRF  , 'min'),1)

Jun <- raster(MeanRF_ALL[6])
Jun_Mask <- mask(x = Jun, mask = UNIT_X[[u]])
Jun_CropRF <- crop(x = Jun_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Jun_CropRF  <- Jun_CropRF  * 0.0393701}
JunMRF   <- round(cellStats(Jun_CropRF,  'mean'),1)
JunMRFx   <- round(cellStats(Jun_CropRF  , 'max'),1)
JunMRFn   <- round(cellStats(Jun_CropRF  , 'min'),1)

Jul <- raster(MeanRF_ALL[7])
Jul_Mask <- mask(x = Jul, mask = UNIT_X[[u]])
Jul_CropRF <- crop(x = Jul_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Jul_CropRF  <- Jul_CropRF  * 0.0393701}
JulMRF   <- round(cellStats(Jul_CropRF,  'mean'),1)
JulMRFx   <- round(cellStats(Jul_CropRF  , 'max'),1)
JulMRFn   <- round(cellStats(Jul_CropRF  , 'min'),1)

Aug <- raster(MeanRF_ALL[8])
Aug_Mask <- mask(x = Aug, mask = UNIT_X[[u]])
Aug_CropRF <- crop(x = Aug_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Aug_CropRF  <- Aug_CropRF  * 0.0393701}
AugMRF   <- round(cellStats(Aug_CropRF,  'mean'),1)
AugMRFx   <- round(cellStats(Aug_CropRF  , 'max'),1)
AugMRFn   <- round(cellStats(Aug_CropRF  , 'min'),1)

Sep <- raster(MeanRF_ALL[9])
Sep_Mask <- mask(x = Sep, mask = UNIT_X[[u]])
Sep_CropRF <- crop(x = Sep_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Sep_CropRF  <- Sep_CropRF  * 0.0393701}
SepMRF   <- round(cellStats(Sep_CropRF, 'mean'),1)
SepMRFx   <- round(cellStats(Sep_CropRF  , 'max'),1)
SepMRFn   <- round(cellStats(Sep_CropRF  , 'min'),1)

Oct <- raster(MeanRF_ALL[10]) 
Oct_Mask <- mask(x = Oct, mask = UNIT_X[[u]])
Oct_CropRF <- crop(x = Oct_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Oct_CropRF  <- Oct_CropRF  * 0.0393701}
OctMRF   <- round(cellStats(Oct_CropRF, 'mean'),1)
OctMRFx   <- round(cellStats(Oct_CropRF  , 'max'),1)
OctMRFn   <- round(cellStats(Oct_CropRF  , 'min'),1)

Nov <- raster(MeanRF_ALL[11]) 
Nov_Mask <- mask(x = Nov , mask = UNIT_X[[u]])
Nov_CropRF <- crop(x = Nov_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Nov_CropRF  <- Nov_CropRF  * 0.0393701}
NovMRF   <- round(cellStats(Nov_CropRF,  'mean'),1)
NovMRFx   <- round(cellStats(Nov_CropRF  , 'max'),1)
NovMRFn   <- round(cellStats(Nov_CropRF  , 'min'),1)

Dec <- raster(MeanRF_ALL[12]) 
Dec_Mask <- mask(x = Dec, mask = UNIT_X[[u]])
Dec_CropRF <- crop(x = Dec_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {Dec_CropRF  <- Dec_CropRF  * 0.0393701}
DecMRF   <- round(cellStats(Dec_CropRF,  'mean'),1)
DecMRFx   <- round(cellStats(Dec_CropRF  , 'max'),1)
DecMRFn   <- round(cellStats(Dec_CropRF  , 'min'),1)

ANN <- raster(MeanRF_ALL[13]) 
ANN_Mask <- mask(x = ANN, mask = UNIT_X[[u]])
ANN_CropRF <- crop(x = ANN_Mask, y = extent(UNIT_X[[u]]))
if(RFUnit == " in") {ANN_CropRF  <- ANN_CropRF  * 0.0393701}
ANNMRFM   <- round(cellStats(ANN_CropRF, 'mean'),0)
ANNUP   <- round(cellStats(ANN_CropRF  , 'max'),0)
ANNLO   <- round(cellStats(ANN_CropRF  , 'min'),0)

##########   Add to Table

MEANRF2 <- c(JanMRF,FebMRF,MarMRF,AprMRF,MayMRF,JunMRF,JulMRF,AugMRF,SepMRF,OctMRF,NovMRF,DecMRF,ANNMRFM)
Cell.CL_Year[1,2:14] <- MEANRF2  

########## Get Threshold for figures

RFUP <- max(JanMRFx,FebMRFx,MarMRFx,AprMRFx,MayMRFx,JunMRFx,JulMRFx,AugMRFx,SepMRFx,OctMRFx,NovMRFx,DecMRFx)
RFLO <- min(JanMRFn,FebMRFn,MarMRFn,AprMRFn,MayMRFn,JunMRFn,JulMRFn,AugMRFn,SepMRFn,OctMRFn,NovMRFn,DecMRFn)

########## Write to Output Folder 

write.csv(Cell.CL_Year,paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Annual Climate.csv"),row.names = F)

##########   CLIMO GRAPHS

#build data frame with temperature and precipitation data
df <- as.data.frame(c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
colnames(df) <- c("month")
df$month <- factor(df$month, levels = month.abb)
PREC <- as.numeric(as.vector(Cell.CL_Year[1,2:13]))
TA <- as.numeric(as.vector(Cell.CL_Year[3,2:13]))

df$PREC <- PREC 
df$TA <- TA

#Set X and Y Limits for Graph 20% higher than the max.

XPR <- max(PREC) + (max(PREC)*0.2)
XTA <- max(TA) +   (max(TA)*0.2)
NTA <- min(TA) - min(TA)*0.2

#Make Cliamo Graph and Save to ouptput folder

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Climograph.png"),width=5*dpi,height=5*dpi,res=dpi) 

par(mar=c(4,4,4,4))
my_bar <- barplot(df$PREC, border=F , names.arg=df$month , 
                  las=2 , 
                  col="darkblue" , 
                  ylim=c(0,XPR) ,
                  ylab = paste0("Rainfall [",RFUnit2,"]"))#,
title(paste0("Monthly Climate: ",UNIT_Ns[u]), line = 0.8,cex.main = 1.5)
# text(my_bar, df$PREC+13 ,df$PREC ,cex=1) 
par(new = TRUE)
plot(df$month,df$TA ,pch=15 , lty = 0, axes = FALSE, xlab = "", ylab = "",col="red",ylim=c(NTA,XTA ))
lines(x = df$month, y = df$TA,lwd=2,col="red")
points(x = df$month, y = df$TA,col="red",pch=16,cex=2)
mtext(paste0("Temperature [",TUnit,"]"), side=4, line=2.5,col="red")
axis(4, ylim=c(NTA,XTA), col="red",col.axis="red",las=1)

#Legend
legend("topleft", legend = c("Rainfall", "Temperature"), 
       col = c("darkblue" , "red") , 
       bty = "n", pch=c(15,20) , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.05, 0.05))


dev.off()


########## Figure for Other Cliamte Variables 
# Decide on a Break for Rainfall 

#Mean Rainfall
# Decide on a Break for Rainfall 
RNGERF <- RFUP-RFLO 
if (RNGERF > 8.99){
  BI_brksRF<-round(seq(RFLO, RFUP, length = 10),0);
  colfuncRF<-colorRampPalette(brewer.pal(10,"YlGnBu"))(50)
}
if (RNGERF < 9 && RNGERF > 7.99 ){
  BI_brksRF<-round(seq(RFLO, RFUP, length = 9),0);
  colfuncRF<-colorRampPalette(brewer.pal(9,"YlGnBu"))(50)
}
if (RNGERF < 8 && RNGERF > 6.99 ){
  BI_brksRF<-round(seq(RFLO, RFUP, length = 8),0);
  colfuncRF<-colorRampPalette(brewer.pal(8,"YlGnBu"))(50)
}
if (RNGERF < 7 && RNGERF > 5.99 ){
  BI_brksRF<-round(seq(RFLO, RFUP, length = 7),0);
  colfuncRF<-colorRampPalette(brewer.pal(7,"YlGnBu"))(50)
}
if (RNGERF < 6 && RNGERF > 4.99 ){
  BI_brksRF<-round(seq(RFLO, RFUP, length = 5),0);
  colfuncRF<-colorRampPalette(brewer.pal(5,"YlGnBu"))(50)
}
if (RNGERF < 5 && RNGERF > 3.99 ){
  BI_brksRF<-round(seq(RFLO, RFUP, length = 5),0);
  colfuncRF<-colorRampPalette(brewer.pal(5,"YlGnBu"))(50)
}
if (RNGERF < 4 && RNGERF > 2.99 ){
  BI_brksRF<-round(seq(0, RFUP, length = 4),0);
  colfuncRF<-colorRampPalette(brewer.pal(4,"YlGnBu"))(50)
}
if (RNGERF < 3 && RNGERF > 1.99 ){
  BI_brksRF<-round(seq(0, RFUP, length = 3),0);
  colfuncRF<-colorRampPalette(brewer.pal(3,"YlGnBu"))(50)
}
if (RNGERF < 2 && RNGERF > 0.99 ){
  BI_brksRF<-round(seq(0, RFUP, length = 3),0);
  colfuncRF<-colorRampPalette(brewer.pal(3,"YlGnBu"))(50)
}
if (RNGERF < 2 && RNGERF > 0 ){
  BI_brksRF<-round(seq(0, 3, length = 3),0);
  colfuncRF<-colorRampPalette(brewer.pal(3,"YlGnBu"))(50)
}


#Mean Rainfall
# Decide on a Break for Rainfall 
RNGERFA <- ANNUP-ANNLO 
if (RNGERFA >= 10){
  BI_brksRFA<-round(seq(ANNLO, ANNUP, length = 11),0);
  colfuncRFA<-colorRampPalette(brewer.pal(9,"YlGnBu"))(50)
}
if (RNGERFA < 10 && RNGERFA > 8.99 ){
BI_brksRFA<-round(seq(ANNLO, ANNUP, length = 10),0);
colfuncRFA<-colorRampPalette(brewer.pal(9,"YlGnBu"))(50)
}
if (RNGERFA < 9 && RNGERFA > 7.99 ){
  BI_brksRFA<-round(seq(ANNLO, ANNUP, length = 9),0);
  colfuncRFA<-colorRampPalette(brewer.pal(9,"YlGnBu"))(50)
}
if (RNGERFA < 8 && RNGERFA > 6.99 ){
  BI_brksRFA<-round(seq(ANNLO, ANNUP, length = 8),0);
  colfuncRFA<-colorRampPalette(brewer.pal(8,"YlGnBu"))(50)
}
if (RNGERFA < 7 && RNGERFA > 5.99 ){
  BI_brksRFA<-round(seq(ANNLO, ANNUP, length = 7),0);
  colfuncRFA<-colorRampPalette(brewer.pal(7,"YlGnBu"))(50)
}
if (RNGERFA < 6 && RNGERFA > 4.99 ){
  BI_brksRFA<-round(seq(0, ANNUP, length = 6),0);
  colfuncRFA<-colorRampPalette(brewer.pal(6,"YlGnBu"))(50)
}
if (RNGERFA < 5 && RNGERFA > 3.99 ){
  BI_brksRFA<-round(seq(0, ANNUP, length = 5),0);
  colfuncRFA<-colorRampPalette(brewer.pal(5,"YlGnBu"))(50)
}
if (RNGERFA < 4 && RNGERFA > 0.99 ){
  BI_brksRFA<-round(seq(0, ANNUP, length = 3),0);
  colfuncRFA<-colorRampPalette(brewer.pal(4,"YlGnBu"))(50)
}

#if (RNGERFA < 3){
#  BI_brksRFA<-round(seq(RFLO, RFUP, length = 3),0);
#  colfuncRFA<-colorRampPalette(brewer.pal(3,"YlGnBu"))(50)
#}

#Use Same Scale for TA 
BI_brksTA<-round(seq(TnLO, TxUP, length = 9),0)
colfuncTA <-colorRampPalette(brewer.pal(9,"YlOrRd"))(50)

#For RH 
RNGERH <- RHUP-RHLO 
if (RNGERH > 8){
BI_brksRH<-round(seq(RHLO, RHUP, length = 9),0)
colfuncRH <-colorRampPalette(brewer.pal(9,"BuPu"))(50)
}
# RH
if (RNGERH < 9 && RNGERH > 7.99 ){
  BI_brksRH<-round(seq(RHLO, RHUP, length = 8),0)
  colfuncRH <-colorRampPalette(brewer.pal(8,"BuPu"))(50)
}
if (RNGERH < 8 && RNGERH > 6.99 ){
  BI_brksRH<-round(seq(RHLO, RHUP, length = 7),0)
  colfuncRH <-colorRampPalette(brewer.pal(7,"BuPu"))(50)
}
if (RNGERH < 5 && RNGERH > 3.99 ){
  BI_brksRH<-round(seq(RHLO, RHUP, length = 4),0)
  colfuncRH <-colorRampPalette(brewer.pal(4,"BuPu"))(50)
}
if (RNGERH < 4 && RNGERH > 2.99 ){
  BI_brksRH<-round(seq(RHLO, RHUP, length = 4),0)
  colfuncRH <-colorRampPalette(brewer.pal(4,"BuPu"))(50)
}

#For SM
BI_brksSM<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
colfuncSM <-colorRampPalette(brewer.pal(9,"YlGn"))(50)
#For KD (Solar)
BI_brksKD<-round(seq(KDLO, KDUP, length = 9),0)
colfuncKD <-colorRampPalette(brewer.pal(9,"OrRd"))(50)
#For CF
BI_brksCF<-c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
colfuncCF <-colorRampPalette(brewer.pal(9,"PuRd"))(50)
#For ET
BI_brksET<-round(seq(ETLO, ETUP, length = 9),0)
colfuncET <-colorRampPalette(brewer.pal(9,"PuBu"))(50)

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Climate.png"),width=5*dpi,height=5*dpi,res=dpi)    

title1=textGrob(paste("Average Annual Climate:", UNIT_Ns[u]), gp=gpar(col="darkred",fontface="bold",fontsize=15))  
grid.arrange(top = title1,
             spplot(ANN_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRFA,
                    main=list(label=paste0("RF (",ANNMRFM ,RFUnit,")"),cex=0.8),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])),
             spplot(Tair_P_Crop, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("Mean TA (",Tair_P_M ,TUnit2,")"),cex=0.8),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])), 
             spplot(CF_P_Crop, col.regions = colfuncCF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksCF,
                    main=list(label=paste0("CF (",CF_P_M ,")"),cex=0.8),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])),
             spplot(RH_P_Crop, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("RH (",RH_P_M ," %)"),cex=0.8),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])),
             spplot(Tmin_P_Crop, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("Min TA (",Tmin_P_M ,TUnit2,")"),cex=0.8),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])),
             spplot(KD_P_Crop, col.regions = colfuncKD, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksKD,
                    main=list(label=paste0("S (",KD_P_M ," W/m2)"),cex=0.8),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])),
             spplot(SM_P_Crop, col.regions = colfuncSM, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksSM,
                    main=list(label=paste0("SM (",SM_P_M ,")"),cex=0.8),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])),
             spplot(Tmax_P_Crop, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("Max TA (",Tmax_P_M ,TUnit2,")"),cex=0.8),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])),
             spplot(ET_P_Crop, col.regions = colfuncET, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksET,
                    main=list(label=paste0("ET (",ET_P_M ,RFUnit,")"),cex=0.8),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])))

dev.off()

##########  Temperature Maps Figure 

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," TA12.png"),width=5*dpi,height=5*dpi,res=dpi)    

title1=textGrob(paste("Monthly Temperature:", UNIT_Ns[u]),gp=gpar(col="darkred",fontface="bold",fontsize=15))  
grid.arrange(top = title1,
             spplot(Jan_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("JAN (",JanMTA ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) , 
             spplot(Feb_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("FEB (",FebMTA ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             spplot(Mar_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("MAR (",MarMTA  ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             spplot(Apr_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("APR (",AprMTA  ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             spplot(May_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("MAY (",MayMTA ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             spplot(Jun_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("JUN (",JunMTA  ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             spplot(Jul_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("JUL (",JulMTA  ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             spplot(Aug_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("AUG (",AugMTA ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             spplot(Sep_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("SEP (",SepMTA  ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             spplot(Oct_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("OCT (",OctMTA  ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             spplot(Nov_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("NOV (",NovMTA  ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             spplot(Dec_CropTA, col.regions = colfuncTA, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksTA,
                    main=list(label=paste0("DEC (",DecMTA ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])))

dev.off()

##########   MEAN ANNUAL Relative Humidity 12-maps

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," RH12.png"),width=5*dpi,height=5*dpi,res=dpi)  

title1=textGrob(paste("Monthly Relative Humidity:", UNIT_Ns[u]),gp=gpar(col="darkred",fontface="bold",fontsize=15))   
grid.arrange(top = title1,
             spplot(Jan_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("JAN (",JanMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Feb_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("FEB (",FebMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Mar_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("MAR (",MarMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Apr_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("APR (",AprMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(May_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("MAY (",MayMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Jun_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("JUN (",JunMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Jul_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("JUL (",JulMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Aug_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("AUG (",AugMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Sep_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("SEP (",SepMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Oct_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("OCT (",OctMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Nov_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("NOV (",NovMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Dec_CropRH, col.regions = colfuncRH, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRH,
                    main=list(label=paste0("DEC (",DecMRH ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])))
dev.off()

dpi=300
##########   MEAN ANNUAL Rainfall 12-maps 
#BI_brksRF<-round(seq(0, RFUP+0.3, length = 4),0) #For dry areas.

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," RF12.png"),width=5*dpi,height=5*dpi,res=dpi)  

title1=textGrob(paste("Monthly Rainfall:", UNIT_Ns[u]),gp=gpar(col="darkred",fontface="bold",fontsize=15))   
grid.arrange(top = title1,
             spplot(Jan_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("JAN (",JanMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) , 
             
             spplot(Feb_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("FEB (",FebMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Mar_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("MAR (",MarMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Apr_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("APR (",AprMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(May_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("MAY (",MayMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Jun_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("JUN (",JunMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Jul_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("JUL (",JulMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Aug_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("AUG (",AugMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Sep_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("SEP(",SepMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Oct_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("OCT (",OctMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Nov_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("NOV(",NovMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Dec_CropRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksRF,
                    main=list(label=paste0("DEC (",DecMRF ,RFUnit,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])))

dev.off()



#########   Seasonal Rainfall Maps 

DrySeasonRF <- (May_CropRF+Jun_CropRF+Jul_CropRF+Aug_CropRF+Sep_CropRF+Oct_CropRF)
WetSeasonRF <- (Nov_CropRF+Dec_CropRF+Jan_CropRF+Feb_CropRF+Mar_CropRF+Apr_CropRF)

DryUP   <- round(cellStats(DrySeasonRF, 'max'),1)
DryLO   <- round(cellStats(DrySeasonRF, 'min'),1)
WetUP   <- round(cellStats(WetSeasonRF, 'max'),1)
WetLO   <- round(cellStats(WetSeasonRF, 'min'),1)

DryUPMO   <- round(DryUP/6,1)
DryLOMO   <- round(DryLO/6,1)
WetUPMO   <- round(WetUP/6,1)
WetLOMO   <- round(WetLO/6,1)

SEAUP <- max(DryUP,WetUP)
SEALO <-min(DryLO,WetLO)

BI_brksSEA <- round(seq(SEALO, SEAUP, length = 10),0)

RNGESEA <- SEAUP-SEALO 
if (RNGESEA > 8.99){
  BI_brksSEA<-round(seq(SEALO, SEAUP, length = 11),1);
  colfuncRF<-colorRampPalette(brewer.pal(9,"YlGnBu"))(50)
}
if (RNGESEA < 9 && RNGESEA > 7.99 ){
  BI_brksSEA<-round(seq(SEALO, SEAUP, length = 9),0);
  colfuncRF<-colorRampPalette(brewer.pal(9,"YlGnBu"))(50)
}
if (RNGESEA < 8 && RNGESEA > 6.99 ){
  BI_brksSEA<-round(seq(SEALO, SEAUP, length = 8),0);
  colfuncRF<-colorRampPalette(brewer.pal(8,"YlGnBu"))(50)
}
if (RNGESEA < 7 && RNGESEA > 5.99 ){
  BI_brksSEA<-round(seq(SEALO, SEAUP, length = 7),0);
  colfuncRF<-colorRampPalette(brewer.pal(7,"YlGnBu"))(50)
}
if (RNGESEA < 6 && RNGESEA > 4.99 ){
  BI_brksSEA<-round(seq(SEALO, SEAUP, length = 6),0);
  colfuncRF<-colorRampPalette(brewer.pal(6,"YlGnBu"))(50)
}
if (RNGESEA < 5 && RNGESEA > 3.99 ){
  BI_brksSEA<-round(seq(SEALO, SEAUP, length = 5),0);
  colfuncRF<-colorRampPalette(brewer.pal(5,"YlGnBu"))(50)
}
if (RNGESEA < 4){
  BI_brksSEA<-round(seq(SEALO, SEAUP, length = 4),0);
  colfuncRF<-colorRampPalette(brewer.pal(4,"YlGnBu"))(50)
}



DSeaMRF   <- round(cellStats(DrySeasonRF, 'mean'),1)
WSeaMRF   <- round(cellStats(WetSeasonRF, 'mean'),1)

DryMEMO   <- as.character(round(DSeaMRF/6,1))
DryUPMO   <- as.character(round(DryUP/6,1))
DryLOMO   <- as.character(round(DryLO/6,1))
WetMEMO   <- as.character(round(WSeaMRF/6,1))
WetUPMO   <- as.character(round(WetUP/6,1))
WetLOMO   <- as.character(round(WetLO/6,1))

Cell.DataCLR[1,12] <- DSeaMRF
Cell.DataCLR[2,12] <- DryUP
Cell.DataCLR[3,12] <- DryLO
Cell.DataCLR[1,13] <- WSeaMRF
Cell.DataCLR[2,13] <- WetUP
Cell.DataCLR[3,13] <- WetLO


png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," SeaRF.png"),width=5*dpi,height=5*dpi,res=dpi)  

title1=textGrob(paste("Seasonal Rainfall:", UNIT_Ns[u]),gp=gpar(col="darkred",fontface="bold",fontsize=15)) 
grid.arrange(top = title1,
             spplot(DrySeasonRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksSEA,
                    main=list(label=paste0("Dry Season (MAY-OCT) ", DSeaMRF  ,RFUnit),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) , 
             
             spplot(WetSeasonRF, col.regions = colfuncRF, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksSEA,
                    main=list(label=paste0("Wet Season (NOV-APR) ", WSeaMRF  ,RFUnit),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])))

dev.off()

########## Add Data To tABLE 

Cell.DataCL[u,5:13] <- c(ANNMRFM ,AnnMTA, AnnMTX,AnnMTN,ANNMRH, SM_P_M,KD_P_M,ET_P_M,CF_P_M)

Cell.DataCLR[1,3:11] <- c(ANNMRFM ,AnnMTA, AnnMTX,AnnMTN,ANNMRH, SM_P_M,KD_P_M,ET_P_M,CF_P_M)
Cell.DataCLR[2,3:11] <- c(ANNUP,AnnMTAx,AnnMTXx,AnnMTNx,ANNMRHx, SMUP, KDUP,ETUP,CFUP )
Cell.DataCLR[3,3:11] <- c(ANNLO,AnnMTAn,AnnMTXn,AnnMTNn,ANNMRHn, SMLO, KDLO,ETLO,CFLO)

Cell.RF_Year[u,3:15] <-  MEANRF2

write.csv( Cell.DataCLR, paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u],"Mean Climate.csv"),row.names = F)
Cell.DataCLR


#######################################################################################################################





##########  Downscaling 

##########  Dynamical and Statistical Downscaling 

print("Downscaling")

PWW_T <- UNIT_X[[u]] #Name Variable UNIT

##########  Create a matrix for each cell

Cell.matrix <- matrix(nrow = 36,ncol = 2)
Cell.Data_DS <-data.frame(Cell.matrix)
colnames(Cell.Data_DS) <- c("Metric","Value")

########## Percent Change DyDs RCP 4.5 & 8.5

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

Dy4.5A_MM <- round(cellStats(Dy4.5A_C ,'mean'),0)
Cell.Data_DS[1,1:2]<-c("Dy4.5A_MM",Dy4.5A_MM)

Dy4.5D_MM <- round(cellStats(Dy4.5D_C  ,'mean'),0)
Cell.Data_DS[2,1:2]<-c("Dy4.5D_MM",Dy4.5D_MM)

Dy4.5W_MM <- round(cellStats(Dy4.5W_C ,'mean'),0)
Cell.Data_DS[3,1:2]<-c("Dy4.5W_MM",Dy4.5W_MM)

Dy8.5A_MM <- round(cellStats(Dy8.5A_C ,'mean'),0)
Cell.Data_DS[4,1:2]<-c("Dy8.5A_MM",Dy8.5A_MM)

Dy8.5D_MM <- round(cellStats(Dy8.5D_C ,'mean'),0)
Cell.Data_DS[5,1:2]<-c("Dy8.5D_MM",Dy8.5D_MM)

Dy8.5W_MM <- round(cellStats(Dy8.5W_C ,'mean'),0)
Cell.Data_DS[6,1:2]<-c("Dy8.5W_MM",Dy8.5W_MM)

Dy4.5A_x<- round(cellStats(Dy4.5A_C ,'max'),0)
Dy4.5D_x<- round(cellStats(Dy4.5D_C ,'max'),0)
Dy4.5W_x<- round(cellStats(Dy4.5W_C ,'max'),0)
Dy8.5A_x<- round(cellStats(Dy8.5A_C ,'max'),0)
Dy8.5D_x<- round(cellStats(Dy8.5D_C ,'max'),0)
Dy8.5W_x<- round(cellStats(Dy8.5W_C ,'max'),0)

Dy4.5A_n<- round(cellStats(Dy4.5A_C ,'min'),0)
Dy4.5D_n<- round(cellStats(Dy4.5D_C ,'min'),0)
Dy4.5W_n<- round(cellStats(Dy4.5W_C ,'min'),0)
Dy8.5A_n<- round(cellStats(Dy8.5A_C ,'min'),0)
Dy8.5D_n<- round(cellStats(Dy8.5D_C ,'min'),0)
Dy8.5W_n<- round(cellStats(Dy8.5W_C ,'min'),0)

DyRF100UP <- max(c(Dy4.5A_x, Dy4.5D_x, Dy4.5W_x, Dy8.5A_x, Dy8.5D_x, Dy8.5W_x))
DyRF100LO <- min(c(Dy4.5A_n, Dy4.5D_n, Dy4.5W_n, Dy8.5A_n, Dy8.5D_n, Dy8.5W_n))

########## DyDs Temperature 2100 

#Percent Change DyDs RCP 4.5 & 8.5
DyDs_P_4.5_T_ChgF  <- raster(DyDs_files[37])
if(TUnit == "\u00B0F") {DyDs_P_4.5_T_ChgF  <- raster(DyDs_files[39])}
Dy4.5A_MT <- mask(x = DyDs_P_4.5_T_ChgF, mask = PWW_T)
Dy4.5A_CT <- crop(x = Dy4.5A_MT, y = extent(PWW_T))

DyDs_P_8.5_T_ChgF  <- raster(DyDs_files[38])
if(TUnit == "\u00B0F") {DyDs_P_4.5_T_ChgF  <- raster(DyDs_files[40])}
Dy8.5A_MT <- mask(x = DyDs_P_8.5_T_ChgF, mask = PWW_T)
Dy8.5A_CT <- crop(x = Dy8.5A_MT, y = extent(PWW_T))

###########   Calculate Mean 

Dy4.5A_MMT <- round(cellStats(Dy4.5A_CT ,'mean'),1)
Cell.Data_DS[7,1:2]<-c("Dy4.5A_MMT",Dy4.5A_MMT)
Dy8.5A_MMT <- round(cellStats(Dy8.5A_CT ,'mean'),1)
Cell.Data_DS[8,1:2]<-c("Dy8.5A_MMT",Dy8.5A_MMT)

Dy4.5A_x <- round(cellStats(Dy4.5A_CT ,'max'),1)
Dy8.5A_x <- round(cellStats(Dy8.5A_CT ,'max'),1)
Dy4.5A_n <- round(cellStats(Dy4.5A_CT ,'min'),1)
Dy8.5A_n <- round(cellStats(Dy8.5A_CT ,'min'),1)

##########   StDs Rainfall 2100 

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

St4.5A_MM<- round(cellStats(St4.5A_C ,'mean'),0)
Cell.Data_DS[9,1:2]<-c("St4.5A_MM",St4.5A_MM)

St4.5D_MM<- round(cellStats(St4.5D_C ,'mean'),0)
Cell.Data_DS[10,1:2]<-c("St4.5D_MM",St4.5D_MM)

St4.5W_MM<- round(cellStats(St4.5W_C ,'mean'),0)
Cell.Data_DS[11,1:2]<-c("St4.5W_MM",St4.5W_MM)

St8.5A_MM<- round(cellStats(St8.5A_C ,'mean'),0)
Cell.Data_DS[12,1:2]<-c("St8.5A_MM",St8.5A_MM)

St8.5D_MM<- round(cellStats(St8.5D_C ,'mean'),0)
Cell.Data_DS[13,1:2]<-c("St8.5D_MM",St8.5D_MM)

St8.5W_MM<- round(cellStats(St8.5W_C ,'mean'),0)
Cell.Data_DS[14,1:2]<-c("St8.5W_MM",St8.5W_MM)

St4.5A_x<- round(cellStats(St4.5A_C ,'max'),0)
St4.5D_x<- round(cellStats(St4.5D_C ,'max'),0)
St4.5W_x<- round(cellStats(St4.5W_C ,'max'),0)
St8.5A_x<- round(cellStats(St8.5A_C ,'max'),0)
St8.5D_x<- round(cellStats(St8.5D_C ,'max'),0)
St8.5W_x<- round(cellStats(St8.5W_C ,'max'),0)

St4.5A_n<- round(cellStats(St4.5A_C ,'min'),0)
St4.5D_n<- round(cellStats(St4.5D_C ,'min'),0)
St4.5W_n<- round(cellStats(St4.5W_C ,'min'),0)
St8.5A_n<- round(cellStats(St8.5A_C ,'min'),0)
St8.5D_n<- round(cellStats(St8.5D_C ,'min'),0)
St8.5W_n<- round(cellStats(St8.5W_C ,'min'),0)

StRF100UP <- max(c(St4.5A_x, St4.5D_x, St4.5W_x, St8.5A_x, St8.5D_x, St8.5W_x))
StRF100LO <- min(c(St4.5A_n, St4.5D_n, St4.5W_n, St8.5A_n, St8.5D_n, St8.5W_n))


##########   StDs Rainfall 2040-2070  

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

St4.5A_MM40 <- round(cellStats(St4.5A_C40 ,'mean'),0)
Cell.Data_DS[15,1:2]<-c("St4.5A_MM40",St4.5A_MM40)
St4.5D_MM40<- round(cellStats(St4.5D_C40 ,'mean'),0)
Cell.Data_DS[16,1:2]<-c("St4.5D_MM40",St4.5D_MM40)
St4.5W_MM40<- round(cellStats(St4.5W_C40 ,'mean'),0)
Cell.Data_DS[17,1:2]<-c("St4.5W_MM40",St4.5W_MM40)
St8.5A_MM40<- round(cellStats(St8.5A_C40 ,'mean'),0)
Cell.Data_DS[18,1:2]<-c("St8.5A_MM40",St8.5A_MM40)
St8.5D_MM40<- round(cellStats(St8.5D_C40 ,'mean'),0)
Cell.Data_DS[19,1:2]<-c("St8.5D_MM40",St8.5D_MM40)
St8.5W_MM40<- round(cellStats(St8.5W_C40 ,'mean'),0)
Cell.Data_DS[20,1:2]<-c("St8.5W_MM40",St8.5W_MM40)

##########   StDs Temperature 2100 ################################

##########   Percent Change StDs RCP 4.5 & 8.5

StDs_P_4.5_T_ChgF  <- raster(DyDs_files[116])
if(TUnit == "\u00B0F") {StDs_P_4.5_T_ChgF  <- raster(DyDs_files[120])}
St4.5A_MT <- mask(x = StDs_P_4.5_T_ChgF, mask = PWW_T)
St4.5A_CT <- crop(x = St4.5A_MT, y = extent(PWW_T))

StDs_P_8.5_T_ChgF  <- raster(DyDs_files[118])
if(TUnit == "\u00B0F") {StDs_P_8.5_T_ChgF  <- raster(DyDs_files[122])}
St8.5A_MT <- mask(x = StDs_P_8.5_T_ChgF , mask = PWW_T)
St8.5A_CT <- crop(x = St8.5A_MT, y = extent(PWW_T))

#########   Calculate Mean 
St4.5A_MMT <- round(cellStats(St4.5A_CT ,'mean'),1)
Cell.Data_DS[21,1:2]<-c("St4.5A_MMT",St4.5A_MMT)
St8.5A_MMT <- round(cellStats(St8.5A_CT ,'mean'),1)
Cell.Data_DS[22,1:2]<-c("St8.5A_MMT",St8.5A_MMT)

St4.5A_x <- round(cellStats(St4.5A_CT ,'max'),1)
St8.5A_x <- round(cellStats(St8.5A_CT ,'max'),1)
St4.5A_n <- round(cellStats(St4.5A_CT ,'min'),1)
St8.5A_n <- round(cellStats(St8.5A_CT ,'min'),1)

DsTAUP <- max(Dy4.5A_x, Dy8.5A_x, Dy4.5A_n, Dy8.5A_n, St4.5A_x, St8.5A_x, St4.5A_n, St8.5A_n)
DsTALO <- min(Dy4.5A_x, Dy8.5A_x, Dy4.5A_n, Dy8.5A_n, St4.5A_x, St8.5A_x, St4.5A_n, St8.5A_n)
DsTaup <- max(abs(DsTAUP ),abs(DsTALO ))
DsTAlo <-  DsTaup*-1

#Percent Change StDs RCP 4.5 & 8.5
StDs_P_4.5_T_ChgF40  <- raster(DyDs_files[115])
if(TUnit == "\u00B0F") { StDs_P_4.5_T_ChgF40  <- raster(DyDs_files[119])}
St4.5A_MT40 <- mask(x = StDs_P_4.5_T_ChgF40 , mask = PWW_T)
St4.5A_CT40 <- crop(x = St4.5A_MT40, y = extent(PWW_T))

StDs_P_8.5_T_ChgF40  <- raster(DyDs_files[117])
if(TUnit == "\u00B0F") {StDs_P_8.5_T_ChgF40  <- raster(DyDs_files[121])}
St8.5A_MT40 <- mask(x = StDs_P_8.5_T_ChgF40 , mask = PWW_T)
St8.5A_CT40 <- crop(x = St8.5A_MT40, y = extent(PWW_T))

#Calculate Mean 
St4.5A_MMT40 <- round(cellStats(St4.5A_CT40 ,'mean'),1)
Cell.Data_DS[23,1:2]<-c("St4.5A_MMT40",St4.5A_MMT40)
St8.5A_MMT40 <- round(cellStats(St8.5A_CT40 ,'mean'),1)
Cell.Data_DS[24,1:2]<-c("St8.5A_MMT40",St8.5A_MMT40)


############################################
###########################################
###########################################




########## LULIN ##########################

#Don't do anything with this data yet 
########## Percent Change Lulin DyDs 8.5

D2_P_8.5_ANN <- raster(D2_files[13])
D28.5A_M <- mask(x = D2_P_8.5_ANN, mask = PWW_T)
D28.5A_C <- crop(x = D28.5A_M, y = extent(PWW_T))

D2_P_8.5_Dry <- raster(D2_files[14])
D2.5D_M <- mask(x = D2_P_8.5_Dry, mask = PWW_T)
D2.5D_C <- crop(x = D2.5D_M, y = extent(PWW_T))

D2_P_8.5_Wet <- raster(D2_files[15])
D28.5W_M <- mask(x = D2_P_8.5_Wet, mask = PWW_T)
D28.5W_C <- crop(x = D28.5W_M, y = extent(PWW_T))

D28.5A_MM <- round(cellStats(Dy8.5A_C ,'mean'),0)
Cell.Data_DS[25,1:2]<-c("D28.5A_MM",D28.5A_MM)

D28.5D_MM <- round(cellStats(Dy8.5D_C ,'mean'),0)
Cell.Data_DS[26,1:2]<-c("D28.5D_MM",D28.5D_MM)

D28.5W_MM <- round(cellStats(Dy8.5W_C ,'mean'),0)
Cell.Data_DS[27,1:2]<-c("D28.5W_MM",D28.5W_MM)

D28.5A_x<- round(cellStats(Dy8.5A_C ,'max'),0)
D28.5D_x<- round(cellStats(Dy8.5D_C ,'max'),0)
D28.5W_x<- round(cellStats(Dy8.5W_C ,'max'),0)

D28.5A_n<- round(cellStats(Dy8.5A_C ,'min'),0)
D28.5D_n<- round(cellStats(Dy8.5D_C ,'min'),0)
D28.5W_n<- round(cellStats(Dy8.5W_C ,'min'),0)

#####################################################

D2RF100UP <- max(c(D28.5A_x, D28.5D_x, D28.5W_x))
D2RF100LO <- min(c(D28.5A_n, D28.5D_n, D28.5W_n))

########## D2Ds Temperature 2100 HERE

#Percent Change DyDs RCP 4.5 & 8.5
#HERE
D2Ds_P_8.5_T_ChgF  <- raster(D2_files[28])
if(TUnit == "\u00B0F") {D2Ds_P_8.5_T_ChgF  <- raster(D2_files[31])}
#crs(D2Ds_P_8.5_T_ChgF)<-crs(EXAMP) 
#D2Ds_P_8.5_T_ChgF <- spTransform(D2Ds_P_8.5_T_ChgF, crs(EXAMP))  
D28.5A_MT <- mask(x = D2Ds_P_8.5_T_ChgF, mask = PWW_T)
D28.5A_CT <- crop(x = D28.5A_MT, y = extent(PWW_T))
plot(D28.5A_MT)
###########   Calculate Mean 

D28.5A_MMT <- round(cellStats(D28.5A_CT ,'mean'),1)
Cell.Data_DS[28,1:2]<-c("D28.5A_MMT",D28.5A_MMT)

D28.5A_x <- round(cellStats(D28.5A_CT ,'max'),1)
D28.5A_n <- round(cellStats(D28.5A_CT ,'min'),1)


StRF100UP <- max(c(St4.5A_x, St4.5D_x, St4.5W_x, St8.5A_x, St8.5D_x, St8.5W_x))
StRF100LO <- min(c(St4.5A_n, St4.5D_n, St4.5W_n, St8.5A_n, St8.5D_n, St8.5W_n))
RFdsup <- max(c(abs(DyRF100UP),abs(DyRF100LO),abs(StRF100UP),abs(StRF100LO),abs(D2RF100UP),abs(D2RF100LO)))
RFdslo <- RFdsup *-1


##########   Write all DS results 

#DY - Dynamical
#ST - Statistical
#A - Annual
#W - Wet season
#D - Dry Season
#8.5 is RCP 8.5
#4.5 is RCP 4.5

write.csv(Cell.Data_DS,paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Downscaling.csv"),row.names = F)

##########   FIGURES

##########   Dynamical and Statistical Downscaling 


############################################################################

#These graphs make use of North arrow and scale bar from RF Maps
colfunc2 <- colorRampPalette(brewer.pal(11,"RdBu"))(100)
BI_brks2<-round(seq(RFdslo, RFdsup , length = 9),0)



png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," DS_RF_8.5.png"),width=6.5*dpi,height=4*dpi,res=dpi)

title1=textGrob("Dynamical and Statistical DS RCP 8.5, 2100", gp=gpar(col="darkred",fontface="bold",fontsize=15))  

grid.arrange(top = title1,
             spplot(Dy8.5A_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("DyDs ANN (",Dy8.5A_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St8.5A_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs ANN (",St8.5A_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Dy8.5D_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("DyDs Dry (",Dy8.5D_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St8.5D_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs Dry (",St8.5D_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) , 
             
             spplot(Dy8.5W_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("DyDs Wet (",Dy8.5W_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St8.5W_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs Wet (",St8.5W_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])))          

dev.off()

########## Dynamical and Statistical Downscaling 

########## Downscaling Compare RF RCP 4.5  

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," DS_RF_2100_4.5.png"),width=6.5*dpi,height=4*dpi,res=dpi)

title1=textGrob("Dynamical and Statistical DS RCP 4.5, 2100", gp=gpar(col="darkred",fontface="bold",fontsize=15))  


grid.arrange(top = title1,
             spplot(Dy4.5A_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("DyDs ANN (",Dy4.5A_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St4.5A_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs ANN (",Dy4.5A_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Dy4.5D_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("DyDs DRY (",Dy4.5D_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St4.5D_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs Dry (",St4.5D_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) , 
             
             spplot(Dy4.5W_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("DyDs WET (",Dy4.5W_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St4.5W_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs WET (",St4.5W_MM ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])))         

dev.off()        

########### Statistical 2040-2070

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," StDsRF2040.png"),width=6.5*dpi,height=4*dpi,res=dpi)

title1=textGrob("Statistical DS RCP 8.5 & RCP 4.5, 2040-2070)", gp=gpar(col="darkred",fontface="bold",fontsize=15))  

grid.arrange(top = title1,
             spplot(St4.5A_C40, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs ANN (",St4.5A_MM40 ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St8.5A_C40, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs ANN (",St8.5A_MM40 ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St4.5D_C40, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs DRY (",St4.5D_MM40 ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St8.5D_C40, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs DRY (",St8.5D_MM40 ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St4.5W_C40, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs WET (",St4.5W_MM40 ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St8.5W_C, col.regions = colfunc2, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brks2,
                    main=list(label=paste0("StDs WET (",St8.5W_MM40 ,"%)"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])))           
dev.off()

##########   Downscaling Compare TEMP RCP 8.5 & 4.5 2100

colfunc3 <- colorRampPalette(brewer.pal(9,"Reds"))(100)

if(TUnit == "°C") {BI_brksF<-round(seq(0, 7 , length = 7),0)}
if(TUnit == "°F") {BI_brksF<-round(seq(0, 9 , length = 9),0)}

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," DS_Temp2100.png"),width=6.5*dpi,height=4*dpi,res=dpi)

title4=textGrob("Dynamical & Statistical DS Compare TEMP (Yr 2100)", gp=gpar(col="darkred",fontface="bold",fontsize=15))
grid.arrange(top = title4,
             spplot(Dy4.5A_CT, col.regions = colfunc3, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksF,
                    main=list(label=paste0("DyDs ANN RCP 4.5 (",Dy4.5A_MMT ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St4.5A_CT, col.regions = colfunc3, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksF,
                    main=list(label=paste0("StDs ANN RCP 4.5 (",St4.5A_MMT ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(Dy8.5A_CT, col.regions = colfunc3, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksF,
                    main=list(label=paste0("DyDs ANN RCP 8.5 (",Dy8.5A_MMT ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) ,
             
             spplot(St8.5A_CT, col.regions = colfunc3, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksF,
                    main=list(label=paste0("StDs ANN RCP 8.5 (",St8.5A_MMT ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])))
dev.off()

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," StDs_Temp2040.png"),width=6.5*dpi,height=4*dpi,res=dpi)  

title4=textGrob("Statistical DS Compare TEMP (Yr 2040-2070)", gp=gpar(col="darkred",fontface="bold",fontsize=15))

grid.arrange(top = title4,
             spplot(St4.5A_CT40, col.regions = colfunc3, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksF,
                    main=list(label=paste0("StDs ANN RCP 4.5 (",St4.5A_MMT40 ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])) , 
             
             spplot(St8.5A_CT40, col.regions = colfunc3, equal=FALSE,
                    axes = TRUE,las = 1, cex.axis=0.7,at=BI_brksF,
                    main=list(label=paste0("StDs ANN RCP 4.5 (",St8.5A_MMT40 ,TUnit2,")"),cex=0.9),
                    colorkey = list(space = "right", height = 1, labels=list(cex=0.6)),
                    sp.layout = list(UNIT_X[u])))

dev.off()

##########################################################################################################################




##########  Create A Monthly Rainfall Time Series

print("Rainfall Extract")  
print("Frazier et al 2016")

#########   Extract RF Data From Monthly Maps (2 DATASETS)
#########   Call in the Shapefile
PWW_T <- UNIT_X[[u]]

#########   Load RF MAPS

RF_Tif_files = dir(RF_Map_Path_A, pattern="*.tif", recursive=T, full.names=T)  #Monthly RF
nfiles <- length(RF_Tif_files)

#Create a matrix for each cell
Cell.AF_Maps <-data.frame(matrix(nrow = 1116,ncol = 4))
colnames(Cell.AF_Maps) <- c("Date","Year","Month","RF")

for (i in 1:nfiles)        {
  if(i == 1){print("AF")}
  if(i == 250){print(i)} #Just to see where Im at in the process 
  if(i == 500){print(i)}
  if(i == 750){print(i)}
  if(i == 1000){print(i)}
  
  #Open Map as raster and change projection
  RF_Map <- RF_Tif_files[i]
  RF_Map2 <- raster(RF_Map)
  #crs(RF_Map2) <- " +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  plot(RF_Map2)
  
  #Get the Name (date) from each map (derek's version)
  name <- basename(RF_Map)
  nameSplit <- unlist(strsplit(name,"_"))
  nameSplit2 <- unlist(strsplit(nameSplit[8] ,".t"))
  Year <- nameSplit[7]
  Month <- nameSplit2[1]
  MY <- paste0(Year,"/",Month)
  
  
  projection(RF_Map2) <- projection(EXAMP)
  RFM_Mask <- mask(x = RF_Map2, mask = PWW_T)
  RFM_Crop <- crop(x = RFM_Mask, y = extent(PWW_T))
  RFM_M   <- round(cellStats(RFM_Crop, 'mean'),1)
  
  Cell.AF_Maps[i,1:4]<-c(MY,Year,Month,RFM_M)
  
}

head(Cell.AF_Maps)
tail(Cell.AF_Maps)

##########   Matty's Maps

print("Lucas")


#Load Daily RF MAPS
RF_Tif_files = dir(RF_Map_Path, pattern="*.tif", recursive=T, full.names=T)  #Monthly RF
nfiles <- length(RF_Tif_files)

#Create a matrix for each cell
Cell.ML_Maps <-data.frame(matrix(nrow = 360,ncol = 4))
colnames(Cell.ML_Maps) <- c("Date","Year","Month","RF")

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
  nameSplit <- unlist(strsplit(name,"_"))
  nameSplit2 <- unlist(strsplit(nameSplit[6] ,".t"))
  D1 <- nameSplit2[1]
  Year <- nameSplit[1]
  Month <- nameSplit[2]
  MY <- paste0(Year,"/",Month)
  
  projection(RF_Map2) <- projection(EXAMP)
  RFM_Mask <- mask(x = RF_Map2, mask = PWW_T)
  RFM_Crop <- crop(x = RFM_Mask, y = extent(PWW_T))
  RFM_M   <- round(cellStats(RFM_Crop, 'mean'),1)

  ##########  ADD Data To Matrix
  
  Cell.ML_Maps[i,1:4]<-c(MY,Year,Month,RFM_M)
  
} #End Month RF Extract

head(Cell.ML_Maps)
tail(Cell.ML_Maps)

##########  Do a comparision with the 23-year overalap period and generate a figure

##########  Compare Monthly Rainfall fromthe two map products          

########## Get Common Period 1990-2012 (23 years)

###

###

# ### If working with exported RF dataframes ###
#     Cell.AF_Maps2<-read.csv("Haleakala National Park/MINI_Phase2 outputs 1/Cell.AF_Maps_9_6_2022.csv")
#     Cell.AF_Maps2<-data.frame(Cell.AF_Maps2[2:ncol(Cell.AF_Maps2)])
#     head(Cell.AF_Maps2)
# 
#     # summary(Cell.AF_Maps2[which(Cell.AF_Maps2$Year<2020),]$RF)
# 
#     Cell.ML_Maps2<-read.csv("Haleakala National Park/MINI_Phase2 outputs 1/Cell.ML_Maps_9_6_2022.csv")
#     Cell.ML_Maps2<-Cell.ML_Maps2[2:ncol(Cell.ML_Maps2)]
#     head(Cell.ML_Maps2)

###
    
###
    
###

### Abby's data
  # subset for 1990 - 2012
  MRF_A2 =  Cell.AF_Maps[c(841:1116),]
  head(MRF_A2)
  
  # just keep RF values
  MRF_A3 =  as.numeric( MRF_A2[,4])
  
  # convert from mm to inches
  if(RFUnit == " in") {MRF_A3 <-  MRF_A3 * 0.0393701} 
  summary(MRF_A3)

### Matty's data
  # subet for 1990 - 2012
MRF_N2 =  Cell.ML_Maps[c(1:276),]
head(MRF_N2)

# just keep RF values
MRF_N3 =  as.numeric(MRF_N2[,4])

# convert mm to inches
if(RFUnit == " in") {MRF_N3  <-  MRF_N3  * 0.0393701} 
summary(MRF_N3)



MBE <- round(mean(MRF_A3 - MRF_N3),1)
MAE <- round(mean(abs(MRF_A3 - MRF_N3)),1)

D_Comp <- cbind(MRF_A3,MRF_N3)
Mx <- max(D_Comp)


FNAME <- paste0("RF_Compare_23_",UNIT_Ns[u],".csv")

##########   Comparision Figures For the two datasets

LM1 <- lm(MRF_A3~MRF_N3)
LM1P <- round(coefficients(summary(LM1))[2,4],4)
LM1R <- round(summary(LM1)$r.squared,2)
TITLE = paste(UNIT_Ns[u], " 23-yr RF Compare (",RFUnit2,")")

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u], " 23yr_RF_Compare.png"),width=5*dpi,height=5*dpi,res=dpi)

plot(MRF_A3~MRF_N3,ylim=c(0,Mx),xlim=c(0,Mx),main=TITLE,
     ylab = "Frazier et al. (2016)",xlab="Lucas et al. (In Review)")
abline(0,1)
legend("topleft",c(paste("R2 = ",LM1R),paste("MBE = ",MBE),paste("MAE = ",MAE)))

dev.off()

##########   Make the same figure in a different folder.
# Derek added "SFOLDER" (same location as RFOLDER)
SFOLDER<-paste0("E:/PDKE/CCVD/MINI_Phase2/", UNIT_N[u],"/RF_Compare/")

png(paste0(SFOLDER, UNIT_N[u]," 23yr_RF_Compare.png"),width=5*dpi,height=5*dpi,res=dpi)

plot(MRF_A3~MRF_N3,ylim=c(0,Mx),xlim=c(0,Mx),main=TITLE,
     ylab = "ABBY",xlab="NEW")
abline(0,1)
legend("topleft",c(paste("R2 = ",LM1R),paste("MBE = ",MBE),paste("MAE = ",MAE)))

dev.off()

##########   Create a 100-Year timeseries and do analyais

##########   Merge Datasets to create common 100-year period 

MRF_AD3 =  Cell.AF_Maps[c(1:840),]
head(MRF_AD3)
tail(MRF_AD3)
MRF_ND3 =  Cell.ML_Maps[c(1:360),]
head(MRF_ND3)

colnames(MRF_AD3) <- c("Date","Year","Month","RF")
colnames(MRF_ND3) <- c("Date","Year","Month","RF")

MRF100 <- rbind(MRF_AD3,MRF_ND3)
if(RFUnit == " in") {MRF100$RF <-  as.numeric(MRF100$RF) * 0.0393701} 

write.csv(MRF100,paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Monthly Rainfall_",RFUnit2,".csv"),row.names = F)

# ### Can read in the csv from above to start script here
#       MRF100<-read.csv("Haleakala National Park/Haleakala National Park Monthly Rainfall_in.csv")
#       # MRF100<-MRF100[2:ncol(MRF100)]


head(MRF100)
tail(MRF100)

RF_IN <- as.numeric(MRF100[,4]) 
summary(RF_IN)

MRF100$Day <- 1
MRF100$Ndate <- as.Date(with( MRF100, paste(Year, Month, Day,sep="-")), "%Y-%m-%d")

MRF_Max <- round(max(RF_IN,na.rm=T),1)
MRF_Min <- round(min(RF_IN,na.rm=TRUE),1)
MRF_MED <- round(median(RF_IN,na.rm=TRUE),1)
MRF_MEAN <- round(mean(RF_IN,na.rm=TRUE),1)

MoRF.ts <- ts(RF_IN, c(1920,1), end = c(2019,12), frequency = 12) 

myts1 <- as.vector(window(MoRF.ts, start=c(1920, 1), end=c(2019, 12)))
myts2 <- as.vector(window(MoRF.ts, start=c(1940, 1), end=c(2019, 12)))
myts3 <- as.vector(window(MoRF.ts, start=c(1960, 1), end=c(2019, 12)))
myts4 <- as.vector(window(MoRF.ts, start=c(1980, 1), end=c(2019, 12)))
myts5 <- as.vector(window(MoRF.ts, start=c(2000, 1), end=c(2019, 12)))
myts6 <- as.vector(window(MoRF.ts, start=c(2010, 1), end=c(2019, 12)))

DateT1 <- seq(as.Date("1920-01-01"), as.Date("2019-12-31"), by="months")
DateT2 <- seq(as.Date("1940-01-01"), as.Date("2019-12-31"), by="months")
DateT3 <- seq(as.Date("1960-01-01"), as.Date("2019-12-31"), by="months")
DateT4 <- seq(as.Date("1980-01-01"), as.Date("2019-12-31"), by="months")
DateT5 <- seq(as.Date("2000-01-01"), as.Date("2019-12-31"), by="months")
DateT6 <- seq(as.Date("2010-01-01"), as.Date("2019-12-31"), by="months")

LM1 <- lm(myts1~DateT1)
LM2 <- lm(myts2~DateT2)
LM3 <- lm(myts3~DateT3)
LM4 <- lm(myts4~DateT4)
LM5 <- lm(myts5~DateT5)
LM6 <- lm(myts6~DateT6)

T1 <- round(coefficients(summary(LM1))[2,1],2)
T2 <- round(coefficients(summary(LM2))[2,1],2)
T3 <- round(coefficients(summary(LM3))[2,1],2)
T4 <- round(coefficients(summary(LM4))[2,1],2)
T5 <- round(coefficients(summary(LM5))[2,1],2)
T6 <- round(coefficients(summary(LM6))[2,1],2)

LM1P <- coefficients(summary(LM1))[2,4]
LM2P <- coefficients(summary(LM2))[2,4]
LM3P <- coefficients(summary(LM3))[2,4]
LM4P <- coefficients(summary(LM4))[2,4]
LM5P <- coefficients(summary(LM5))[2,4]
LM6P <- coefficients(summary(LM6))[2,4]

LM1R <- summary(LM1)$r.squared
LM2R <- summary(LM2)$r.squared
LM3R <- summary(LM3)$r.squared
LM4R <- summary(LM4)$r.squared
LM5R <- summary(LM5)$r.squared
LM6R <- summary(LM6)$r.squared

########## 2003-2019   

MA02 <- round(mean(myts6),1)
ME02 <- round(median(myts6),1)
MX02 <- round(max(myts6),1)
MI02 <- round(min(myts6),1)

##########   Aggregate Month to Year 

##########   Aggregate From Daily to Daily Average  

short.date_M = strftime(MRF100$Ndate, "%Y/%m")

Mean_M_RF = aggregate(as.numeric(MRF100$RF) ~ short.date_M, FUN = mean)
colnames(Mean_M_RF) <- c("Date","RF")
short.date_Y = strftime(as.Date(MRF100$Ndate), "%Y")

Mean_Y_RF = aggregate(as.numeric(MRF100$RF) ~ short.date_Y, FUN = mean)
colnames(Mean_Y_RF) <- c("Date","RF")

##########   Plot Annual RF  

YrRF.ts <- ts(Mean_Y_RF$RF, c(1920), end = c(2019), frequency = 1) 
myts1Y <- as.vector(window(YrRF.ts, start=c(1920), end=c(2019)))
myts2Y <- as.vector(window(YrRF.ts, start=c(1940), end=c(2019)))
myts3Y <- as.vector(window(YrRF.ts, start=c(1960), end=c(2019)))
myts4Y <- as.vector(window(YrRF.ts, start=c(1980), end=c(2019)))
myts5Y <- as.vector(window(YrRF.ts, start=c(2000), end=c(2019)))
myts6Y <- as.vector(window(YrRF.ts, start=c(2010), end=c(2019)))

YDateT1 <- seq(as.Date("1920-01-01"), as.Date("2019-12-31"), by="years")
YDateT2 <- seq(as.Date("1940-01-01"), as.Date("2019-12-31"), by="years")
YDateT3 <- seq(as.Date("1960-01-01"), as.Date("2019-12-31"), by="years")
YDateT4 <- seq(as.Date("1980-01-01"), as.Date("2019-12-31"), by="years")
YDateT5 <- seq(as.Date("2000-01-01"), as.Date("2019-12-31"), by="years")
YDateT6 <- seq(as.Date("2010-01-01"), as.Date("2019-12-31"), by="years")

LM1Y <- lm(myts1Y~YDateT1)
LM2Y <- lm(myts2Y~YDateT2)
LM3Y <- lm(myts3Y~YDateT3)
LM4Y <- lm(myts4Y~YDateT4)
LM5Y <- lm(myts5Y~YDateT5)
LM6Y <- lm(myts6Y~YDateT6)

T1Y <- round(coefficients(summary(LM1Y))[2,1],3)
T1Y<- round(coefficients(summary(LM1Y))[2,1], 3)
T2Y <- round(coefficients(summary(LM2Y))[2,1],3)
T3Y <- round(coefficients(summary(LM3Y))[2,1],3)
T4Y <- round(coefficients(summary(LM4Y))[2,1],3)
T5Y <- round(coefficients(summary(LM5Y))[2,1],3)
T6Y <- round(coefficients(summary(LM6Y))[2,1],3)

LM1PY <- round(coefficients(summary(LM1Y))[2,4],2)
LM2PY <- round(coefficients(summary(LM2Y))[2,4],2)
LM3PY <- round(coefficients(summary(LM3Y))[2,4],2)
LM4PY <- round(coefficients(summary(LM4Y))[2,4],2)
LM5PY <- round(coefficients(summary(LM5Y))[2,4],2)
LM6PY <- round(coefficients(summary(LM6Y))[2,4],2)

LM1RY <- round(summary(LM1Y)$r.squared,2)
LM2RY <- round(summary(LM2Y)$r.squared,2)
LM3RY <- round(summary(LM3Y)$r.squared,2)
LM4RY <- round(summary(LM4Y)$r.squared,2)
LM5RY <- round(summary(LM5Y)$r.squared,2)
LM6RY <- round(summary(LM6Y)$r.squared,2)

##########   Seasonal RF 

# Wet Season
WET_RF <- subset(MRF100, Month==c('01','02','03','04','11','12'))
head(WET_RF)
WET_RF2 <- rbind(c(NA,"12",NA,NA,NA), WET_RF)
WET_RF3 <- rbind(c(NA,"11",NA,NA,NA), WET_RF2)
WET_RF4 <- as.numeric(WET_RF3$RF)

# WET_RF<-MRF100[MRF100$Month !="5" & MRF100$Month !="6" & MRF100$Month != "7" & 
#                  MRF100$Month != "8" & MRF100$Month != "9" & MRF100$Month != "10",]
# WET_RF5<-as.numeric(WET_RF$RF)

#Get seasonal average 
WET_RF5 <-  as.vector(tapply(WET_RF4, gl(length(WET_RF4)/6, 6), mean,na.rm=T))
DRY_RF <-MRF100[MRF100$Month != "01" & MRF100$Month  != "02" & 
                  MRF100$Month  != "03" & MRF100$Month  != "04" & 
                  MRF100$Month  != 11 & MRF100$Month  != 12,]  

########## Season Stats

DRY_RF2 <- as.numeric(DRY_RF$RF)

##########   Wet
W_MRF_Max <- round(max(WET_RF5,na.rm=T),0)
W_MRF_Min <- round(min(WET_RF5,na.rm=TRUE),0)
W_MRF_MED <- round(median(WET_RF5,na.rm=TRUE),0)
W_MRF_MEAN <- round(mean(WET_RF5,na.rm=TRUE),0)
##########   Dry
D_MRF_Max <- round(max(DRY_RF2,na.rm=T),0)
D_MRF_Min <- round(min(DRY_RF2,na.rm=TRUE),0)
D_MRF_MED <- round(median(DRY_RF2,na.rm=TRUE),0)
D_MRF_MEAN <- round(mean(DRY_RF2,na.rm=TRUE),0)




short.date_Y = strftime(DRY_RF$Ndate, "%Y")

Dry_RF = aggregate(as.numeric(DRY_RF$RF) ~ short.date_Y, FUN = mean)
colnames(Dry_RF) <- c("Date","RF")

####### Seasonal Trends ########################

#WET Season 
YrRF.tsW <- ts(WET_RF5, c(1920), end = c(2019), frequency = 1) 
myts1YW <- as.vector(window(YrRF.tsW, start=c(1920), end=c(2019)))
myts2YW <- as.vector(window(YrRF.tsW, start=c(1940), end=c(2019)))
myts3YW <- as.vector(window(YrRF.tsW, start=c(1960), end=c(2019)))
myts4YW <- as.vector(window(YrRF.tsW, start=c(1980), end=c(2019)))
myts5YW <- as.vector(window(YrRF.tsW, start=c(2000), end=c(2019)))
myts6YW <- as.vector(window(YrRF.tsW, start=c(2010), end=c(2019)))

#DRy Season 
YrRF.tsD <- ts(DRY_RF2, c(1920), end = c(2019), frequency = 1) 
myts1YD <- as.vector(window(YrRF.tsD, start=c(1920), end=c(2019)))
myts2YD <- as.vector(window(YrRF.tsD, start=c(1940), end=c(2019)))
myts3YD<- as.vector(window(YrRF.tsD, start=c(1960), end=c(2019)))
myts4YD <- as.vector(window(YrRF.tsD, start=c(1980), end=c(2019)))
myts5YD <- as.vector(window(YrRF.tsD, start=c(2000), end=c(2019)))
myts6YD <- as.vector(window(YrRF.tsD, start=c(2010), end=c(2019)))

# WET and DRy # Annual 
YDateT1 <- seq(as.Date("1920-01-01"), as.Date("2019-12-31"), by="years")
YDateT2 <- seq(as.Date("1940-01-01"), as.Date("2019-12-31"), by="years")
YDateT3 <- seq(as.Date("1960-01-01"), as.Date("2019-12-31"), by="years")
YDateT4 <- seq(as.Date("1980-01-01"), as.Date("2019-12-31"), by="years")
YDateT5 <- seq(as.Date("2000-01-01"), as.Date("2019-12-31"), by="years")
YDateT6 <- seq(as.Date("2010-01-01"), as.Date("2019-12-31"), by="years")

#WET Regression
LM1YW <- lm(myts1YW~YDateT1)
LM2YW <- lm(myts2YW~YDateT2)
LM3YW <- lm(myts3YW~YDateT3)
LM4YW <- lm(myts4YW~YDateT4)
LM5YW <- lm(myts5YW~YDateT5)
LM6YW <- lm(myts6YW~YDateT6)

# Dry Regression 
LM1YD <- lm(myts1YD~YDateT1)
LM2YD <- lm(myts2YD~YDateT2)
LM3YD <- lm(myts3YD~YDateT3)
LM4YD <- lm(myts4YD~YDateT4)
LM5YD <- lm(myts5YD~YDateT5)
LM6YD <- lm(myts6YD~YDateT6)

#WET SLOPE
T1YW <- round(coefficients(summary(LM1YW))[2,1],3)
T2YW <- round(coefficients(summary(LM2YW))[2,1],3)
T3YW <- round(coefficients(summary(LM3YW))[2,1],3)
T4YW <- round(coefficients(summary(LM4YW))[2,1],3)
T5YW <- round(coefficients(summary(LM5YW))[2,1],3)
T6YW <- round(coefficients(summary(LM6YW))[2,1],3)

#DRY SLOPE
T1YD <- round(coefficients(summary(LM1YD))[2,1],3)
T2YD <- round(coefficients(summary(LM2YD))[2,1],3)
T3YD <- round(coefficients(summary(LM3YD))[2,1],3)
T4YD <- round(coefficients(summary(LM4YD))[2,1],3)
T5YD <- round(coefficients(summary(LM5YD))[2,1],3)
T6YD <- round(coefficients(summary(LM6YD))[2,1],3)

#Wet Pvale
LM1PYW <- round(coefficients(summary(LM1YW))[2,4],2)
LM2PYW <- round(coefficients(summary(LM2YW))[2,4],2)
LM3PYW <- round(coefficients(summary(LM3YW))[2,4],2)
LM4PYW <- round(coefficients(summary(LM4YW))[2,4],2)
LM5PYW <- round(coefficients(summary(LM5YW))[2,4],2)
LM6PYW <- round(coefficients(summary(LM6YW))[2,4],2)

#Dry Pvale
LM1PYD <- round(coefficients(summary(LM1YD))[2,4],2)
LM2PYD <- round(coefficients(summary(LM2YD))[2,4],2)
LM3PYD <- round(coefficients(summary(LM3YD))[2,4],2)
LM4PYD <- round(coefficients(summary(LM4YD))[2,4],2)
LM5PYD <- round(coefficients(summary(LM5YD))[2,4],2)
LM6PYD <- round(coefficients(summary(LM6YD))[2,4],2)

#Wet Regression 
LM1RYW <- round(summary(LM1YW)$r.squared,2)
LM2RYW <- round(summary(LM2YW)$r.squared,2)
LM3RYW <- round(summary(LM3YW)$r.squared,2)
LM4RYW <- round(summary(LM4YW)$r.squared,2)
LM5RYW <- round(summary(LM5YW)$r.squared,2)
LM6RYW <- round(summary(LM6YW)$r.squared,2)

#Dry Regression 
LM1RYD <- round(summary(LM1YD)$r.squared,2)
LM2RYD <- round(summary(LM2YD)$r.squared,2)
LM3RYD <- round(summary(LM3YD)$r.squared,2)
LM4RYD <- round(summary(LM4YD)$r.squared,2)
LM5RYD <- round(summary(LM5YD)$r.squared,2)
LM6RYD <- round(summary(LM6YD)$r.squared,2)


dpi<-300
png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," RF_Trend.png"),width=6.3*dpi,height=7*dpi,res=dpi)

##########   Annual and Seasonal Plot 

par(mfrow=c(3,1))
par(mar=c(4,4,4,2))

MLIM1 <- max(c(myts1Y,myts1YW,myts1YD)) 
YLIM <-  min(myts1Y) 
MLIM <-  (MLIM1 + (MLIM1 *0.45))

par(mai=c(0.3,0.6,0.2,0.2))
plot(myts1Y~YDateT1,ylab = paste0("Average Rainfall (",RFUnit2,"/month)"),type="l",col="blue",xlab="",xaxt="n",ylim=c(YLIM,MLIM),cex.axis =1.3,las=1)
title(paste("Rainfall Trend 1920-2019:",UNIT_Ns[u]), line = 0.5,cex.main = 1.5)

legend("topright", c(paste0("1920-2019 R2 = ",LM1RY, " p = ",LM1PY),
                     #paste0("1940-2019 Trend =", T2Y,", R2 =",LM2RY, " p = ",LM2PY),
                     paste0("1960-2019 R2 = ",LM3RY, " p = ",LM3PY),
                     #paste0("1980-2019 Trend =", T4Y,", R2 =",LM4RY, " p = ",LM4PY),
                     paste0("2000-2019 R2 = ",LM5RY, " p = ",LM5PY)),
                     #paste0("2010-2019 Trend =", T6Y,", R2 =",LM6RY, " p = ",LM6PY)),
       #lty = 1, col = c("darkred","darkorange","darkgreen","darkblue","purple","darkcyan"), lwd = 3)
       lty = 1, col = c("grey70","grey30","grey1"), lwd = 3)
       
legend("topleft",c("Annual"),cex=1.5,text.font=2, bty = "n")

ablineclip(lm(myts1Y~YDateT1),x1=-19000,x2=18000,col=alpha("grey70",0.9),lwd=3)
#ablineclip(lm(myts2Y~YDateT2),x1=-11000,x2=18000,col="darkorange",lwd=3)
ablineclip(lm(myts3Y~YDateT3),x1=-4500,x2=18000,col=alpha("grey30",0.7),lwd=3)
#ablineclip(lm(myts4Y~YDateT4),x1=3500,x2=18000,col="darkblue",lwd=3)
ablineclip(lm(myts5Y~YDateT5),x1=11000,x2=18000,col=alpha("grey1",0.7),lwd=3)
#ablineclip(lm(myts6Y~YDateT6),x1=14000,x2=18000,col="darkcyan",lwd=3)

####### Wet Season ###########
par(mai=c(0.3,0.6,0.2,0.2))
YLIM <-  min(myts1YW) 
plot(myts1YW~YDateT1,ylab = paste0("Average Rainfall (",RFUnit2,"/month)"),type="l",col="blue",xlab="",xaxt="n",ylim=c(YLIM,MLIM),cex.axis =1.3,las=1)

legend("topright", c(paste0("1920-2019 R2 = ",LM1RYW, " p = ",LM1PYW),
                     #paste0("1940-2019 Trend =", T2YW,", R2 =",LM2RYW, " p = ",LM2PYW),
                     paste0("1960-2019 R2 = ",LM3RYW, " p = ",LM3PYW),
                     #paste0("1980-2019 Trend =", T4YW,", R2 =",LM4RYW, " p = ",LM4PYW),
                     paste0("2000-2019 R2 = ",LM5RYW, " p = ",LM5PYW)),
                     #paste0("2010-2019 Trend =", T6YW,", R2 =",LM6RYW, " p = ",LM6PYW)),
       #lty = 1, col = c("darkred","darkorange","darkgreen","darkblue","purple","darkcyan"), lwd = 3,cex=1)
       lty = 1, col = c("grey70","grey30","grey1"), lwd = 3)

legend("topleft",c("Wet Season"),cex=1.5,text.font=2, bty = "n")

ablineclip(lm(myts1YW~YDateT1),x1=-19000,x2=18000,col="grey70",lwd=3)
#ablineclip(lm(myts2YW~YDateT2),x1=-11000,x2=18000,col="darkorange",lwd=3)
ablineclip(lm(myts3YW~YDateT3),x1=-4500,x2=18000,col="grey30",lwd=3)
#ablineclip(lm(myts4YW~YDateT4),x1=3500,x2=18000,col="darkblue",lwd=3)
ablineclip(lm(myts5YW~YDateT5),x1=11000,x2=18000,col="grey1",lwd=3)
#ablineclip(lm(myts6YW~YDateT6),x1=14000,x2=18000,col="darkcyan",lwd=3)

########## Dry Season 

par(mai=c(0.3,0.6,0.2,0.2))
YLIM <-  min(myts1YD) 
plot(myts1YD~YDateT1,ylab = paste0("Average Rainfall (",RFUnit2,"/month)"),type="l",col="blue",xlab="",ylim=c(YLIM,MLIM),cex.axis =1.3,las=1)
# axis(1, labels = T)
# title(main = "Average Wet Season Rainfall Pu'u Wa'awa'a (1920-2019)", line = 1)
legend("topright", c(paste0("1920-2019 R2 = ",LM1RYD, " p = ",LM1PYD),
                     #paste0("1940-2019 Trend =", T2YD,", R2 =",LM2RYD, " p = ",LM2PYD),
                     paste0("1960-2019 R2 = ",LM3RYD, " p = ",LM3PYD),
                     #paste0("1980-2019 Trend =", T4YD,", R2 =",LM4RYD, " p = ",LM4PYD),
                     paste0("2000-2019 R2 = ",LM5RYD, " p = ",LM5PYD)),
                     #paste0("2010-2019 Trend =", T6YD,", R2 =",LM6RYD, " p = ",LM6PYD)),
       #lty = 1, col = c("darkred","darkorange","darkgreen","darkblue","purple","darkcyan"), lwd = 3)
       lty = 1, col = c("grey70","grey30","grey1"), lwd = 3)
       
legend("topleft",c("Dry Season"),cex=1.5,text.font=2, bty = "n")

ablineclip(lm(myts1YD~YDateT1),x1=-19000,x2=18000,col="grey70",lwd=3)
#ablineclip(lm(myts2YD~YDateT2),x1=-11000,x2=18000,col="darkorange",lwd=3)
ablineclip(lm(myts3YD~YDateT3),x1=-4500,x2=18000,col="grey30",lwd=3)
#ablineclip(lm(myts4YD~YDateT4),x1=3500,x2=18000,col="darkblue",lwd=3)
ablineclip(lm(myts5YD~YDateT5),x1=11000,x2=18000,col="grey1",lwd=3)
#ablineclip(lm(myts6YD~YDateT6),x1=14000,x2=18000,col="darkcyan",lwd=3)

dev.off() 

############################################################################# 
# Count All Drought Events Long-term 12 and Short term 12 and 3. 
Cell.SPICNT <-data.frame(matrix(ncol = 4,nrow=4))
colnames(Cell.SPICNT) <- c("Total", "SPI_12_Long","SPI_3_Short", "SPI_12_Short")
Cell.SPICNT[1:4,1] <- c("Drought Events", "Moderate","Severe", "Extreme")

### Ryan's code
print("SPI FIG 12")
RF <- as.numeric(MRF100$RF)
library(SPEI)
RFDATA <- cbind(MRF100[,c(2,3)])

#SPI-12 Calculation
SPI12 <- spi(RF, scale = 12, distribution = 'Gamma')
#SPI-3 Calculation
SPI3 <- spi(RF, scale = 3, distribution = 'Gamma')

#PlOT ALL SPI
png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," SPI.png"),width=6.5*dpi,height=4*dpi,res=dpi)

plot(spi(ts(RF,freq=12,start=c(1920,1)),scale = 12,distribution = 'Gamma'),main = paste0("SPI-12 1920-2019: ",UNIT_Ns[u]))

dev.off()
# 
# print("SPI-12-TAB") 
# SPI_12M <- SPI12
# SPI_3M <- SPI3
# 
# SPIVEC <- as.vector(SPI_12M$fitted)
# SPIVEC3 <- as.vector(SPI_3M$fitted)
# 
# SPIVEC[is.na(SPIVEC )] <- 0
# SPIVEC3[is.na(SPIVEC3 )] <- 0
# 
# DATE5 <-   MRF100$Ndate
# SPI_Data <- SPIVEC
# SPI_Data3 <- SPIVEC3 
# SPID <- cbind(DATE5,SPI_Data,SPI_Data3)
# SPID[1:20]
# 
# Nrow_SPI <- sum(!is.na(SPI_Data))
# SM<-0  
# PH2<-0 
# SMVec<-vector()

##########   Table Metrics SPI 12

Cell.DataSPI <-data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI) <- c("Start","End", "Duration","A Intensity", "P Intensity", "Magnitude")
Cell.DataSPI

##### Derek's Drought Code from Guam
    # load rainfall dataset created around line 2828 above
    wd<-paste0(RFOLDER,UNIT_N,"/")
    setwd(wd)
    csv<-paste0(UNIT_N[u]," Monthly Rainfall_",RFUnit2,".csv")
    RF_Data <- read.csv(csv)
    head(RF_Data)
    
    # calculate SPI and SPEI
    library(SPEI)
    
    spi<-spi(ts(RF_Data$RF, freq=12, start=c(1920, 1)), scale = 12)
    plot(spi)
    spi3<-spi(ts(RF_Data$RF, freq=12, start=c(1920,1)), scale = 3)
    plot(spi3)
    
    # write spi to dataframe
    spi_val<-spi$fitted
    spi_val<-data.frame(spi=as.matrix(spi_val), date=time(spi_val))
    spi_val$m.scale<-12
    spi_val$date<-sub("\\..*","",spi_val$date)
    
    head(spi_val, 20)
    
    spi3_val<-spi3$fitted
    spi3_val<-data.frame(spi=as.matrix(spi3_val), date=time(spi3_val))
    spi3_val$m.scale<-3
    spi3_val$date<-sub("\\..*","",spi3_val$date)
    
    head(spi3_val, 20)
    
    # combine dataframes
    SPI_ALL<-rbind(spi_val,spi3_val)
    head(SPI_ALL, 20)
    
    # fix date column
    library(zoo)
    SPI_ALL$date2<-rep(c(1:12))
    SPI_ALL$date<-as.Date(paste0(SPI_ALL$date,"/",SPI_ALL$date2, "/01"), format = "%Y/%m/%d")
    
    SPI_ALL<-SPI_ALL[1:3]
    colnames(SPI_ALL)<-c("SPI","date","m.scale")
    
    # sort by date
    SPI_ALL<-SPI_ALL[order(as.Date(SPI_ALL$date)),]
    
    # Convert drought events (negative values) to positive
    SPI_ALL$spi_negs<-ifelse(SPI_ALL$SPI>0,0,SPI_ALL$SPI)
    SPI_ALL$spi_negs<-abs(SPI_ALL$spi_negs)
    
    summary(SPI_ALL$spi_negs)
    head(SPI_ALL, 50)

    ### Drought Event Count for SPI-12 ###
    
        spi2<-SPI_ALL[which(SPI_ALL$m.scale == 12),]
      
        # create binary column of drought (SPI>0) yes or no
        spi2$drought<-ifelse(spi2$spi_negs>0, 1, 0)
        head(spi2, 50)
        
        # make column of consecutive drought months
        library(data.table)
        
        spi2$DRGHT_ct<-with(spi2, (drought == 1)*
                              ave(drought, rleid(drought == 1),
                                  FUN=seq_along))
        
        summary(spi2$DRGHT_ct)
        
        # subset data for drought event months only
        spi3<-spi2[which(spi2$DRGHT_ct >= 1),]
        head(spi3, 20)
        
        # create empty event count column
        spi3$event_ct<-0
      
        ### fill in event_ct column
        for (r in 1:nrow(spi3)) {
          
          SPI_M<-spi3[r,]
          SPI_M
          
          if(r == 1) {spi3[r,]$event_ct<-r}
          if(SPI_M$DRGHT_ct > 1) {spi3[r,]$event_ct<-spi3[r-1,]$event_ct}
          if(SPI_M$DRGHT_ct == 1 && r > 1) {spi3[r,]$event_ct<-spi3[r-1,]$event_ct + 1}
        }
        
        head(spi3,100)
        summary(spi3$event_ct)
        
        ### keep only events with peak spi_negs >=1
        event<-unique(spi3$event_ct)
        
        for (x in event) {
          
          sub<-spi3[which(spi3$event_ct == x),]
          
          if(max(sub$spi_negs)<1) {spi3[which(spi3$event_ct == x),]$event_ct<-0}
        }
        
        head(spi3, 100)
        
        ### make final event_ct column
          # subset for events only
          spi4<-spi3[which(spi3$event_ct>0),]
          spi4$event_ct2<-NA
          spi4
          
          event<-data.frame(unique(spi4$event_ct))
          event
          
          # fill event_Ct2 column with event number
          for (x in 1:nrow(event)){
            
            # get event number
            e<-event[x,]
            
            # for each event, add event_ct2 value
            spi4[which(spi4$event_ct == e),]$event_ct2<-x
          }
        
        spi4
        
        # merge event_ct back into full SPI dataset
        spi5<-merge(spi2, spi4, all=T)
        
        head(spi5, 10)
        summary(spi5$event_ct2)
        
        # get rid of wrong event column
        spi5<-subset(spi5, select=-c(event_ct))
        
        colnames(spi5)[which(names(spi5) == "event_ct2")]<-"event_ct"
        
        # sort by date
        spi5<-spi5[order(as.Date(spi5$date)),]
        head(spi5,20)

        ### create column which labels each drought event by intensity
        # create new copy of dataframe
        spi6<-spi5
        
        # get total events count as list
        events<-as.list(1:max(spi6$event_ct, na.rm=T))
        
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
        
        head(spi6)
        head(spi6,100)
        summary(spi6$intensity)
        summary(spi6$peak)
        summary(spi6$months)
        
        ### make table of start and end dates for each drought intensity event

        Cell.DataSPI <-data.frame(matrix(ncol = 6))
        colnames(Cell.DataSPI) <- c("Start","End", "Duration","A Intensity", "P Intensity", "Magnitude")
        Cell.DataSPI$Start<-as.Date(Cell.DataSPI$Start)
        Cell.DataSPI$End<-as.Date(Cell.DataSPI$End)
        Cell.DataSPI
        colnames(Cell.DataSPI)

        ### loop through rows and fill table
        
        for (y in 1:nrow(spi6)) {

          # get row
          row<-spi6[y,]
          row

          # get previous and next row
          prev<-spi6[y-1,]

          nex<-spi6[y+1,]

          # get event number
          number<-row$event_ct

          if(!is.na(row$intensity) && is.na(prev$intensity)) {Cell.DataSPI[number,]$Start<-as.Date(nex$date)}
          #if(!is.na(row$intensity) && is.na(prev$intensity)) {Cell.DataSPI[number,]$intensity<-row$intensity}
          if(!is.na(row$intensity) && is.na(nex$intensity)) {Cell.DataSPI[number,]$End<-as.Date(nex$date)}
          if(!is.na(row$intensity)) {Cell.DataSPI[number,]$Duration<-row$months}
          #if(!is.na(row$intensity)) {Cell.DataSPI[number,]$event_ct<-row$event_ct}
          if(!is.na(row$intensity)) {Cell.DataSPI[number,]$`P Intensity`<-row$peak}
          if(!is.na(row$intensity)) {Cell.DataSPI[number,]$`A Intensity`<-row$mean}
          if(!is.na(row$intensity)) {Cell.DataSPI[number,]$Magnitude<-row$mag}
        }
        
        # fix date formats
        # library(zoo)
        # Cell.DataSPI$start<-as.yearmon(sub(" .*","", Cell.DataSPI$start), "%Y-%m-%d")
        # Cell.DataSPI$stop<-as.yearmon(sub(" .*","", Cell.DataSPI$stop), "%Y-%m-%d")
        
Cell.DataSPI  

###### Back to Ryan's code #####

write.csv(Cell.DataSPI,paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Drought History.csv"),row.names = F)

##########   Count Droughts For Figure - Derek's edits
Cell.DataSPI$`P Intensity`<-as.numeric(Cell.DataSPI$`P Intensity`)
    
EX_Cnt <- sum(Cell.DataSPI[5] > 2)
SV_Cnt <- sum(Cell.DataSPI[5] > 1.5)
MO_Cnt <- sum(Cell.DataSPI[5] > 1)

D_Cnt <- MO_Cnt
SV_Cnt2 <- SV_Cnt - EX_Cnt
MO_Cnt2 <- D_Cnt -  SV_Cnt2 - EX_Cnt

# create SPI count dataframe
Cell.SPICNT <-data.frame(matrix(ncol = 4,nrow=4))
colnames(Cell.SPICNT) <- c("Total", "SPI_12_Long","SPI_3_Short", "SPI_12_Short")
Cell.SPICNT[1:4,1] <- c("Drought Events", "Moderate","Severe", "Extreme")

Cell.SPICNT[1:4,2] <- c(D_Cnt,MO_Cnt2,SV_Cnt2,EX_Cnt)
Cell.SPICNT

    # ##########   Count Droughts For Figure - Ryan's version
    # 
    # EX_Cnt <- sum(Cell.DataSPI[5] > -1.99)
    # SV_Cnt <- sum(Cell.DataSPI[5] > -1.49)
    # MO_Cnt <- sum(Cell.DataSPI[5] > -0.99)
    # D_Cnt <- MO_Cnt
    # SV_Cnt2 <- SV_Cnt - EX_Cnt
    # MO_Cnt2 <- D_Cnt -  SV_Cnt2 - EX_Cnt
    # 
    # Cell.SPICNT[1:4,2] <- c(D_Cnt,MO_Cnt2,SV_Cnt2,EX_Cnt)

##########   Remove Postive SPI and make absloute values 

SPIVEC[SPIVEC  > 0 ] <- 0
head(SPI_ALL)

### Derek's code
SPIVEC<-SPI_ALL[which(SPI_ALL$m.scale == 12),]$SPI
SPIVEC[SPIVEC > 0] <- 0
SPIVEC_Abs <- as.vector(abs(SPIVEC))

M_SPI.ts <- ts(SPIVEC_Abs, c(1920,1), end = c(2019,12), frequency = 12)
myts66 <- as.vector(window(M_SPI.ts, start=c(1920, 1), end=c(2019, 12)))
DateT1 <- as.Date(seq(as.Date("1920-01-01"), as.Date("2019-12-31"), by="months"))
#short.date_M = strftime(MonthlyRF$Ndate, "%Y-%m")

xx <- data.frame(DateT1,myts66)
colnames(xx)<- c("DT","SP")
xx$DT <- as.Date(xx$DT)

write.csv(xx,          paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u],"SPI_12.csv"),row.names = F)

dpi=300
png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u],"Drought_History.png"),width=6.5*dpi,height=4*dpi,res=dpi)

print(ggplot(xx, aes(x = DT, y = SP)) +
        geom_area(fill="darkorange", color="black") +
        xlab("") +
        
        labs(title = paste0("SPI-12 Drought Events 1920 -2019: ",UNIT_Ns[u]),
             x = "",
             y = "Drought Intensity") +
        geom_hline(yintercept=2, linetype="dashed", color = "darkred", size = 1) + 
        geom_hline(yintercept=1.5, linetype="dashed", color = "red", size = 1) +
        geom_hline(yintercept=1, linetype="dashed", color = "orange", size = 1)  +
        annotate("text", x = xx$DT[350], y = 3.4, label = paste0("Drought Events = ",  D_Cnt),size=5, fontface="bold") +
        annotate("text", x = xx$DT[350], y = 3.1, label = paste0("Moderate Droughts = ", MO_Cnt2),size=5, fontface="bold",colour = "orange") +
        annotate("text", x = xx$DT[350], y = 2.8, label = paste0("Severe Droughts = ", SV_Cnt2),size=5, fontface="bold",colour = "red") +
        annotate("text", x = xx$DT[350], y = 2.5, label = paste0("Extreme Droughts = ", EX_Cnt),size=5, fontface="bold",colour = "darkred"))


dev.off()





##########################################################################################################################################

#Short (_S) timescales - Derek's version

#### SPI-3

head(SPI_ALL)
spi2 <- SPI_ALL[which(SPI_ALL$m.scale == 3),]
spi2 <- spi2[which(spi2$date >= as.Date("1990-01-01")),]
head(spi2)

# create binary column of drought (SPI>=1) yes or no
spi2$drought<-ifelse(spi2$spi_negs>0, 1, 0)
head(spi2, 30)

# make column of consecutive drought months
library(data.table)

spi2$DRGHT_ct<-with(spi2, (drought == 1)*
                      ave(drought, rleid(drought == 1),
                          FUN=seq_along))

# subset data for drought event months only
spi3<-spi2[which(spi2$drought >= 1),]
head(spi3, 20)

# create empty event count column
spi3$event_ct<-0

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
event<-unique(spi3$event_ct)

for (x in event) {
  
  sub<-spi3[which(spi3$event_ct == x),]
  
  if(max(sub$spi_negs)<1) {spi3[which(spi3$event_ct == x),]$event_ct<-0}
}

head(spi3, 100)

### make final event_ct column
# subset for events only
spi4<-spi3[which(spi3$event_ct>0),]
spi4$event_ct2<-NA
spi4

event<-data.frame(unique(spi4$event_ct))
event

# fill event_Ct2 column with event number
for (x in 1:nrow(event)){
  
  # get event number
  e<-event[x,]
  
  # for each event, add event_ct2 value
  spi4[which(spi4$event_ct == e),]$event_ct2<-x
}

spi4

# merge event_ct back into full SPI dataset
spi5<-merge(spi2, spi4, all=T)

# sort by date
spi5<-spi5[order(as.Date(spi5$date)),]
head(spi5, 10)
summary(spi5$event_ct2)

# get rid of wrong event column
spi5<-subset(spi5, select=-c(event_ct))
colnames(spi5)[which(names(spi5) == "event_ct2")]<-"event_ct"

### create column which labels each drought event by intensity
# create new copy of dataframe
spi6<-spi5

# get total events count as list
events<-as.list(1:max(spi6$event_ct, na.rm=T))

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

head(spi6,100)
summary(spi6$intensity)
summary(spi6$peak)
summary(spi6$months)

### make table of start and end dates for each drought intensity event

Cell.DataSPI3_S <-data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI3_S) <- c("Start","End", "Duration","A Intensity", "P Intensity", "Magnitude")

Cell.DataSPI3_S$Start<-as.POSIXct(Cell.DataSPI3_S$Start)
Cell.DataSPI3_S$End<-as.POSIXct(Cell.DataSPI3_S$End)
Cell.DataSPI3_S
colnames(Cell.DataSPI3_S)

### loop through rows and fill table

for (y in 1:nrow(spi6)) {
  
  # get row
  row<-spi6[y,]
  row
  
  # get previous and next row
  prev<-spi6[y-1,]
  
  nex<-spi6[y+1,]
  
  # get event number
  number<-row$event_ct
  
  if(!is.na(row$intensity) && is.na(prev$intensity)) {Cell.DataSPI3_S[number,]$Start<-as.Date(nex$date)}
  # if(!is.na(row$intensity) && is.na(prev$intensity)) {event.t[number,]$intensity<-row$intensity}
  if(!is.na(row$intensity) && is.na(nex$intensity)) {Cell.DataSPI3_S[number,]$End<-as.Date(nex$date)}
  if(!is.na(row$intensity)) {Cell.DataSPI3_S[number,]$Duration<-row$months}
  # if(!is.na(row$intensity)) {event.t[number,]$event_ct<-row$event_ct}
  if(!is.na(row$intensity)) {Cell.DataSPI3_S[number,]$`P Intensity`<-row$peak}
  if(!is.na(row$intensity)) {Cell.DataSPI3_S[number,]$`A Intensity`<-row$mean}
  if(!is.na(row$intensity)) {Cell.DataSPI3_S[number,]$Magnitude<-row$mag}
}

# fix date formats
library(zoo)
Cell.DataSPI3_S$Start<-as.yearmon(sub(" .*","", Cell.DataSPI3_S$Start), "%Y-%m-%d")
Cell.DataSPI3_S$End<-as.yearmon(sub(" .*","", Cell.DataSPI3_S$End), "%Y-%m-%d")

Cell.DataSPI3_S

##########   Count Droughts For Figure - Derek's edits
Cell.DataSPI3_S$`P Intensity`<-as.numeric(Cell.DataSPI3_S$`P Intensity`)

EX_Cnt <- sum(Cell.DataSPI3_S[5] > 2)
SV_Cnt <- sum(Cell.DataSPI3_S[5] > 1.5)
MO_Cnt <- sum(Cell.DataSPI3_S[5] > 1)

D_Cnt <- MO_Cnt
SV_Cnt2 <- SV_Cnt - EX_Cnt
MO_Cnt2 <- D_Cnt -  SV_Cnt2 - EX_Cnt

Cell.SPICNT[1:4,3] <- c(D_Cnt,MO_Cnt2,SV_Cnt2,EX_Cnt)
Cell.SPICNT

##########   Remove Postive SPI and make absloute values 
### Derek's Version ###

# use full SPI (3 and 12) dataframe
head(SPI_ALL)

# subset for SPI3 and 1990 - 2020 time period
SPIVEC<-SPI_ALL[which(SPI_ALL$m.scale == 3),]
SPIVEC<-SPIVEC[which(SPIVEC$date >= as.Date("1990-01-01")),]$SPI

# change positive SPI values to 0
SPIVEC[SPIVEC > 0] <- 0

# convert negatives to absolute
SPIVEC_Abs <- as.vector(abs(SPIVEC))

# make dataframe of date and drought (absolute SPI) value
M_SPI.ts <- ts(SPIVEC_Abs, c(1990,1), end = c(2019,12), frequency = 12)
myts66 <- as.vector(window(M_SPI.ts, start=c(1990, 1), end=c(2019, 12)))
DateT1 <- as.Date(seq(as.Date("1990-01-01"), as.Date("2019-12-31"), by="months"))
xx <- data.frame(DateT1,myts66)
colnames(xx)<- c("DT","SP")
xx$DT <- as.Date(xx$DT)
xx

# plot and write figure
dpi=300
png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u],"Drought_HistoryS_3.png"),width=6.5*dpi,height=4*dpi,res=dpi)

print(ggplot(xx, aes(x = DT, y = SP)) +
        geom_area(fill="darkorange", color="black") +
        xlab("") +
        
        labs(title = paste0("SPI-3 Drought Events 1990-2019: ",UNIT_Ns[u]),
             x = "",
             y = "Drought Intensity") +
        geom_hline(yintercept=2, linetype="dashed", color = "darkred", size = 1) + 
        geom_hline(yintercept=1.5, linetype="dashed", color = "red", size = 1) +
        geom_hline(yintercept=1, linetype="dashed", color = "orange", size = 1)  +
        annotate("text", x = xx$DT[70], y = 3.4, label = paste0("Drought Events = ",  D_Cnt),size=5, fontface="bold") +
        annotate("text", x = xx$DT[70], y = 3.1, label = paste0("Moderate Droughts = ", MO_Cnt2),size=5, fontface="bold",colour = "orange") +
        annotate("text", x = xx$DT[70], y = 2.8, label = paste0("Severe Droughts = ", SV_Cnt2),size=5, fontface="bold",colour = "red") +
        annotate("text", x = xx$DT[70], y = 2.5, label = paste0("Extreme Droughts = ", EX_Cnt),size=5, fontface="bold",colour = "darkred")+
        annotate("text", x = xx$DT[300], y = 3.1, label = "SPI-3",size=12, fontface="bold",colour = "black"))

dev.off()


######################################################################################################################################    




#### SPI-12 short

head(SPI_ALL)
spi2 <- SPI_ALL[which(SPI_ALL$m.scale == 12),]
spi2 <- spi2[which(spi2$date >= as.Date("1990-01-01")),]
head(spi2)

# create binary column of drought (SPI>=1) yes or no
spi2$drought<-ifelse(spi2$spi_negs>0, 1, 0)
head(spi2, 30)

# make column of consecutive drought months
library(data.table)

spi2$DRGHT_ct<-with(spi2, (drought == 1)*
                      ave(drought, rleid(drought == 1),
                          FUN=seq_along))

# subset data for drought event months only
spi3<-spi2[which(spi2$DRGHT_ct >= 1),]
head(spi3, 20)

# create empty event count column
spi3$event_ct<-0

### fill in event_ct column
for (r in 1:nrow(spi3)) {
  
  SPI_M<-spi3[r,]
  SPI_M
  
  if(r == 1) {spi3[r,]$event_ct<-r}
  if(SPI_M$DRGHT_ct > 1) {spi3[r,]$event_ct<-spi3[r-1,]$event_ct}
  if(SPI_M$DRGHT_ct == 1 && r > 1) {spi3[r,]$event_ct<-spi3[r-1,]$event_ct + 1}
}

head(spi3,100)
summary(spi3$event_ct)

### keep only events with peak spi_negs >=1
event<-unique(spi3$event_ct)

for (x in event) {
  
  sub<-spi3[which(spi3$event_ct == x),]
  
  if(max(sub$spi_negs)<1) {spi3[which(spi3$event_ct == x),]$event_ct<-0}
}

head(spi3, 100)

### make final event_ct column
# subset for events only
spi4<-spi3[which(spi3$event_ct>0),]
spi4$event_ct2<-NA
spi4

event<-data.frame(unique(spi4$event_ct))
event

# fill event_Ct2 column with event number
for (x in 1:nrow(event)){
  
  # get event number
  e<-event[x,]
  
  # for each event, add event_ct2 value
  spi4[which(spi4$event_ct == e),]$event_ct2<-x
}

spi4

# merge event_ct back into full SPI dataset
spi5<-merge(spi2, spi4, all=T)

# sort by date
spi5<-spi5[order(as.Date(spi5$date)),]
head(spi5, 10)
summary(spi5$event_ct2)

# get rid of wrong event column
spi5<-subset(spi5, select=-c(event_ct))
colnames(spi5)[which(names(spi5) == "event_ct2")]<-"event_ct"

### create column which labels each drought event by intensity
# create new copy of dataframe
spi6<-spi5

# get total events count as list
events<-as.list(1:max(spi6$event_ct, na.rm=T))

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

head(spi6,100)
summary(spi6$intensity)
summary(spi6$peak)
summary(spi6$months)

### make table of start and end dates for each drought intensity event

Cell.DataSPI12_S <-data.frame(matrix(ncol = 6))
colnames(Cell.DataSPI12_S) <- c("Start","End", "Duration","A Intensity", "P Intensity", "Magnitude")

Cell.DataSPI12_S$Start<-as.POSIXct(Cell.DataSPI12_S$Start)
Cell.DataSPI12_S$End<-as.POSIXct(Cell.DataSPI12_S$End)
Cell.DataSPI12_S
colnames(Cell.DataSPI12_S)

### loop through rows and fill table

for (y in 1:nrow(spi6)) {
  
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
  if(!is.na(row$intensity)) {Cell.DataSPI12_S[number,]$Duration<-row$months}
  # if(!is.na(row$intensity)) {event.t[number,]$event_ct<-row$event_ct}
  if(!is.na(row$intensity)) {Cell.DataSPI12_S[number,]$`P Intensity`<-row$peak}
  if(!is.na(row$intensity)) {Cell.DataSPI12_S[number,]$`A Intensity`<-row$mean}
  if(!is.na(row$intensity)) {Cell.DataSPI12_S[number,]$Magnitude<-row$mag}
}

# fix date formats
library(zoo)
Cell.DataSPI12_S$Start<-as.yearmon(sub(" .*","", Cell.DataSPI12_S$Start), "%Y-%m-%d")
Cell.DataSPI12_S$End<-as.yearmon(sub(" .*","", Cell.DataSPI12_S$End), "%Y-%m-%d")

Cell.DataSPI12_S
Cell.DataSPI12_S[which(Cell.DataSPI12_S$`P Intensity` > 1.5),]


##########   Count Droughts For Figure - Derek's edits
Cell.DataSPI12_S$`P Intensity`<-as.numeric(Cell.DataSPI12_S$`P Intensity`)

EX_Cnt <- sum(Cell.DataSPI12_S[5] > 2)
SV_Cnt <- sum(Cell.DataSPI12_S[5] > 1.5)
MO_Cnt <- sum(Cell.DataSPI12_S[5] > 1)

D_Cnt <- MO_Cnt
SV_Cnt2 <- SV_Cnt - EX_Cnt
MO_Cnt2 <- D_Cnt -  SV_Cnt2 - EX_Cnt

Cell.SPICNT[1:4,4] <- c(D_Cnt,MO_Cnt2,SV_Cnt2,EX_Cnt)
Cell.SPICNT
Cell.SPICNT

write.csv(Cell.SPICNT,paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," Drought Count.csv"),row.names = F)

##########   Remove Postive SPI and make absloute values 
### Derek's Version ###

# use full SPI (3 and 12) dataframe
head(SPI_ALL)

# subset for SPI12 and 1990 - 2020 time period
SPIVEC<-SPI_ALL[which(SPI_ALL$m.scale == 12),]
SPIVEC<-SPIVEC[which(SPIVEC$date >= as.Date("1990-01-01")),]$SPI

# change positive SPI values to 0
SPIVEC[SPIVEC > 0] <- 0

# convert negatives to absolute
SPIVEC_Abs <- as.vector(abs(SPIVEC))

# make dataframe of date and drought (absolute SPI) value
M_SPI.ts <- ts(SPIVEC_Abs, c(1990,1), end = c(2019,12), frequency = 12)
myts66 <- as.vector(window(M_SPI.ts, start=c(1990, 1), end=c(2019, 12)))
DateT1 <- as.Date(seq(as.Date("1990-01-01"), as.Date("2019-12-31"), by="months"))
xx <- data.frame(DateT1,myts66)
colnames(xx)<- c("DT","SP")
xx$DT <- as.Date(xx$DT)
xx

# plot and write figure
dpi=300
png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u],"Drought_HistoryS_12.png"),width=6.5*dpi,height=4*dpi,res=dpi)

print(ggplot(xx, aes(x = DT, y = SP)) +
        geom_area(fill="darkorange", color="black") +
        xlab("") +
        
        labs(title = paste0("SPI-12 Drought Events 1990-2019: ",UNIT_Ns[u]),
             x = "",
             y = "Drought Intensity") +
        geom_hline(yintercept=2, linetype="dashed", color = "darkred", size = 1) + 
        geom_hline(yintercept=1.5, linetype="dashed", color = "red", size = 1) +
        geom_hline(yintercept=1, linetype="dashed", color = "orange", size = 1)  +
        annotate("text", x = xx$DT[70], y = 3.4, label = paste0("Drought Events = ",  D_Cnt),size=5, fontface="bold") +
        annotate("text", x = xx$DT[70], y = 3.1, label = paste0("Moderate Droughts = ", MO_Cnt2),size=5, fontface="bold",colour = "orange") +
        annotate("text", x = xx$DT[70], y = 2.8, label = paste0("Severe Droughts = ", SV_Cnt2),size=5, fontface="bold",colour = "red") +
        annotate("text", x = xx$DT[70], y = 2.5, label = paste0("Extreme Droughts = ", EX_Cnt),size=5, fontface="bold",colour = "darkred")+
        annotate("text", x = xx$DT[300], y = 3.1, label = "SPI-12",size=12, fontface="bold",colour = "black"))

dev.off()


################################################################################################################################################     

##########  MEI 

##########   UNIT

print("MEI")

Cell.MEI <- data.frame(matrix(nrow = 5,ncol = 7))
colnames(Cell.MEI) <- c("Phase","W-Mean","D-Mean","W-Max","D-Max","W-Min","D-Min")
Cell.MEI[1:5,1] <- c("Strong EL","Weak EL", "Neutral","Weak LA","Strong LA")

##########   All MEI 
#Cell.MEI_All[u:1] <- UNIT_N[u] 
#Cell.MEI_All[u,8] <- W_MRF_MEAN
#Cell.MEI_All[UNIT_C+u,8] <- D_MRF_MEAN

MEI_W <- MEI$MEI_W
MEI_D <- MEI$MEI_D

MRF  =  Cell.AF_Maps[c(1:1116),]
head(MRF)

##########  Need to remove Rows to get seasons correct

MRF2 =  MRF[-c(1115:1116),]
MRF3 =  MRF2[-c(1:4),]
head(MRF3)
tail(MRF3)
#Aggregate every 6 months
colnames(MRF3)[4] <- "RANCH_RF"
RF6 <- as.numeric(MRF3$RANCH_RF)
SixMo <-colMeans(matrix(RF6, nrow=6))  #Will Need to change This 

if(RFUnit == " in") {SixMo <-  SixMo * 0.0393701} 

SixMo[seq(1,length(SixMo),2)]

DryRF <- SixMo[c(TRUE, FALSE)]
WetRF <- SixMo[c(FALSE, TRUE)]

##########   Bind MEI and Data

L0_W <-cbind(WetRF,MEI_W)
L0_D <-cbind(DryRF,MEI_D)

##########   Seperate by ENSO Phase 

ELN_W  <- subset(L0_W, MEI_W > 0.5)
LAN_W  <- subset(L0_W, MEI_W < -0.5)
NUT_Wx <- subset(L0_W, MEI_W > -0.5)
MEI_W2 <- NUT_Wx[,2]
NUT_W  <- subset(NUT_Wx,MEI_W2 < 0.5)

##########   Strong and Weak EL Wet Season 
W0_X_EL <- ELN_W[,2]
ELN_W_Weak <- subset(ELN_W , W0_X_EL  < 1)
ELN_W_Strong <- subset(ELN_W , W0_X_EL  > 1)

##########   Strong and Weak La Wet Season 
W0_X_LA <- LAN_W[,2]
LAN_W_Weak <- subset(LAN_W , W0_X_LA  > -1)
LAN_W_Strong <- subset(LAN_W , W0_X_LA  < -1)

##########   DRY Season 
ELN_D <- subset(L0_D, MEI_D > 0.5)
LAN_D <- subset(L0_D, MEI_D < -0.5)
NUT_Dx <- subset(L0_D, MEI_D > -0.5)
MEI_D2 <- NUT_Dx[,2]
NUT_D  <- subset(NUT_Dx,MEI_D2 < 0.5)

##########   Strong and Weak EL Dry Season 
D0_X_EL <- ELN_D[,2]
ELN_D_Weak <- subset(ELN_D , D0_X_EL  < 1)
ELN_D_Strong <- subset(ELN_D , D0_X_EL  > 1)

##########   Strong and Weak La Dry Season 
D0_X_LA <- LAN_D[,2]
LAN_D_Weak <- subset(LAN_D , D0_X_LA  > -1)
LAN_D_Strong <- subset(LAN_D , D0_X_LA  < -1)

##########   SUBSET 
EL_W_S <- ELN_W_Strong[,1]
EL_W_W <- ELN_W_Weak[,1]
LA_W_S <- LAN_W_Strong[,1]
LA_W_W <- LAN_W_Weak[,1]
NU_W <- NUT_W[,1]

EL_D_S <- ELN_D_Strong[,1]
EL_D_W <- ELN_D_Weak[,1]
LA_D_S <- LAN_D_Strong[,1]
LA_D_W <- LAN_D_Weak[,1]
NU_D <- NUT_D[,1]

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

##########   Mean

Me_EL_W_S <- round(mean(EL_W_S,na.rm=T),1) 
Me_EL_W_W <- round(mean(EL_W_W,na.rm=T),1) 
Me_LA_W_S <- round(mean(LA_W_S,na.rm=T),1)
Me_LA_W_W <- round(mean(LA_W_W,na.rm=T),1)
Me_NU_W <-   round(mean(NU_W,na.rm=T),1)
Me_EL_D_S <- round(mean(EL_D_S,na.rm=T),1)
Me_EL_D_W <- round(mean(EL_D_W,na.rm=T),1)
Me_LA_D_S <- round(mean(LA_D_S,na.rm=T),1)
Me_LA_D_W <- round(mean(LA_D_W,na.rm=T),1)
Me_NU_D <-   round(mean(NU_D,na.rm=T),1)

Cell.MEI[1:5,2] <- c( Me_EL_W_S, Me_EL_W_W, Me_NU_W,Me_LA_W_W, Me_LA_W_S)
Cell.MEI[1:5,3] <- c( Me_EL_D_S, Me_EL_D_S, Me_NU_D,Me_LA_D_W, Me_LA_D_S)

##########   MAX

Mx_EL_W_S <- round(max(EL_W_S,na.rm=T),1) 
Mx_EL_W_W <- round(max(EL_W_W,na.rm=T),1) 
Mx_LA_W_S <- round(max(LA_W_S,na.rm=T),1)
Mx_LA_W_W <- round(max(LA_W_W,na.rm=T),1)
Mx_NU_W <-   round(max(NU_W,na.rm=T),1)
Mx_EL_D_S <- round(max(EL_D_S,na.rm=T),1)
Mx_EL_D_W <- round(max(EL_D_W,na.rm=T),1)
Mx_LA_D_S <- round(max(LA_D_S,na.rm=T),1)
Mx_LA_D_W <- round(max(LA_D_W,na.rm=T),1)
Mx_NU_D <-   round(max(NU_D,na.rm=T),1)

Cell.MEI[1:5,4] <- c( Mx_EL_W_S, Mx_EL_W_W, Mx_NU_W,Me_LA_W_W, Mx_LA_W_S)
Cell.MEI[1:5,5] <- c( Mx_EL_D_S, Mx_EL_D_S, Mx_NU_D,Me_LA_D_W, Mx_LA_D_S)

##########   MIN

Mn_EL_W_S <- round(min(EL_W_S,na.rm=T),1) 
Mn_EL_W_W <- round(min(EL_W_W,na.rm=T),1) 
Mn_LA_W_S <- round(min(LA_W_S,na.rm=T),1)
Mn_LA_W_W <- round(min(LA_W_W,na.rm=T),1)
Mn_NU_W <-   round(min(NU_W,na.rm=T),1)
Mn_EL_D_S <- round(min(EL_D_S,na.rm=T),1)
Mn_EL_D_W <- round(min(EL_D_W,na.rm=T),1)
Mn_LA_D_S <- round(min(LA_D_S,na.rm=T),1)
Mn_LA_D_W <- round(min(LA_D_W,na.rm=T),1)
Mn_NU_D <-   round(min(NU_D,na.rm=T),1)

Cell.MEI[1:5,6] <- c( Mn_EL_W_S, Mn_EL_W_W, Mn_NU_W, Mn_LA_W_W, Mn_LA_W_S)
Cell.MEI[1:5,7] <- c( Mn_EL_D_S, Mn_EL_D_S, Mn_NU_D, Mn_LA_D_W, Mn_LA_D_S)

##########   DRY SEASON   

##########   Scaler
MAXA <- max(c(Mx_EL_D_S,Mx_EL_D_W,Mx_LA_D_S,Mx_LA_D_W,Mx_NU_D))
MAXB <- MAXA * 0.05
MAXB2 <- MAXA * 0.10
MAXC <- MAXA * 0.15
MAXD <- (MAXC * 1.2) + MAXA 
MAXD2 <- MAXD - MAXB2
MAXD22 <- MAXD - MAXB
MAXD3 <- MAXD - MAXC 

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u],"MEI_DRY.png"),width=5*dpi,height=5*dpi,res=dpi)

DRY_ENSO <- c(EL_D_S,EL_D_W,NU_D,LA_D_S,LA_D_W)

boxplot(c(EL_D_S),(EL_D_W),(NU_D),(LA_D_W),(LA_D_S),col=c("darkred","red","grey","blue","darkblue"),
        names = c("SEL" ,"WEL","NUT","WLA", "SLA"),ylab = paste0("Avg. Monthly Rainfall (",RFUnit2,")"),
        ylim=c(0,MAXD))


title(paste0("Dry Sea. RF: ",UNIT_Ns[u]), line = 1,cex.main = 1.1,col.main="gold4")

text(1,MAXD, paste0("Count = ",C_EL_D_S),cex=0.7)
text(1,MAXD2, paste0("Mean = ",Me_EL_D_S),cex=0.7)
text(1,MAXD3, paste0("Max = ",Mx_EL_D_S),cex=0.7)
text(1,MAXD22, paste0("Min = ",Mn_EL_D_S),cex=0.7)

text(2,MAXD, paste0("Count = ",C_EL_D_W),cex=0.7)
text(2,MAXD2, paste0("Mean = ",Me_EL_D_W),cex=0.7)
text(2,MAXD3, paste0("Max = ",Mx_EL_D_W),cex=0.7)
text(2,MAXD22, paste0("Min = ",Mn_EL_D_W),cex=0.7)

text(3,MAXD, paste0("Count = ",C_NU_D),cex=0.7)
text(3,MAXD2, paste0("Mean = ",Me_NU_D),cex=0.7)
text(3,MAXD3, paste0("Max = ",Mx_NU_D),cex=0.7)
text(3,MAXD22, paste0("Min = ",Mn_NU_D),cex=0.7)

text(4,MAXD, paste0("Count = ",C_LA_D_W),cex=0.7)
text(4,MAXD2, paste0("Mean = ",Me_LA_D_W),cex=0.7)
text(4,MAXD3, paste0("Max = ", Mx_LA_D_W),cex=0.7)
text(4,MAXD22, paste0("Min = ",Mn_LA_D_W),cex=0.7)

text(5,MAXD, paste0("Count = ",C_LA_D_S),cex=0.7)
text(5,MAXD2, paste0("Mean = ",Me_LA_D_S),cex=0.7)
text(5,MAXD3, paste0("Max = ",Mx_LA_D_S),cex=0.7)
text(5,MAXD22, paste0("Min = ",Mn_LA_D_S),cex=0.7)

dev.off() 

##########   Wet Season      
##########   Scaler

MAXA <- max(c(Mx_EL_W_S,Mx_EL_W_W,Mx_LA_W_S,Mx_LA_W_W,Mx_NU_W))
MAXB <- MAXA * 0.05
MAXB2 <- MAXA * 0.10
MAXC <- MAXA * 0.15
MAXD <- (MAXC * 1.2) + MAXA 
MAXD2 <- MAXD - MAXB2
MAXD22 <- MAXD - MAXB
MAXD3 <- MAXD - MAXC 

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u],"MEI_WET.png"),width=5*dpi,height=5*dpi,res=dpi)

boxplot(c(EL_W_S),(EL_W_W),(NU_W),(LA_W_W),(LA_W_S),col=c("darkred","red","grey","blue","darkblue"),
        names = c("SEL" ,"WEL","NUT","WLA", "SLA"),ylab=paste0("Avg. Monthly Rainfall (",RFUnit2,")"),
        ylim=c(0,MAXD))
title(paste0("Wet Sea RF: ",UNIT_Ns[u]), line = 1,cex.main = 1.1,col.main="darkgreen")

text(1,MAXD, paste0("Count = ",C_EL_W_S),cex=0.7)
text(1,MAXD2, paste0("Mean = ",Me_EL_W_S),cex=0.7)
text(1,MAXD3, paste0("Max = ",Mx_EL_W_S),cex=0.7)
text(1,MAXD22, paste0("Min = ",Mn_EL_W_S),cex=0.7)

text(2,MAXD, paste0("Count = ",C_EL_W_W),cex=0.7)
text(2,MAXD2, paste0("Mean = ",Me_EL_W_W),cex=0.7)
text(2,MAXD3, paste0("Max = ",Mx_EL_W_W),cex=0.7)
text(2,MAXD22, paste0("Min = ",Mn_EL_W_W),cex=0.7)

text(3,MAXD, paste0("Count = ",C_NU_W),cex=0.7)
text(3,MAXD2, paste0("Mean = ",Me_NU_W),cex=0.7)
text(3,MAXD3, paste0("Max = ",Mx_NU_W),cex=0.7)
text(3,MAXD22, paste0("Min = ",Mn_NU_W),cex=0.7)

text(4,MAXD, paste0("Count = ",C_LA_W_W),cex=0.7)
text(4,MAXD2, paste0("Mean = ",Me_LA_W_W),cex=0.7)
text(4,MAXD3, paste0("Max = ", Mx_LA_W_W),cex=0.7)
text(4,MAXD22, paste0("Min = ",Mn_LA_W_W),cex=0.7)

text(5,MAXD, paste0("Count = ",C_LA_W_S),cex=0.7)
text(5,MAXD2, paste0("Mean = ",Me_LA_W_S),cex=0.7)
text(5,MAXD3, paste0("Max = ",Mx_LA_W_S),cex=0.7)
text(5,MAXD22, paste0("Min = ",Mn_LA_W_S),cex=0.7)

dev.off()

write.csv(Cell.MEI,paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u]," MEI_A.csv"),row.names = F)
#}

#####################################################
##### Air temperature graph

maps<-AT_Map_Path_A

# make list of years
years<-as.list(list.files(maps))

# create empty table for raster stats
table<-data.frame(matrix(ncol=5,nrow=0, 
                         dimnames=list(NULL, c("year","month", "min", "max", "mean"))))

# set row
no<-1

### for each year, loop through months and extract min, max, mean air temp value
### from within the study area

for (i in years) {

  setwd(paste0(maps,i,"/"))
  
  # make list of monthly files
  f<-as.list(list.files())
  
  # loop through monthly files and fill in table
  for (x in f) {
 
    # set year and month
    y<-substr(x, nchar(x) - 10,nchar(x))
    table[no,]$year<-substr(y,1,4)
    m<-substr(x, nchar(x) - 5,nchar(x))
    table[no,]$month<-substr(m,1,2)
    
    # load raster
    map<-raster(x)
    
    # crop to study area
    m2 <- crop(map, extent(HALE))
    m3 <- mask(m2, HALE)
    
    # calculate min, max, mean air temp
    table[no,]$min<-round(cellStats(m3, stat="min"), digits=2)
    table[no,]$max<-round(cellStats(m3, stat="max"), digits=2)
    table[no,]$mean<-round(cellStats(m3, stat="mean"), digits=2)
    
    no<-no+1
  }
}

table

# write to csv
write.csv(table, paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u],"_monthly_airtemp.csv"))

###### Monthly air temperature time series
dat<-table

head(dat)
tail(dat)

# convert all celcius to farenheit
dat$min<-(dat$min*(9/5)) + 32
dat$max<-(dat$max*(9/5)) + 32
dat$mean<-(dat$mean * (9/5)) + 32

# date column
dat$date<-as.Date(paste0(dat$year,"-",dat$month,"-01"))

##### add annual mean values (from code below this one)
dat2<-merge(dat,dat.y[,c("year","mean")], by="year", all.x=T)
head(dat2, 20)

# set y-axis limits
ylow<-min(dat2$mean.x)*.95
yhi<-max(dat2$mean.x)*1.0005

# set slope
slope<-formatC((coef(lm(dat2$mean.x~dat2$date))[2]), format="e", digits=2)
slope

dpi=300

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u],"_monthly_airtemp.png"),width=6*dpi,height=3*dpi,res=dpi)

ggplot(dat2, aes(x=date,y=mean.x)) +
  geom_line(color="grey60") +
  geom_smooth(span=0.2, se=F, size=1.3, color="orange") +
  geom_smooth(method=lm, se=F, color="black") +
  stat_cor(method="pearson", label.x=as.Date("2012-01-01"), label.y=58) +
  scale_x_date(date_breaks = "4 years", labels = date_format(format="%Y")) +
  ylim(ylow,yhi) +
  labs(title="Monthly Air Temperature",
       y="Temperature (F)", x="") +
  geom_hline(yintercept=0) +
  annotate("text", x=as.Date("2015-07-01"), y=56.8, 
           label=paste0("Slope = ",slope)) +
  theme(panel.background=element_rect(fill=NA, color="black"),
        panel.grid.major=element_line(color="grey90"),
        panel.grid.minor=element_blank())

dev.off()

#############################
### Annual air temp trend

# aggregate monthly to annual air temp
dat.y<-aggregate(mean ~ year, dat, mean)
dat.y$mean<-round(dat.y$mean, digit=2)
dat.y$meanmin<-round((aggregate(mean ~ year,dat, min))$mean, digits=2)
dat.y$meanmax<-round((aggregate(mean ~ year,dat, max))$mean, digits=2)
dat.y$min<-round((aggregate(min ~ year, dat, min))$min, digits=2)
dat.y$max<-round((aggregate(max ~ year, dat, max))$max, digits=2)

head(dat.y)

# make date column
dat.y$date<-as.Date(paste0(dat.y$year,"-01-01"))

slope<-formatC((coef(lm(dat.y$mean~dat.y$date))[2]), format="e", digits=2)
slope

png(paste0(RFOLDER,UNIT_N[u],"/",UNIT_N[u],"_annual_airtemp.png"),width=6*dpi,height=3*dpi,res=dpi)

ggplot(dat.y, aes(x=date, y=mean)) +
  geom_line(size=1.2, color="orange") +
  geom_smooth(method=lm, se=F, size=1, color="black") +
  stat_cor(method="pearson", label.x=as.Date("1991-01-01"), label.y=58) +
  scale_x_date(date_breaks = "2 years", labels = date_format(format="%Y")) +
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




#END OF LOOP 
write.csv(Cell.DataCL,paste0(SFOLDER,UNIT_N[u],"_CL.csv"),row.names = F)

Cell.RF_Year[2] <- Cell.DataCL[3]
write.csv(Cell.RF_Year,paste0(SFOLDER,UNIT_N[u],"_RF.csv"),row.names = F)

#write.csv(Cell.MEI_All,paste0(SFOLDER,Shape_Name,"_MEI.csv"),row.names = F)
###  END  #############


