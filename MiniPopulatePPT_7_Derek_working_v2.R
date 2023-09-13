library(magrittr)
library(tidyverse)
library(rvg)
library(dplyr)
library(ggplot2)
library(magrittr)
library(knitr)
library(xtable)
library(flextable)
library(officer)
library(mschart)
library(purrr)
library(imager)
library(pdftools)

setwd("F:/PDKE/CCVD/MINI_Phase2/")
#Path = "C:/Users/BUNNY1/Documents/Work From Home/PAPERS/Beef_Production/Specific Ranch/"
I_FOLDER <- "F:/PDKE/CCVD/IMAGE/"        # Folder with images 
R_FOLDER <- "F:/PDKE/CCVD/MINI_Phase2/"  # File with your site specific folders 
P_FOLDER <- "F:/PDKE/CCVD/MINI_PPT/"     # Output folder 
RANL <- list.files(R_FOLDER)
NF <- length(RANL)

### SET SITE FOLDER AND VERSION ###
RANL
# LOOP 
f<-1
RANL[f]

# VERSION
ver<-5

#

#TUnit = "\u00B0C"
#TUnit2 = " \u00B0C"
#RFUnit = " mm"
#RFUnit2 = "mm"
#RFUnit3 = " millimeters"
#ELUnit = " m"
#ELUnit2 = "m"

TUnit = "\u00B0F"
TUnit2 = " \u00B0F"
RFUnit = " in"
RFUnit2 = "in"
RFUnit3 = "°F"
ELUnit = " ft"
ELUnit2 = "ft"

############################ LOGOS 
## Read in images from the Image folder 

PDKE_Short <- paste0(I_FOLDER,"PDKE_Logo_Color_Rounded_Type-03.jpg")
PDKE_Short
PDKE_S <- external_img(src = PDKE_Short, height = 1,width = 1) 

PDKE_Long <- paste0(I_FOLDER,"PDKE_Logo_Color_Horizontal_Type-05.png")
PDKE_L <- external_img(src = PDKE_Long, height = 1,width = 1) 

EWCfile <- paste0(I_FOLDER,"EWC No Boarder.png")
EWCimg <- external_img(src = EWCfile, height = 1,width = 1) 

EWClfile <- paste0(I_FOLDER,"ewc_long.png")
EWClimg <- external_img(src = EWClfile, height = 1,width = 1)

CASCfile <- paste0(I_FOLDER,"PICASC_NoBoarder.jpg")
CASCimg <- external_img(src = CASCfile , height = 1, width = 1) 

FSfile <- paste0(I_FOLDER,"FS.png")
FSimg <- external_img(src = FSfile, height = 1, width = 1) 

UHfile <- paste0(I_FOLDER,"UH.png")
UHimg <- external_img(src = UHfile, height = 1, width = 1) 

USGSfile <- paste0(I_FOLDER,"USGS.png")
USGSimg <- external_img(src = USGSfile, height = 1, width = 1) 

# RISAfile <- paste0(I_FOLDER,"PacificRISA.jpg.png")
# RISAimg <- external_img(src = RISAfile, height = 1, width = 1) 

NIDISfile <- paste0(I_FOLDER,"NIDIS.png")
NIDISimg <- external_img(src = NIDISfile, height = 1, width = 1) 

NOAAfile <- paste0(I_FOLDER,"NOAA.png")
NOAAimg <- external_img(src = NOAAfile, height = 1, width = 1) 

HAVOfile <- paste0(I_FOLDER,"HAVO.png")
HAVOimg <- external_img(src = HAVOfile, height = 1, width = 1) 

USGSfile <- paste0(I_FOLDER,"USGS.png")
USGSimg <- external_img(src = USGSfile, height = 1, width = 1) 

DOFAWfile <- paste0(I_FOLDER,"DOFAW.png")
DOFAWimg <- external_img(src = DOFAWfile, height = 1, width = 1) 

DLNRfile <- paste0(I_FOLDER,"DLNR.jpg")
DLNRimg <- external_img(src = DLNRfile, height = 1, width = 1)

DPlanfile <- paste0(I_FOLDER,"DPlan.png")
DPlanimg <- external_img(src = DPlanfile, height = 1, width = 1) 

PIDPfile <- paste0(I_FOLDER,"PIDP.png")
PIDPimg <- external_img(src = PIDPfile, height = 1, width = 1)

CLfile <- paste0(I_FOLDER,"Clark.png")
CLimg <- external_img(src = CLfile, height = 1, width = 1)

############################ Pictures

WRRCfile <- paste0(I_FOLDER,"WRRC-200x200.png")
WRRCimg <- external_img(src = WRRCfile, height = 1, width = 1) 

### original "Part 1: Climate Characteristics" images
    SUNfile <- paste0(I_FOLDER,"SUN.jpg")
    SUNimg <- external_img(src = SUNfile, height = 1,width = 1) 
    
    RFfile <- paste0(I_FOLDER,"RF.jpg")
    RFimg <- external_img(src = RFfile , height = 1, width = 1) 
    
    RBfile <- paste0(I_FOLDER,"Rainbow.jpg")
    RBimg <- external_img(src = RBfile, height = 1,width = 1) 
    
### Derek's "Part 1: Climate Characteristics" images
    RGfile <- paste0(I_FOLDER,"climstation.jpg")
    RGimg <- external_img(src = RGfile, height = 1,width = 1) 

    RFfile <- paste0(I_FOLDER,"Mean Annual Rainfall.jpg")
    RFimg <- external_img(src = RFfile, height = 1,width = 1) 
    
    MODfile <- paste0(I_FOLDER,"modis.jpg")
    MODimg <- external_img(src = MODfile, height = 1,width = 1) 

### Derek's slide 6.1 "Landcover"
    LClegfile <- paste0(I_FOLDER,"LCMAP_legend.jpg")
    LClegimg <- external_img(src = LClegfile, height = 2, width = 1)
    
Theromfile <- paste0(I_FOLDER,"therom.jpg")
Theromimg <- external_img(src = Theromfile , height = 1, width = 1) 

RAtlasfile <- paste0(I_FOLDER,"RFAtlas.png")
RAtlasimg <- external_img(src = RAtlasfile, height = 1,width = 1) 

CoHfile <- paste0(I_FOLDER,"CoH.png")
CoHmimg <- external_img(src = CoHfile , height = 1, width = 1) 

Droughtfile <- paste0(I_FOLDER,"Drought.jpg")
Droughtimg <- external_img(src = Droughtfile , height = 1, width = 1) 

Firefile <- paste0(I_FOLDER,"Fire.jpg")
Fireimg <- external_img(src = Firefile , height = 1, width = 1) 

Downfile <- paste0(I_FOLDER,"Down.png")
Downimg <- external_img(src = Downfile , height = 1, width = 1)

D_AGfile <- paste0(I_FOLDER,"AG_D.png")
D_AGimg <- external_img(src = D_AGfile , height = 1, width = 1)

D_METfile <- paste0(I_FOLDER,"MET_D.png")
D_METimg <- external_img(src = D_METfile , height = 1, width = 1)

D_HYDfile <- paste0(I_FOLDER,"HYD_D.png")
D_HYDimg <- external_img(src = D_HYDfile , height = 1, width = 1)

D_SOCfile <- paste0(I_FOLDER,"SOC_D.png")
D_SOCimg <- external_img(src = D_SOCfile , height = 1, width = 1)

D_ECOfile <- paste0(I_FOLDER,"ECO_D.png")
D_ECOimg <- external_img(src = D_ECOfile , height = 1, width = 1)

D_fivefile <- paste0(I_FOLDER,"types_of_drought2.png")
D_fiveimg <- external_img(src = D_fivefile , height = 1, width = 1)

#media
D_AGfileM <- paste0(I_FOLDER,"Media1_SPI.MP4")
D_AGimgM <- external_img(src = D_AGfileM , height = 1, width = 1)

#ahupuaa image
D_AHfile <- paste0(I_FOLDER,"ahupuaa_art.png")
D_AHimg <- external_img(src = D_AHfile , height = 1, width = 1)

#Hawaiian land division figure boxes
D_HLDfile <- paste0(I_FOLDER,"HLD_fig.jpg")
D_HLDimg <- external_img(src = D_HLDfile , height = 1, width = 1)

######


#for(f in 1:NF) {
#for(f in 79:82) {

#create path to folder 
RAN_F <- paste0(R_FOLDER,RANL[f],"/")
# list the output from Code 1 
RANL2 <- list.files(RAN_F)
#Unit Name
NameF <- basename(RAN_F)

##### Read in CSV Files 
## Read in Mean Climate File 
CLIM <- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF,"MEAN Climate.csv"),sep=",")
CLIM
  #Short Name
  SNameF <- CLIM[1,16]
  SNameF

## Read in landcover file
LAND <- read.csv(paste0(R_FOLDER,NameF,"/",NameF," Landcover.csv"),sep=",")

## Hawaiian land division files
MOKU <- read.csv(paste0(R_FOLDER,NameF,"/",NameF," Moku.csv"),sep=",")
AHU <- read.csv(paste0(R_FOLDER,NameF,"/",NameF," Ahupuaa.csv"),sep=",")
AHU5 <- read.csv(paste0(R_FOLDER,NameF,"/",NameF," Ahupuaa.csv"),sep=",")
  
#Font styles 
  #Suggetions Type_Color_Size  ex B_DR_40
  
  fp_BR <- fp_text(bold = TRUE, color = "orangered3", font.size = 40)
  fp_BR2 <- fp_text(bold = TRUE, color = "black", font.size = 40)
  fp_BR3 <- fp_text(bold = TRUE, color = "darkgreen", font.size = 40)
  fp_NM <- fp_text(bold = TRUE, font.size = 20,underlined = TRUE)
  fp_NMS <- fp_text(bold = TRUE, font.size = 8,underlined = TRUE)
  fp_NM2 <- fp_text(font.size = 18, color = "black")
  fp_NM2a <- fp_text(font.size = 17, color = "black")
  fp_NM3 <- fp_text(font.size = 18, color = "darkgreen")
  fp_NM4 <- fp_text(font.size = 18, color = "darkred")
  fp_NM5 <- fp_text(font.size = 18, color = "darkred",underlined = TRUE,bold = TRUE)
  fp_NM6<- fp_text(font.size = 18, color = "black",bold = TRUE)
  fp_NM7<- fp_text(font.size = 18, color = "black")
  fp_NM8<- fp_text(font.size = 16, color = "black", bold = TRUE)
  fp_Fig <- fp_text(font.size = 12, color = "black")
  fp_Fig2 <- fp_text(font.size = 10, color = "black")
  fp_Fig3 <- fp_text(font.size = 14, color = "black")
  fp_Fig4<- fp_text(bold = TRUE, color = "darkred",font.size = 15,underlined = TRUE)
  fp_Fig5 <-fp_text(font.size=11, color = "black")
  fp_ft <- fp_text(font.size=12, font.family = "Calibri (Body)", color = "darkgrey")
  
  FTXTT <-  fp_text(bold = TRUE,color = "darkblue", font.size = 40)
  FTXTT2 <- fp_text(italic = TRUE, color = "black", font.size = 40)
  FTXTT3 <- fp_text(italic = TRUE, color = "darkred", font.size = 40)
  FTXTT4 <- fp_text(color = "black", font.size = 38)
  
  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 20)
  
  #May or may not be used
  fp_TxB <- fp_text(italic = TRUE, color = "black", font.size = 20) 
  fp_TxBs <- fp_text(italic = TRUE, color = "black", font.size = 15) 
  
  
  ############# Slide 1 TITLE 
  ISL<- CLIM[1,15]
  
  #Slide one Title 
  #Standardize S1_TIT
  TIT <- fpar(ftext("Climate Change, Climate Variability, & Drought Portfolio", fp_BR),fp_p = fp_par(text.align = "center") )
  SUB <- fpar(ftext(paste0(NameF,", ",ISL),fp_BR2),fp_p = fp_par(text.align = "center") )
  
################ Slide 2 PDKE
  # start slide/page numbers 
  p<-1
  
  p<-p+1
  p1<-p
  #S2_TIT<- block_list(
    #fpar(ftext("Pacific Drought Knowledge Exchange", FTXTT),fp_p = fp_par(text.align = "center")))
  
  # Determine if the AOI full and short name are different and use the correct text
  if((nchar(NameF))>(nchar(SNameF))) {
  HDKE<- paste0("Climate change, climate variability, and drought (CCVD) will exert a growing impact on Hawaii's ecosystems, ",
                "agriculture and communities in the future. While resource managers are tasked with preparing for this with the ",
                "best available information, it is hard to know what data, research and recommendations are available. ",
                "The Pacific Drought Knowledge Exchange (PDKE) program ",
                "focuses on facilitating knowledge exchange between the research community and resource managers and stakeholders, ",
                "thereby expanding the utility of climate and drought-related scientific products.

This CCVD portfolio is a comprehensive synthesis of climate and drought information developed specifically for ", NameF," (",SNameF,"). ",
                "It is designed to provide relevant climate and drought information needed to inform land management and guide future ",
                "research and extension. While we try to include a wide range of useful site-specific data products, we also recognize ",
                "that every site is unique and PDKE is happy to collaborate on producing additional drought products beyond the CCVD portfolio ",
                "to meet stakeholder needs.
                
The PDKE program was piloted in November of 2019 with funding from the Pacific Islands Climate Adaptation Science Center (PICASC). ",
"Subsequent PDKE activities and updates to the CCVD portfolio have been funded by PICASC, the East-West Center, ",
"and the National Integrated Drought Information System (NIDIS). ") }

  if((nchar(NameF)) == (nchar(SNameF))) {
    HDKE<- paste0("Climate change, climate variability, and drought (CCVD) will exert a growing impact on Hawaii's ecosystems, ",
                  "agriculture and communities in the future. While resource managers are tasked with preparing for this with the ",
                  "best available information, it is hard to know what data, research and recommendations are available. ",
                  "The Pacific Drought Knowledge Exchange (PDKE) program ",
                  "focuses on facilitating knowledge exchange between the research community and resource managers and stakeholders, ",
                  "thereby expanding the utility of climate and drought-related scientific products.

This CCVD portfolio is a comprehensive synthesis of climate and drought information developed specifically for ", NameF,". ",
                  "It is designed to provide relevant climate and drought information needed to inform land management and guide future ",
                  "research and extension. While we try to include a wide range of useful site-specific data products, we also recognize ",
                  "that every site is unique and PDKE is happy to collaborate on producing additional drought products beyond the CCVD portfolio ",
                  "to meet stakeholder needs.
                
The PDKE program was piloted in November of 2019 with funding from the Pacific Islands Climate Adaptation Science Center (PICASC). ",
                  "Subsequent PDKE activities and updates to the CCVD portfolio have been funded by PICASC, the East-West Center, ",
                  "and the National Integrated Drought Information System (NIDIS). ") }

# Text for slide as a variable 
fp_HDKE <- fpar(ftext(HDKE, fp_TxBs))


# ###############   Slide 3 CCVD
#   
#   #not sure why I use "block list" check this. 
# 
#   S3_TIT<- block_list(
#   fpar(ftext("CCVD Portfolio", FTXTT),fp_p = fp_par(text.align = "center")))
# 
#   CCVD<- paste0("The climate change, climate variability, and drought (CCVD) portfolio is a comprehensive synthesis ",
#               "of climate and drought information developed specifically for ", NameF, " (",SNameF,"). This portfolio is designed to ",
#               "provide both researchers and land managers with relevant climate and drought information ",
#               "needed to inform land management and guide future research and extension.") 
# 
#   fp_CCVD <- fpar(ftext(CCVD, fp_Tx))
# 
# #Map Figure 
# MAPfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Map.png")
# MAPimg <- external_img(src = MAPfile) 
# 
# plot(load.image(MAPimg))
# FIG_1 <- block_list(fpar(ftext(paste0("Figure 1. Map of ",ISL," with ",SNameF," in red."), fp_Fig)))

################ Slide 3 PART 1
  p<-p+1
  p2<-p
  S3_TIT<- block_list(
  fpar(ftext("Part 1: Describing the Area", FTXTT),fp_p = fp_par(text.align = "center")))

  CCVD<- paste0("In describing any area of management in Hawaii, it is important to present both traditional and contemporary knowledge. ",
                "Traditional Hawaiian landscape divisions are well-documented and were established largely following geological features and ",
                "the natural flow of resources throughout the landscape. 

This portfolio provides a brief description of ",SNameF," in context of traditional Hawaiian landscapes, as well as contemporary knowledge on ",
                "elevation and current landcover.")

  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 19) 
  fp_CCVD <- fpar(ftext(CCVD, fp_Tx))
  FIG_1 <- block_list(fpar(ftext(paste0("Credit: Matt Foster"), fp_Fig2)))
  

################ Slide 4 Hawaiian Land Divisions
  p<-p+1
  p3<-p
  
  S4a_TIT<- block_list(
    fpar(ftext("Hawaiian Land Divisions", FTXTT),fp_p = fp_par(text.align = "center")))
  
  HLD1<- paste0("There are three types of traditional Hawaiian landscape divisions available as GIS layers and ",
  "presented here for ",SNameF,". The land divisions are shown in red and ",SNameF," is shown in blue.")
  
  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 19) 
  fp_HLD1 <- fpar(ftext(HLD1, fp_Tx))
  
  # Make moku label based on how many moku there are
  MO<-MOKU[1,]$moku
  if(nrow(MOKU)==1) {MO1<-MO}
  if(nrow(MOKU)==2) {MO1<-paste0(MO, " and ",MOKU[2,]$moku)}
  if(nrow(MOKU)>2) {MO1<-paste0(MO,", ",MOKU[2,]$moku," and ",MOKU[3,]$moku)}
  MO1
  
  # Ahupuaa count and names
  AHc<-nrow(AHU)
  if(AHc==1) {AHa<-paste0(AHU5[1,]$ahupuaa2,".")}
  if(AHc==2) {AHa<-paste0(AHU5[1,]$ahupuaa2," and ",AHU5[2,]$ahupuaa2,".")}
  if(AHc>2) {AHa<-paste0(AHU5[1,]$ahupuaa2,"; ",AHU5[2,]$ahupuaa2,"; and ",AHU5[3,]$ahupuaa2,".")}
  AHa

  # Figure text
  MP <- block_list( 
    fpar(ftext(paste0("Mokupuni is the largest land division and refers to entire islands. ",SNameF," is in the mokupuni of ",ISL,"."), fp_Fig5)))
  MP
  
  if(nrow(MOKU)==1){MO <- block_list(
    fpar(ftext(paste0("Within mokupuni are smaller divisions called moku. ",SNameF," is in the moku of ",MO1,"."), fp_Fig5)))}  
  if(nrow(MOKU)>1){MO <- block_list(
    fpar(ftext(paste0("Within mokupuni are smaller divisions called moku. ",SNameF," is primarily in the moku of ",MO1,"."), fp_Fig5)))}
  MO
  
  if(AHc<3) {AH <- block_list(
      fpar(ftext(paste0("Within each moku are several ahupuaa which commonly extend from uplands to the sea. ",SNameF,
                        " is situated within ",AHc," ahupuaa -  ",AHa), fp_Fig5)))}
  if(AHc>2) {AH <- block_list(
    fpar(ftext(paste0("Within each moku are several ahupuaa which commonly extend from uplands to the sea. ",SNameF,
    " is situated within ",AHc," ahupuaa - primarily ",AHa), fp_Fig5)))}
  AH
  
  # Load land division maps
  MPmfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Mokupuni.png") 
  MPmimg <- external_img(src = MPmfile)
  
  MOmfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Moku.png") 
  MOmimg <- external_img(src = MOmfile)
  
  AHmfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Ahupuaa.png") 
  AHmimg <- external_img(src = AHmfile)
  
  ft1 <- block_list(fpar(ftext(paste0("https://planning.hawaii.gov/gis"), fp_ft)))

###############   Slide 5 Elevation
  p<-p+1
  p4<-p
  
S5_TIT<- block_list(
  fpar(ftext("Elevation", FTXTT),fp_p = fp_par(text.align = "center")))

#Total land area (acres)
ar<-round(sum(LAND$acres),0)

#Pulling from table that was read in
El_Mean <- round(CLIM[1,2],0)
El_Max<-  round(CLIM[2,2],0)
El_Min <-  round(CLIM[3,2],0)
#calculat Elevation Difference
EL_Dif <- El_Max - El_Min

ELEV_T<- paste0(SNameF," is located on the Island of ", ISL, ". It covers ",ar," acres and a vertical elevation range of ", EL_Dif, ELUnit, ". In Hawaii, climate gradients can change ",
                      "significantly over short distances due to changes in elevation, topography, and orientation to the prevailing winds.")

M_Elev <-  block_list(
  fpar(ftext(ELEV_T, fp_Tx)),
  fpar(ftext(" ", fp_NM2)),
  fpar(ftext(paste0("Elevation ", SNameF), fp_NM5)),
  fpar(ftext       ("          ", fp_NMS)),
  fpar(ftext(paste0("Minimum = ", El_Min, ELUnit), fp_NM3)),
  fpar(ftext(paste0("Mean = ", El_Mean, ELUnit), fp_NM3)),
  fpar(ftext(paste0("Maximum = ", El_Max, ELUnit), fp_NM3)))
  
# Figure 2 caption 
FIG_2 <- block_list(
  fpar(ftext(paste0("Figure 2. Elevation for the Island of ",ISL," with ",SNameF," outlined in black."), fp_Fig)))

#Read in elevation figure 
    Elevfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," ELMap.png")
    Elevimg <- external_img(src = Elevfile) 
    
ft2 <- block_list(fpar(ftext(paste0("https://www.pacioos.hawaii.edu/"), fp_ft)))
    

############### Slide 6.1 (NEW) Landcover
p<-p+1
p5<-p

S6.1_TIT<- block_list(
  fpar(ftext("Landcover", FTXTT),fp_p = fp_par(text.align = "center")))

#Pulling from table that was read in
LC1 <- LAND[1,5]
LC2 <- LAND[2,5]
LC3 <- LAND[3,5]

#Body text
LAND_T<- paste0("The three most common types of landcover in ",SNameF," are ",LC1,", ",LC2,", and ",LC3,". Each landcover type exhibits different climate change impacts and management needs.")
LAND_T

M_LAND <-  block_list(
  fpar(ftext(LAND_T, fp_Tx)))

#Read in landcover figures 
LCbarfile <- paste0(R_FOLDER,NameF,"/",NameF," LC_barchart.png")
LCbarimg <- external_img(src = LCbarfile, height = 3,width = 3)

LCmapfile <- paste0(R_FOLDER,NameF,"/",NameF," LCMap.png")
LCmapimg <- external_img(src = LCmapfile)

# Figure 3.1 caption 
FIG_3.1 <- block_list(
  fpar(ftext(paste0("Figure 3. Bar graph showing amount and percent of each landcover type within ",SNameF,"."), fp_Fig)))

# Figure 3.2 caption
FIG_3.2 <- block_list(
  fpar(ftext(paste0("Figure 4. Landcover mapping for the island of ", ISL," with ",SNameF," outlined in red. ",
                    "The maps shown in the following slides will be for the ",SNameF," area only."), fp_Fig)))

ft3 <- block_list(fpar(ftext(paste0("https://www.usgs.gov/special-topics/lcmap"), fp_ft)))

################ Slide 7.1 Water Sources - Aquifers
p<-p+1
p6<-p

S6.2_TIT<- block_list(
  fpar(ftext("Water Sources", FTXTT),fp_p = fp_par(text.align = "center")))

#Aquifer spreadsheet
AQ <- read.csv(paste0(R_FOLDER,NameF,"/",NameF," Aquifer.csv"),sep=",")
AQc<-nrow(AQ)

AQ_T<-flextable(AQ)
AQ_T<-bold(AQ_T, bold=TRUE, part="header")
AQ_T<-fontsize(AQ_T, size=12)
AQ_T<-autofit(AQ_T)
AQ_T

#Body text
AQ_Tt<- paste0("There are ",AQc," aquifers in ",SNameF,", their characteristics are listed below.

In general, basal aquifers are more susceptible to saltwater intrusion than high level aquifers.")
AQ_Tt

M_AQ <-  block_list(
  fpar(ftext(AQ_Tt, fp_Tx)))
M_AQ

# map
AQmapfile <- paste0(R_FOLDER,NameF,"/",NameF," Aquifers.png")
AQmapimg <- external_img(src = AQmapfile)

# Figure caption
FIG_5.b <- block_list(
  fpar(ftext(paste0("Figure 5. Department of Health aquifer mapping labeled by DOH Aquifer number ",
                    "and colored by hydrology type. ",SNameF," outlined in red."), fp_Fig)))
FIG_5.b

# Table caption
TAB_0 <- block_list(
  fpar(ftext(paste0("Table 1. Aquifer characteristics for ", SNameF,"."), fp_Fig)))
TAB_0

ft3.1 <- block_list(fpar(ftext(paste0("https://geoportal.hawaii.gov/datasets"), fp_ft)))

################ Slide 7.2 Water Sources - Streams
p<-p+1
p7<-p

S6.2_TIT<- block_list(
  fpar(ftext("Water Sources", FTXTT),fp_p = fp_par(text.align = "center")))

#Hydrologic features spreadsheet
HF <- read.csv(paste0(R_FOLDER,NameF,"/",NameF," Hydro_features.csv"),sep=",")
HFc<-nrow(HF)

# make text string naming all feature types within AOI
if(HFc>0) {ft<-paste0(HF[1,])}
if(HFc==0) {ft<-"no"}

for (i in 2:nrow(HF)) {
  t<-HF[i,]
  if(i<nrow(HF)) {ft<-paste0(ft,", ",t)}
  if(i==nrow(HF)) {ft<-paste0(ft," and ",t)}  
}
ft

#Body text
HF_T<- paste0("There are ",ft," hydrologic features within ",SNameF,".

Perennial streams are typically reliable water sources, while intermittent ",
"streams stop flowing during dry periods. Well-managed waterways and canal/ditch ",
"systems can reduce the impacts of drought and flooding.")

M_HF <-  block_list(
  fpar(ftext(HF_T, fp_Tx)))
M_HF

# map
HFmapfile <- paste0(R_FOLDER,NameF,"/",NameF," Streams.png")
HFmapimg <- external_img(src = HFmapfile)

# legend
HFlmapfile <- paste0(I_FOLDER,"water_features_leg.jpg")
HFlmapimg <- external_img(src = HFlmapfile)

# Figure caption
FIG_6.b <- block_list(
  fpar(ftext(paste0("Figure 6. Streams and hydrologic features with ",
                    SNameF," outlined in red."), fp_Fig)))

# ft3.1 <- block_list(fpar(ftext(paste0("https://geoportal.hawaii.gov/datasets"), fp_ft)))


################ Slide 7 PART 1
p<-p+1
p8<-p

P1_TIT<- block_list(
  fpar(ftext("Part 2: Climate Characteristics", FTXTT)),
  fpar(ftext(SNameF, FTXTT3),fp_p = fp_par(text.align = "center")))

PART1<- paste0("In developing this Portfolio, we relied on several gridded climate products available for the State of Hawaii. Annual and monthly estimates of rainfall were ",
               "obtained from the Hawaii Climate Data Portal (HCDP). Gridded estimates of other climate variables were ",
               "obtained from the UH Manoa Climate of Hawaii data page. ",
               "We retrieved all the data points that fell within the boundaries of ",SNameF," from our 250 meter resolution state-wide maps to support the presented analyses.")

fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 19) 
fp_PART1 <- fpar(ftext(PART1, fp_Tx))

#NameF_I_Tx <- fp_text(italic = TRUE, color = "darkblue", font.size = 20) 
#NameF_I <- fp_text(italic = TRUE, color = "darkblue", font.size = 20) 

###############   Rainfall Stations Slide 8a
p<-p+1
p9<-p

S8.0_TIT<- block_list(
  fpar(ftext("Rainfall Station Locations", FTXTT),fp_p = fp_par(text.align = "center")))

### Bring in table of stations and network links
RS<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," rain_stations.csv"),sep=",")
RS_T <- flextable(RS)
RS_T <- bold(RS_T, bold = TRUE, part = "header")
RS_T <- fontsize(RS_T, size = 12)
#Creates a table for the PPT
RS_T <- autofit(RS_T)
RS_T

# get closest station and network
s1<-RS[1,1]
snn<-RS[1,2]

# set island name to Big Island if it's ... on Big Island
if(ISL == c("Hawaii")) {ISL2<-c("Big Island")}
if(ISL != c("Hawaii")) {ISL2<-ISL}

RS1<- paste0("Rainfall data for this portfolio were estimated based on measurements made ",
             "by hundreds of stations across the state.

The closest station to ",SNameF," is ",s1,", part of the ",snn," climate station network.")

fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 19) 
fp_RS1 <- fpar(ftext(RS1, fp_Tx))
fp_RS1

# Figure 5 caption 
FIG_5.a <- block_list(
  fpar(ftext(paste0("Figure 6. Rainfall station locations across ",ISL2," with ",SNameF,
                    " outlined in blue and three closest stations in orange."), fp_Fig)))

# Figure 6 caption 
FIG_6.a <- block_list(
  fpar(ftext(paste0("Figure 7. Three closest stations to ",SNameF," with links to more ",
                    "information on the station network. If there are more than three stations ",
                    "at the site, only three will be listed here."), fp_Fig)))

# get map of stations and AOI
RSfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," rf_stations.png")
RSimg <- external_img(src = RSfile) 

###############   Slide 6 Mean Climate
p<-p+1
p10<-p

    S6_TIT<- block_list(
      fpar(ftext("Annual Climate Characteristics", FTXTT),fp_p = fp_par(text.align = "center")))
    
    ANNCLIM_T<- paste0("Climatic conditions in Hawaii can vary greatly across the landscape. ",
    "These maps show the variation in select climate variables across ",SNameF," and the table below has min. and max. values taken from the maps.")
    M_ANNCLIM <-  block_list(
      fpar(ftext(ANNCLIM_T, fp_Tx)))
    
    # Make a table of min and maxes
    T<-data.frame(matrix(0, ncol=3, nrow=7))
    T$X1<-c("Rainfall (in.)","Air Temperature (°F)","Relative Humidity (%)",
             "Solar Radiation (W/m2)","Soil Moisture (Ratio)","Evapotranspiration (in.)","Windspeed (mph)")
    T[1,]$X2<-CLIM[3,3]
    T[1,]$X3<-CLIM[2,3]
    T[2,]$X2<-CLIM[3,4]
    T[2,]$X3<-CLIM[2,4]
    T[3,]$X2<-CLIM[3,7]
    T[3,]$X3<-CLIM[2,7]
    T[4,]$X2<-CLIM[3,9]
    T[4,]$X3<-CLIM[2,9]
    T[5,]$X2<-CLIM[3,8]
    T[5,]$X3<-CLIM[2,8]
    T[6,]$X2<-CLIM[3,10]
    T[6,]$X3<-CLIM[2,10]
    T[7,]$X2<-CLIM[3,12]
    T[7,]$X3<-CLIM[2,12]
    colnames(T)<-c("Climate Variable","Min","Max")

    # format table
    Tt<-flextable(T)
    Tt<-bold(Tt, bold = TRUE, part = "header")
    # Tt<-width(Tt, j=1, 2, unit="in")

    #Creates a table for the PPT
    Tt <- autofit(Tt, add_h=-0.5, part=c("body","header"), unit = "in")
    Tt
 
    FIG_5.1 <- block_list(
      fpar(ftext(paste0("Table 2. Minimum and maximum average annual values for selected climate variables from within ",SNameF,"."), fp_Fig)))

    FIG_6.1 <- block_list(
  fpar(ftext(paste0("Figure 8. Mean annual climate of ",SNameF,
                    " with area average shown in heading of each plot."), fp_Fig)))

#Climate Variable figures 
RF2file <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climate_less_RF.png")
RF2img <- external_img(src = RF2file, height = 3,width = 3) 
TAAfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climate_less_TA.png")
TAAimg <- external_img(src = TAAfile, height = 3,width = 3) 
RHfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climate_less_RH.png")
RHimg <- external_img(src = RHfile, height = 3,width = 3) 
SRfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climate_less_SR.png")
SRimg <- external_img(src = SRfile, height = 3,width = 3) 
SMfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climate_less_SM.png")
SMimg <- external_img(src = SMfile, height = 3,width = 3) 
ETfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climate_less_ET.png")
ETimg <- external_img(src = ETfile, height = 3,width = 3) 
WSfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climate_less_WS.png")
WSimg <- external_img(src = WSfile, height = 3,width = 3)

si<-1.65
###############   Slide 8 Average monthly rainfall 
p<-p+1
p11<-p

S7_TIT<- block_list(
  fpar(ftext("Average Monthly Rainfall", FTXTT),fp_p = fp_par(text.align = "center")))
  # fpar(ftext("and Temperature", FTXTT),fp_p = fp_par(text.align = "center")))

CLIMA <- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Annual Climate.csv"),sep=",")

RFT <- CLIMA[1,2:13]
MinT<- CLIMA[2,2:13]
TAT <- CLIMA[3,2:13] # Mean Temperature 
MaxT <- CLIMA[4,2:13]
RH <- CLIMA[5,2:13]

MONLIST <- c("January","February","March","April","May","June","July","August","September","October","November","December")

RFx <- max(RFT)
RFxm <- match(RFx,RFT)
RFn <- min(RFT)
RFnm <- match(RFn,RFT)
RFx<-round(RFx, 0)
RFn<-round(RFn, 0)

CLIMO_RF <- paste0("Average monthly rainfall patterns vary over the course of the year. At ", SNameF, ", the highest monthly rainfall is received in ", MONLIST[RFxm], 
                " (",RFx, RFUnit,".) and the lowest monthly rainfall is received in ", MONLIST[RFnm]," (",RFn, RFUnit,".).")

fp_CLIMO_RF <- fpar(ftext(CLIMO_RF, fp_Tx))

FIG_6a <- block_list(
  fpar(ftext(paste0("Figure 9. Mean monthly rainfall at ",SNameF,"."), fp_Fig)))
  
FIG_7a <- block_list(
  fpar(ftext(paste0("Figure 10. Monthly rainfall maps for the wettest (top) and driest (bottom) months."), fp_Fig)))

#Monthly rainfall barchart
RFBfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climograph_RF.png")
RFBimg <- external_img(src = RFBfile, height = 2,width = 1) 

#Month rainfall maps
RFWfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," RF_wm.png")
RFWimg <- external_img(src = RFWfile, height = 3,width = 3) 

RFDfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," RF_dm.png")
RFDimg <- external_img(src = RFDfile, height = 3,width = 3) 

###############   Slide 9 Average monthly air temperature 
p<-p+1
p12<-p

S8a_TIT<- block_list(
  fpar(ftext("Average Monthly Temperature", FTXTT),fp_p = fp_par(text.align = "center")))

RFT <- CLIMA[1,2:13]
MinT<- CLIMA[2,2:13]
TAT <- CLIMA[3,2:13] # Mean Temperature 
MaxT <- CLIMA[4,2:13]
RH <- CLIMA[5,2:13]

MONLIST <- c("January","February","March","April","May","June","July","August","September","October","November","December")
TAx <- max(TAT)
TAxm <- match(TAx,TAT)
TAn <- min(TAT)
TAnm <- match(TAn,TAT)
Tdif <- abs(TAn-TAx)
TAx<-round(TAx, 0)
TAn<-round(TAn, 0)

# Get annual average variation in air temp
min.col <- function(m, ...) max.col(-m, ...)

TAh<-TAT[max.col(TAT)]
TAl<-TAT[min.col(TAT)]
TAv<-round(TAh-TAl, 2)
TAv

CLIMO_TA <- paste0("Average monthly air temperature patterns vary over the course of the year. At ", SNameF, " there is a ",TAv, TUnit2," annual variation in temperature,",
                " with the warmest month of ",MONLIST[TAxm]," (", TAx, TUnit2,") and the coolest month in ",MONLIST[TAnm],
                " (", TAn, TUnit2,").")

fp_CLIMO_TA <- fpar(ftext(CLIMO_TA, fp_Tx))

FIG_8a <- block_list(
  fpar(ftext(paste0("Figure 11. Mean monthly air temperature at ",SNameF," with monthly rainfall in the background."), fp_Fig)))

FIG_9a <- block_list(
  fpar(ftext(paste0("Figure 12. Monthly air temperature maps for the coldest (top) and hottest (bottom) months."), fp_Fig)))

#Monthly temp line chart
TALfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climograph_AT.png")
TALimg <- external_img(src = TALfile, height = 2,width = 1) 

#Month rainfall maps
TACfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," TA_cm.png")
TACimg <- external_img(src = TACfile, height = 3,width = 3) 

TAHfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," TA_hm.png")
TAHimg <- external_img(src = TAHfile, height = 3,width = 3) 
###############   Slide 8 TA And RF
p<-p+1
p13<-p

S8_TIT<- block_list(
  fpar(ftext("Average Monthly Climate", FTXTT),fp_p = fp_par(text.align = "center")))

TAfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," TA12.png")
TAimg <- external_img(src = TAfile , height = 3,width = 3) 

RFFfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," RF12.png")
RFFimg <- external_img(src = RFFfile , height = 3,width = 3) 

FIG_5 <- block_list(
  fpar(ftext(paste0("Figure 13. Mean monthly rainfall  ",SNameF,
                    " with area average shown in heading of each plot."), fp_Fig)))

FIG_6 <- block_list(
  fpar(ftext(paste0("Figure 14. Mean monthly temperature at ",SNameF,
                    " with area average shown in heading of each plot."), fp_Fig)))

###############   Slide 9 Seasonal Rainfall 
p<-p+1
p14<-p

S9_TIT<- block_list(
  fpar(ftext("Average Seasonal Rainfall", FTXTT),fp_p = fp_par(text.align = "center")))

CLIM

#Season 
DSeaMRF  <- round(CLIM[1,13],0)
DSeaMRFx <- round(CLIM[2,13],0)
DSeaMRFn <- round(CLIM[3,13],0)
WSeaMRF  <- round(CLIM[1,14],0)
WSeaMRFx <- round(CLIM[2,14],0)
WSeaMRFn <- round(CLIM[3,14],0)

#per month 
DSeaMRF_M<- round((CLIM[1,13]/6),1)
DSeaMRF_Mx<- round((CLIM[2,13]/6),1)
DSeaMRF_Mn<- round((CLIM[3,13]/6),1)
WSeaMRF_M <- round((CLIM[1,14]/6),1)
WSeaMRF_Mx <- round((CLIM[2,14]/6),1)
WSeaMRF_Mn <- round((CLIM[3,14]/6),1)

#Per Month Range
D_RF_R_M  <- paste0(DSeaMRF_Mn, " to ", DSeaMRF_Mx) 
W_RF_R_M  <- paste0(WSeaMRF_Mn, " to ", WSeaMRF_Mx)

#Season Range 
D_RF_R  <- paste0(DSeaMRFn, " to ", DSeaMRFx) 
W_RF_R  <- paste0(WSeaMRFn, " to ", WSeaMRFx) 

WetM <- WSeaMRF_M
DryM <- DSeaMRF_M
  
SeaD <- round((WetM - DryM),1)

#Season monthly rainfall percentiles
P<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF,"RF percentiles.csv"),sep=",")
D_RF_P<-P[1,2]
W_RF_P<-P[1,1]
RFUnit3
           
SEA <- paste0("Hawaii has two distinct 6-month seasons of rainfall: hooilo (Wet season: November to April) and kau (Dry season: May to October).",
              " Average Wet season monthly rainfall across ",SNameF," is ",WetM," ",RFUnit2," and Dry season is ",DryM," ",RFUnit2,
              ". These monthly values are in the ",W_RF_P," and ",D_RF_P," percentiles for rainfall across the whole state, respectively.
              
Management plans should anticipate and minimize negative impacts of these seasonal rainfall variations.")

fp_Tx19 <- fp_text(italic = TRUE, color = "black", font.size = 19)
fp_SEA <- fpar(ftext(SEA, fp_Tx19))

FIG_12.1 <- block_list(
  fpar(ftext(paste0("Figure 15. Average monthly rainfall maps for the wet (top) and dry (bottom) seasons. ",
  SNameF,"."), fp_Fig)))

#Mean CLIM Figure 
SEAfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," SeaMRF.png")
SEAimg <- external_img(src = SEAfile, height = 3,width = 3) 



###############   Slide 10  Average Climate Table
p<-p+1
p15<-p

ACLIM<- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Annual Climate.csv"),sep=",")
ACLIM2 <- ACLIM[,1:13]
ACLIM_T <- flextable(ACLIM)
ACLIM_T <- bg(ACLIM_T, bg = "yellow", part = "header")
ACLIM_T <- fontsize(ACLIM_T, size = 10)
#Creates a table for the PPT
ACLIM_T <- autofit(ACLIM_T)


S10_TIT<- block_list(
  fpar(ftext("Average Monthly Climate Table", FTXTT),fp_p = fp_par(text.align = "center")),
  fpar(ftext(SNameF, FTXTT3),                      fp_p = fp_par(text.align = "center")))

 
TAB1 <- block_list(
  fpar(ftext(paste0("Table 3. Average monthly climate variables characteristics at ",NameF, ". Where, RF is rainfall; Min TA is average minimum air temperature",
                    " Mean TA is average air temperature; Max TA is average maximum air temperature; RH is relative humidity; CF is cloud frequency; ET is ",
                    "evapotranspiration; SM is soil moisture; S is shortwave downward radiation: ANN, is annual total for rainfall and annual average for all other variables."), fp_Fig)))


################ Slide 12
p<-p+1
p16<-p

  P2_TIT<- block_list(
  fpar(ftext("Part 3: Climate Variability", FTXTT),fp_p = fp_par(text.align = "center")))
  # fpar(ftext(SNameF, FTXTT3),fp_p = fp_par(text.align = "center")))

  fp_Txa <- fp_text(italic = TRUE, color = "black", font.size = 17)
  
  InterAn <-  block_list(
    fpar(ftext(paste0("Rainfall and temperature in Hawaii can vary greatly from year-to-year due to natural climatic systems such as the El Niño-Southern Oscillation (ENSO). ",
                      "ENSO is a periodic fluctuation of ocean temperatures in the tropical Pacific, and this has a strong influence on rainfall variability. ",
                      "ENSO consists of five phases, as shown in the graph below."),fp_Txa)))

fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 17)
# fp_InterAn <- fpar(ftext(InterAn, fp_Tx))

MEIfile<-paste0(I_FOLDER,"ENSO_timeseries.png")
MEIimg<-external_img(src=MEIfile, height = 1, width = 2)

ENSOfile <- paste0(I_FOLDER,"ENSO_New.png")
ENSO2img <- external_img(src = ENSOfile, height = 1.2,width = 1.2) 

FIG_10.1 <- block_list(
  fpar(ftext(paste0("Figure 16. Timeseries of changes in sea surface temperature (SST) and associated ENSO phase from 1950 - 2022. ",SNameF,"."), fp_Fig)))

ft4 <- block_list(fpar(ftext(paste0("https://origin.cpc.ncep.noaa.gov/"), fp_ft)))

################ Slide 13 Seasonal Rainfall and ENSO
p<-p+1
p17<-p

MEIRF<- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," MEI_S.csv"),sep=",")
MEIRF

# get last year
rf<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Annual_RF_in.csv"),sep=",")
yr<-rf[nrow(rf),]$Date

EN_TIT<- block_list(
  fpar(ftext("Seasonal Rainfall and ENSO", FTXTT),fp_p = fp_par(text.align = "center")))

# SEL wet season rainfall
SEL_RFW <- MEIRF[which(MEIRF$Season == "wet"),]
SEL_RFW <- SEL_RFW[which(SEL_RFW$Phase == "Strong El Nino"),]$Mean
SEL_RFW

# SLA dry season rainfall
SLA_RFD <- MEIRF[which(MEIRF$Season == "dry"),]
SLA_RFD <- SLA_RFD[which(SLA_RFD$Phase == "Strong La Nina"),]$Mean
SLA_RFD

# Longterm average seasonal rainfalls (across all ENSO phases)
WetSA <- MEIRF[which(MEIRF$Season == "wet"),]
WetSA <- mean(WetSA$Mean)
WetSA

DrySA <- MEIRF[which(MEIRF$Season == "dry"),]
DrySA <- mean(DrySA$Mean)
DrySA

# SEL wet season rainfall % difference from overall average 
DifWRF <- round(((WetSA - SEL_RFW) / WetSA) * 100,0)
DifWRF

# SLA dry season rainfall % difference from overall average
DifDRF <- round(((DrySA - SLA_RFD) / DrySA) * 100,0)
DifDRF

MEI2 <- paste0("In Hawaii, the Warm (El Niño) phase typically brings below average rainfall during the wet season, and above average rainfall in the dry season. This pattern is reversed for the Cool (La Niña) phase.

At ",SNameF,", the wet season during a Strong El Niño is ",DifWRF,"% dryer than the long-term wet season average, and the dry season during a Strong La Niña is ",DifDRF,"% dryer than average. These patterns influence drought conditions and wildfire susceptibility, and management activities can benefit from incorporating this ENSO-influenced seasonal rainfall variability.")
MEI3 <-  block_list(
  fpar(ftext(MEI2, fp_Txa)))

MEISfile <- paste0(R_FOLDER,NameF,"/",NameF,"ENSO_season_barplot.png")
MEISimg <- external_img(src = MEISfile, height = 1,width = 2.5) 

FIG_14.1 <- block_list(
  fpar(ftext(paste0("Figure 17. Barplot of average monthly rainfall grouped by season and ENSO phase. Numbers above the bars are how many seasons from 1950 to ",yr," fell within each ENSO phase."), fp_Fig)))

#################### Slide 14 ###############################
p<-p+1
p18<-p

S13_TIT<- block_list(
  fpar(ftext("Long-Term Trends in Rainfall", FTXTT),fp_p = fp_par(text.align = "center")))

Trend <- paste0("Linear trends in annual and seasonal rainfall at ",SNameF," have been calculated over three different periods since 1920 ",
                "to show long, mid, and short-term trends. The directions of change for the annual plot are shown below.")
                
fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 20) 
fp_Trend <- fpar(ftext(Trend, fp_Tx))

Trendfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," RF_Trend.png")
Trendimg <- external_img(src = Trendfile, height = 2,width = 2) 

### Bring in table of periods and trend directions
  RFT<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," RF_trend_dirctions.csv"),sep=",")
  RFT_T <- flextable(RFT)
  # RFT_T <- bold(RFT_T, bold = T, part = "header")
  RFT_T <- fontsize(RFT_T, size = 14)
  #Creates a table for the PPT
  RFT_T <- autofit(RFT_T)

 TAB2.1 <- block_list(
    fpar(ftext(paste0("Table 4. Direction of trendline for annual average monthly rainfall over three periods within the record."), fp_Fig2)))
  
FIG_9 <- block_list(
  fpar(ftext(paste0("Figure 18. Rainfall time series (1920-",yr,") at " ,SNameF," with ",
                    "linear trends. Trendlines with p-value < 0.05 are statistically significant."), fp_Fig2)))


#################### Slide 14 (NEW) Air Temperature Monthly and Annual #########################
p<-p+1
p19<-p

S14_TIT<- block_list(
  fpar(ftext("Trends in Air Temperature", FTXTT),fp_p = fp_par(text.align = "center")))

# use monthly air temperature data to figure out values for text
AT<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF,"_monthly_airtemp.csv"),sep=",")

# get average annual air temperature values
AT2<-aggregate(mean ~ year, AT, mean)

# find start and end values for linear trend
s<-AT2[1,]$year
e<-AT2[nrow(AT2),]$year
yrs<-e-s

mod<-lm(mean ~ year, data=AT2)
d<-data.frame(year=c(s,e))
pred<-predict(mod,d)

# determine direction of change
if(pred[[1]] > pred[[2]]) {airchange <- c("decreased")}
if(pred[[1]] < pred[[2]]) {airchange <- c("increased")}

# # calculate % change
# airchange2<-paste0(abs(round((((pred[[2]] - pred[[1]])/pred[[1]])*100), 1)),"%.")

# calculate absolute change
airchange2<-paste0(abs(round(pred[[2]] - pred[[1]])),"°F.")

# calculate average annual range (hottest - coldest month)
AT$range<-AT$max - AT$min
ra<-round(mean(AT$range), 1)

# find month-year with hottest temp
ATd<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF,"_daily_airtemp.csv"),sep=",")
head(ATd)
mt<-max(ATd$max.x)
mt2<-round(mt, 1)

hot<-ATd[which(ATd$max.x == mt),]
if(nrow(hot)>1) {hot<-hot[nrow(hot),]}
hot.y<-hot$year

if(hot$month == 1) {hot.m<-c("January")}
if(hot$month == 2) {hot.m<-c("February")}
if(hot$month == 3) {hot.m<-c("March")}
if(hot$month == 4) {hot.m<-c("April")}
if(hot$month == 5) {hot.m<-c("May")}
if(hot$month == 6) {hot.m<-c("June")}
if(hot$month == 7) {hot.m<-c("July")}
if(hot$month == 8) {hot.m<-c("August")}
if(hot$month == 9) {hot.m<-c("September")}
if(hot$month == 10) {hot.m<-c("October")}
if(hot$month == 11) {hot.m<-c("November")}
if(hot$month == 12) {hot.m<-c("December")}

AirTrend <- paste0("Trends in air temperature have been calculated over a ",yrs,"-year record at ",SNameF,
". From ",s," to ",e," average annual air temperature has ",airchange, " by ", airchange2)
AirTrend

fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 20) 
fp_AirTrend <- fpar(ftext(AirTrend, fp_Tx))

AirTrend2 <- paste0("At this site there is an average range of ",ra,"°F between the hottest and coldest months ",
                    "within a single year. The highest monthly temperature of ",mt2,"°F was recorded in ",hot.m," ",hot.y,".") 
fp_AirTrend2 <- fpar(ftext(AirTrend2, fp_Tx))

Trendfile1 <- paste0(R_FOLDER,"/",NameF,"/",NameF,"_monthly_airtemp.png")
Trendimg1 <- external_img(src = Trendfile1, height = 1,width = 2)

Trendfile2 <- paste0(I_FOLDER,"airtemp_legend2.jpg")
Trendimg2 <- external_img(src = Trendfile2, height = 1,width = 1)

FIG_10_11 <- block_list(
  fpar(ftext(paste0("Figure 19. ",yrs,"-year (",s," to ",e,") monthly air temperature time series at ",SNameF,
                    ". The linear trend is determined to be statistically significant when the ",
                    "p-value is less than 0.05."), fp_Fig2)))

ft5 <- block_list(fpar(ftext(paste0("https://www.hawaii.edu/climate-data-portal/"), fp_ft)))

################ Slide 15 PART 4 Drought and Fire Occurrence ###############
#Today, 
#increased wildland fire activity, 
#and damage to terrestrial and aquatic habitats - all of which can contribute to substantial economic losses# 
p<-p+1
p20<-p

P3_TIT<- block_list(
  fpar(ftext("Part 4: Drought and Fire History", FTXTT)),
  fpar(ftext(SNameF, FTXTT3),fp_p = fp_par(text.align = "center")))

Part3 <- paste0("Drought is a prominent feature of the climate system in Hawaii and can cause severe impacts across multiple sectors. ",
                  "Droughts in Hawaii often result in reduced crop yields, loss of livestock, drying of streams and reservoirs, depletion ",
                   "of groundwater, and increased wildland fire activity. These impacts can cause substantial economic ",
                  "losses as well as long-term damage to terrestrial and aquatic habitats.")

fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 20) 
fp_Part3 <- fpar(ftext(Part3, fp_Tx))



####################  SLIDE 16 FIVE TYPES OF DROUght######################################
p<-p+1
p21<-p

TIT_D <- fpar(ftext("Five Types of Drought", FTXTT),fp_p = fp_par(text.align = "center") )
# fp_D_M <-  block_list(
#   fpar(ftext( "Meteorological Drought.", fp_Fig4)),
#   fpar(ftext( "Defined as a lack of rainfall" , fp_Fig3)))
# 
# fp_D_A <-  block_list(
#   fpar(ftext( "Agriculture Drought", fp_Fig4)),
#   fpar(ftext( "Refers to a period of declining soil moisture and subsequent crop failure", fp_Fig3)))
# 
# fp_D_H <-  block_list( 
#   fpar(ftext( "Hydrological Drought", fp_Fig4)),
#   fpar(ftext( "Expressed as decreased streamflow and sub-surface water storage", fp_Fig3)))
# 
# fp_D_E<-  block_list(
#   fpar(ftext( "Ecological Drought", fp_Fig4)),
#   fpar(ftext( "Includes any impacts to ecosystems including an increase in wildfire occurrence" , fp_Fig3)))
# 
# fp_D_S<-  block_list(   
#   fpar(ftext( "Socio-Economic Drought", fp_Fig4)),
#   fpar(ftext( "Includes impacts to social and economic systems, including increased costs or revenue losses, and impacts on public health and safety", fp_Fig3)))

D_t <- paste0("There are five major types of drought. During droughts there is sometimes a progression from ",
              "one type to the next. Depending on local factors, these drought types can also happen simultaneously ",
              "and at different levels of severity.")

fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 18) 
fp_D_t <- fpar(ftext(D_t, fp_Tx))

#################### Slide 17 ###############################
p<-p+1
p22<-p

S15_TIT<- block_list(
  fpar(ftext("Identifying Droughts Using the", FTXTT),fp_p = fp_par(text.align = "center")),
  fpar(ftext("Standard Precipitation Index", FTXTT),fp_p = fp_par(text.align = "center")))


SPI <- paste0("The Standardized Precipitation Index (SPI) is one of the most widely used indices for meteorological drought. SPI compares current rainfall with ",
              "its multi-year average, so that droughts are defined relative to local average rainfall. This standardized index allows wet and dry climates ",
              "to be represented on a common scale. Here, 100+ years of monthly rainfall are used to used to calculate SPI-12, which compares how ",
              "a 12-month period compares with all 12-month periods in the record. SPI-12 is a good measure of sustained droughts that affect hydrological processes ",
              "at ",SNameF, ".")

  fp_SPI <- fpar(ftext(SPI, fp_Txa))

  SPIfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," SPI.png")
  SPIimg <- external_img(src = SPIfile, height = 2,width = 2) 

  FIG_12 <- block_list(
   fpar(ftext(paste0("Figure 20. SPI-12 time series (1920-",yr,") at ", 
                      SNameF ,". Positive SPI (blue) indicate wet periods, negative ", 
                     "SPI (red) indicate dry periods."), fp_Fig)))


  
  #################### Slide 18 ###############################
  p<-p+1
  p23<-p
  
  S16_TIT<- block_list(
    fpar(ftext("A 100+ Year History of Drought", FTXTT),fp_p = fp_par(text.align = "center")))
   

    #INFO for Slides 16 & 17
  DROT <- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Drought History.csv"),sep=",")

    D_Sum <- sum(!is.na(DROT$A.Intensity))
  #if(DROT$P.Intensity[D_Sum] > 1){DROT <- DROT[-D_Sum,]}
  EX_Cnt <- sum(DROT$P.Intensity > 2)
  SV_Cnt <- sum(DROT$P.Intensity > 1.5)
  MO_Cnt <- sum(DROT$P.Intensity > 1)
  CNTDRT <- sum(!is.na(DROT$A.Intensity))
  MO_D <- abs(MO_Cnt - SV_Cnt)
  SV_D <- abs(SV_Cnt - EX_Cnt)
  EX_D <- EX_Cnt
  SoG <- SV_D + EX_D
  DROT$A.Intensity <- abs(DROT$A.Intensity)
  DROT$P.Intensity <- abs(DROT$P.Intensity)
  DROT$Magnitude <- abs(DROT$Magnitude)
  LongD <- max(DROT$Duration)
  colnames(DROT) <- c("Start Date", "End Date","Duration (months)","Average Intensity","Peak Intensity","Magnitude")

  # count rows in DROT and separate into 12-row tables
  Drot_ct<-nrow(DROT)
  
  # split table up into smaller pieces
  if(Drot_ct <=12) {DROT1<-DROT}
  if(Drot_ct >12) {DROT1<-DROT[1:12,]}
  if(Drot_ct > 12) {DROT2<-DROT[13:nrow(DROT),]}
  if(Drot_ct > 24) {DROT2<-DROT2[1:12,]}
  if(Drot_ct > 24) {DROT3<-DROT[25:nrow(DROT),]}

  DrotT_T1 <- flextable(DROT1[1:6])
  DrotT_T1 <- bg(DrotT_T1, bg = "orange", part = "header")
  DrotT_T1 <- autofit(DrotT_T1)
  DrotT_T1

  if(exists("DROT2")) {DrotT_T2 <- flextable(DROT2[1:6]);
                       DrotT_T2 <- bg(DrotT_T2, bg = "orange", part = "header");
                       DrotT_T2 <- autofit(DrotT_T2)}
  
  if(exists("DROT3")) {DrotT_T3 <- flextable(DROT3[1:6]);
                       DrotT_T3 <- bg(DrotT_T3, bg = "orange", part = "header");
                       DrotT_T3 <- autofit(DrotT_T3)}

  # # make conditional statement to get last slide number (will be used on last slide)
  # if(Drot_ct>12) {ls <- 33}
  # if(Drot_ct>24) {ls <- 34}
  
  DHist1 <- paste0("Negative SPI values (dry periods) are inverted to show a complete drought timeseries at ", SNameF, ". Dashed lines and corresponding color coding indicates instances ",
                 "of Moderate (SPI > 1), Severe (SPI > 1.5), and Extreme (SPI > 2) drought.") 
                 
  DHist2 <- paste0("A total of ",CNTDRT," droughts were observed over the entire record ",
                   "with a total of ", SoG, " drought events of severe strength or greater. The longest drought lasted for a total of ",LongD," months (see Annex III).")
 
  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 18)
  
  fp_DHist<- block_list(
             fpar(ftext(DHist1, fp_Tx)),
             fpar(ftext(" ", fp_Tx)),
             fpar(ftext(DHist2, fp_Tx)))
  
  DHistfile <- paste0(R_FOLDER,"/",NameF,"/",NameF,"Drought_History.png")
  DHistimg <- external_img(src = DHistfile, height = 2,width = 4) 

  FIG_13 <- block_list(
    fpar(ftext(paste0("Figure 21. SPI time series (1920-",yr,") (reversed axis) at ", SNameF ,
                       ". Dashed lines show, moderate (yellow), severe (red), and extreme (dark red), drought thresholds."), fp_Fig)))
  
  
  #################### Slide 20 ############################### 
  p<-p+1
  p24<-p
  
  SPI_NA<- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," SPI_NEGS_ALL.csv"),sep=",")
  SPI_NA<-SPI_NA[which(!is.na(SPI_NA$date)),]
  SPICNT<- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Drought Count.csv"),sep=",")
  
  SP_TIT<- block_list(
    fpar(ftext("Short-term vs Long-term Droughts", FTXTT),fp_p = fp_par(text.align = "center")))
  
  SPI_3_SMfile <- paste0(R_FOLDER,"/",NameF,"/",NameF,"Drought_HistoryS_3.png")
  SPI3_SMimg <- external_img(src = SPI_3_SMfile, height = 2,width = 4)
  
  SPI_12_SMfile <- paste0(R_FOLDER,"/",NameF,"/",NameF,"Drought_HistoryS_12.png")
  SPI12_SMimg <- external_img(src = SPI_12_SMfile, height = 2,width = 4)
  
  FIG_14 <- block_list(
    fpar(ftext(paste0("Figure 22. SPI-3 time series (reversed axis) at ", NameF ,"."), fp_Fig)))
  FIG_15 <- block_list(
    fpar(ftext(paste0("Figure 23. SPI-12 time series (reversed axis) at ", NameF ,"."), fp_Fig)))
  SPI3v12<- block_list(
    fpar(ftext(paste0("SPI-3 provides a comparison of rainfall over a specific 3-month period and reflects short-term ", 
                       "conditions. SPI-12 provides 12-month comparisions and reflects long-term conditions. It is important ",
                       "to consider both timescales for planning."), fp_Txa)))
  
  # bring in monthly rainfall dataset to get dates
  MRF<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Monthly Rainfall_in.csv"),sep=",")
  RFC<-MRF[nrow(MRF),]
  
  RFC.y<-RFC$Year
  
  if(RFC$Month == 1) {RFC.m<-c("January")}
  if(RFC$Month == 2) {RFC.m<-c("February")}
  if(RFC$Month == 3) {RFC.m<-c("March")}
  if(RFC$Month == 4) {RFC.m<-c("April")}
  if(RFC$Month == 5) {RFC.m<-c("May")}
  if(RFC$Month == 6) {RFC.m<-c("June")}
  if(RFC$Month == 7) {RFC.m<-c("July")}
  if(RFC$Month == 8) {RFC.m<-c("August")}
  if(RFC$Month == 9) {RFC.m<-c("September")}
  if(RFC$Month == 10) {RFC.m<-c("October")}
  if(RFC$Month == 11) {RFC.m<-c("November")}
  if(RFC$Month == 12) {RFC.m<-c("December")}
  
  ### determine current drought status
  # SPI3
  # spi3A<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF,"SPI_3.csv"),sep=",")
  head(SPI_NA, 20)
  spi3A<-SPI_NA[which(SPI_NA$m.scale == 3),]
  spi3C<-round(spi3A[nrow(spi3A),]$spi_negs,1)
  spi3C
  tail(spi3A,30)
  
  spi3H<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF,"Drought History SPI_3.csv"),sep=",")
  spi3HC<-spi3H[nrow(spi3H),]
  spi3HC
  
  if(spi3HC$P.Intensity <= 1.5) {spi3HRL <- "moderate"}
  if(spi3HC$P.Intensity <= 2 && spi3HC$P.Intensity > 1.5) {spi3HRL <- "severe"}
  if(spi3HC$P.Intensity >2) {spi3HRL <- "extreme"}
  
  spi3HRs<-spi3HC$Start
  spi3HRe<-spi3HC$End
  
  if(!is.na(spi3HC$End)) {spi3CT <- paste0("SPI-3: Currently not in drought. Most recently there was ",
                                           spi3HRL," drought from ",spi3HRs," - ",spi3HRe,". Current drought intensity is ",
                                           spi3C,".")}
  if(is.na(spi3HC$End)) {spi3CT <- paste0("SPI-3: Currently in ",spi3HRL, " drought since ",spi3HRs,". ",
                                          "Current drought intensity is ",spi3C,".")}
spi3CT

  # SPI12
  # spi12A<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF,"SPI_12.csv"),sep=",")
  spi12A<-SPI_NA[which(SPI_NA$m.scale == 12),]
  spi12C<-round(spi12A[nrow(spi12A),]$spi_negs, 1)
  tail(spi12A,20)
  
  spi12H<-read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Drought History.csv"),sep=",")
  spi12HC<-spi12H[nrow(spi12H),]
  spi12HC
  
  if(spi12HC$P.Intensity <= 1.5) {spi12HRL <- "moderate"}
  if(spi12HC$P.Intensity <= 2 && spi12HC$P.Intensity > 1.5) {spi12HRL <- "severe"}
  if(spi12HC$P.Intensity >2) {spi12HRL <- "extreme"}
  
  spi12HRs<-spi12HC$Start2
  spi12HRe<-spi12HC$End2
  
  if(!is.na(spi12HC$End)) {spi12CT <- paste0("SPI-12: Currently not in drought. Most recently there was ",
                                           spi12HRL," drought from ",spi12HRs," - ",spi12HRe,". Current drought intensity is ",
                                           spi12C,".")}
  if(is.na(spi12HC$End)) {spi12CT <- paste0("SPI-12: Currently in ",spi12HRL, " drought since ",spi12HRs,". ",
                                          "Current drought intensity is ",spi12C,".")}
  spi12CT
  
  SPI3v12.1<- block_list(
    fpar(ftext(paste0("As of ",RFC.m," ",RFC.y," the most recent drought events are as follows:"), fp_NM8)))
  
  SPI3v12.2<- block_list(
    fpar(ftext(paste0(spi3CT), fp_Txa)))
  
  SPI3v12.3<- block_list(
    fpar(ftext(paste0(spi12CT), fp_Txa)))
  
                       # A total of ", SPICNT[1,3], " droughts were observed at ",
                       #  SNameF , ". Over this same time period, only ", SPICNT[1,4] ," droughts were identified when looking at the SPI-12 timeseries. It is important to compare the 3-month SPI with longer time scales. ",
                       # "A relatively normal 3-month period could occur in the middle of a longer-term drought that would only be visible at longer time scales. Looking at longer time scales ",
                       # "would prevent a misinterpretation that a drought might be over."), fp_Txa)))

  #################### Slide 21 ###############################
  p<-p+1
  p25<-p
  
  S18_TIT<- block_list(
    fpar(ftext(paste("Fire Occurrence in ",ISL), FTXTT),fp_p = fp_par(text.align = "center")))

  
  FIRE <- paste0("Ecological drought often drives an increase in wildfire occurrence. In Hawaii, wildfires are most extensive in dry ",
                 "and mesic non-native grass and shrublands. During drought events, wildfire risk in these areas increases rapidly. ",
                 "Currently, agricultural abandonment is resulting in increased grass and shrublands. ",
                 "This combined with recurring incidences of drought is expected to increase the risk of future wildfire in Hawaii.")
                  
  
  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 20) 
  fp_Fire <- fpar(ftext(FIRE, fp_Tx))
  
  FireMfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Fire.png")
  FireMimg <- external_img(src = FireMfile) 
  
  
  FIG_16 <- block_list(
    fpar(ftext(paste0 ("Figure 24. The map shows wildfires that have occurred on the island of ", ISL, " between 1999 and 2020.") , fp_Fig)))
  
  
  ################ Slide 22 PART 5 Downscaling ###############
  p<-p+1
  p26<-p
  
  P4_TIT<- block_list(
    fpar(ftext("Part 5: Future Climate", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext(SNameF, FTXTT3),fp_p = fp_par(text.align = "center")))
  
  Part4<- paste0("Global Climate Models are used to predict future changes in rainfall and temperature, simulating future conditions under different ",
                 "scenarios for how much carbon dioxide we emit into the air. Two common scenarios are RCP 4.5 which assumes we reduce our carbon ",
                 "emissions, and RCP 8.5, which is an increased emissions scenario.")
  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 18) 
  fp_Part4 <- fpar(ftext(Part4, fp_Tx))
  
  Part4.2<- paste0("Data downscaling is used make these models useful at the local management level. In Hawaii, two types of downscaled ",
                   "projections are available:")
  fp_Part4.2 <- fpar(ftext(Part4.2, fp_Tx))
  
  Part4.3<- paste0("Statistical: Available for Mid & End-of-Century Dynamical: Only available for End-of-Century")
  fp_Tx4.3 <- fp_text(italic = TRUE, bold = TRUE, color = "black", font.size = 18) 
  fp_Part4.3 <- fpar(ftext(Part4.3, fp_Tx4.3))
  
  Part4.4<- paste0("Both downscaling projections are presented here. These two projections sometimes agree with each other, ",
                   "and other times they provide conflicting results. When viewing the maps, we can observe where these similarities ",
                   "and differences are, for example which areas show reduced rainfall under both projections.")
  fp_Part4.4 <- fpar(ftext(Part4.4, fp_Tx))

  # M_Down <-  block_list(
  #   fpar(ftext("In Hawaii, two types of downscaled " , fp_NM3)),
  #   fpar(ftext("projections are available.", fp_NM3)),
  #   fpar(ftext("Dynamical Downscaling (End of Century)" , fp_NM3)),
  #   fpar(ftext("Statistical Downscaling (Mid & End of Century)" , fp_NM3)),
  #   fpar(ftext("Results for both types of downscaling ", fp_NM3)),
  #   fpar(ftext("and both scenarios will be shown here.", fp_NM3)))
   
  
  #################### Slide 23 ###############################
  p<-p+1
  p27<-p
  
  S21_TIT<- block_list(
    fpar(ftext("Average Rainfall Change", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext("Mid-Century (2040-2070)", FTXTT),fp_p = fp_par(text.align = "center")))
  
  Down <- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Downscaling.csv"),sep=",")

  # Get PERCENT CHANGE values
  RFA_Thresh40 <-   Down[c(15,18),2]
  RFA_Thresh40_4.5 <-   Down[c(15),2]
  RFA_Thresh40_8.5 <-   Down[c(18),2]

  # RFD_Thresh40 <-   Down[c(16,19),2]
  # RFW_Thresh40 <-   Down[c(17,20),2]
  
  # Get AMMOUNT (inches) change values
  RFA <- CLIMA[1,14]
  RF40_4.5 <- round(RFA * (RFA_Thresh40_4.5*0.01),0)
  RF40_8.5 <- round(RFA * (RFA_Thresh40_8.5*0.01),0)

  AnnualRM_in  <- paste0(RF40_4.5," to ",RF40_8.5," in/year")
  AnnualRM_in
  AnnualRM  <- paste0(RFA_Thresh40_4.5," to ",RFA_Thresh40_8.5,"% change from present")
  AnnualRM
  
  # choose whether text says "increase" or "decrease"
  ifelse(RFA_Thresh40_4.5>=0, RF40_4.5_ch<-"increase", RF40_4.5_ch<-"decrease")
  ifelse(RFA_Thresh40_8.5>=0, RF40_8.5_ch<-"increase", RF40_8.5_ch<-"decrease")
  
  RF_Down4 <-  block_list(
    fpar(ftext("Rainfall for Years 2040-2070" , fp_NM5)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext("Change in Annual Rainfall", fp_NM6)),
    fpar(ftext(AnnualRM_in, fp_NM3)),
    fpar(ftext(AnnualRM,fp_NM3)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext(paste0("These Statistical Downscaling maps show the projected change in rainfall under RCP 4.5 and 8.5 ",
                      "conditions. At ",SNameF,", annual rainfall is projected to ",RF40_4.5_ch," by ",abs(RF40_4.5)," inches (RCP 4.5), ",
                      "or ",RF40_8.5_ch," by ",abs(RF40_8.5)," inches (RCP 8.5) by mid-century."), fp_NM2)))
  RF_Down4
  
  #RF40file <- paste0(R_FOLDER,"/",NameF,"/",NameF," StDsRF2040.png")
  #RF40img <- external_img(src = RF40file, height = 6,width = 4) 
  
  RF40file <- paste0(R_FOLDER,"/",NameF,"/",NameF," StDsRF2040.png")
  RF40img <- external_img(src = RF40file) 
  
  FIG_18<- block_list(
    fpar(ftext(paste0("Figure 25. Downscaled future rainfall projections (% change from present) ",
                      "at ",SNameF," by mid-century (2040-2070) using Statistical Downscaling. RCP 4.5 ",
                      "(left) and RCP 8.5 (right)."), fp_Fig)))
  
  
  ##################### Slide 24 #################################################
  p<-p+1
  p28<-p
  
  S20_TIT<- block_list(
    fpar(ftext("Average Rainfall Change", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext("End-of-Century (2100)", FTXTT),fp_p = fp_par(text.align = "center")))

  # get PERCENT change
  RFA_Thresh100_D4.5 <- Down[1,2]
  RFA_Thresh100_D8.5 <- Down[4,2]
  RFA_Thresh100_S4.5 <- Down[9,2]
  RFA_Thresh100_S8.5 <- Down[12,2]
  
  # Find largest percent change value for each RCP scenario
  ifelse(abs(RFA_Thresh100_D4.5)>=abs(RFA_Thresh100_S4.5), RF100_4.5c<-RFA_Thresh100_D4.5, RF100_4.5c<-RFA_Thresh100_S4.5)
  ifelse(abs(RFA_Thresh100_D8.5)>=abs(RFA_Thresh100_S8.5), RF100_8.5c<-RFA_Thresh100_D8.5, RF100_8.5c<-RFA_Thresh100_S8.5)
  
  RFA_Thresh100 <-  Down[c(1,4,9,12),2]
  # RFD_Thresh100 <-  Down[c(2,5,10,13),2]
  # RFW_Thresh100 <-  Down[c(3,6,11,14),2]
  # RFA_Thresh40 <-   Down[c(15,18),2]
  # RFD_Thresh40 <-   Down[c(16,19),2]
  # RFW_Thresh40 <-   Down[c(17,20),2]
  RFA_100m<-min(RFA_Thresh100)
  # report minimum and maximum percent change
  AnnualR  <- paste(min(RFA_Thresh100), " to ", max(RFA_Thresh100), "% change from present")
  AnnualR
  
  # get minimum and maximum AMMOUNT change
  RFA <- CLIMA[1,14]
  MinRF <- round(RFA * (min(RFA_Thresh100*0.01)),0)
  MaxRF <- round(RFA * (max(RFA_Thresh100*0.01)),0)
  if(MinRF>=0) {MinRF<-paste0("+",MinRF)}
  if(MaxRF>=0) {MaxRF<-paste0("+",MaxRF)}
  
  # get RCP setting AMMOUNT changes
  RF100_4.5a <- round(RFA * (RF100_4.5c*0.01),0)
  RF100_8.5a <- round(RFA * (RF100_8.5c*0.01),0)

  # choose whether text says "increase" or "decrease"
  ifelse(RF100_4.5c>=0, RF100_4.5t<-"increase", RF100_4.5t<-"decrease")
  ifelse(RF100_8.5c>=0, RF100_8.5t<-"increase", RF100_8.5t<-"decrease")

  # choose whether change value has + or - before it
  ifelse(RF100_4.5c>=0, RF100_4.5c<-paste0("+",RF100_4.5c), RF100_4.5c<-RF100_4.5c)
  ifelse(RF100_8.5c>=0, RF100_8.5c<-paste0("+",RF100_8.5c), RF100_8.5c<-RF100_8.5c)

 # determine if the downscaling method outputs agree on direction of change, and choose the appropriate text
 ifelse(RFA_Thresh100_D4.5>=0, a<-1, a<-2)
 ifelse(RFA_Thresh100_S4.5>=0, b<-1, b<-2)
 ifelse(RFA_Thresh100_D8.5>=0, c<-1, c<-2)
 ifelse(RFA_Thresh100_S8.5>=0, d<-1, d<-2)
 
 if(a==b && c==d)
 {RF2100t<-ftext(paste0("These Dynamical and Statistical Downscaling maps show the projected ",
                        "change in rainfall under RCP 4.5 and 8.5 conditions. At ",SNameF,
                        ", annual rainfall is projected to ",RF100_4.5t," by ",abs(RF100_4.5a),
                        " inches (RCP 4.5), or ",RF100_8.5t," by ",abs(RF100_8.5a)," inches (RCP 8.5) by end-of-century."), fp_NM2)}
 if(a!=b && c==d)
 {RF2100t<-ftext(paste0("These Dynamical and Statistical Downscaling maps show the projected ",
                        "change in rainfall under RCP 4.5 and 8.5 conditions. At ",SNameF,
                        ", the models do not agree on the direction of change for RCP 4.5. Annual ",
                        "rainfall is projected to ",RF100_8.5t," by ",abs(RF100_8.5a)," inches (RCP 8.5) by end-of-century."), fp_NM2)}
 if(a==b && c!=d)
 {RF2100t<-ftext(paste0("These Dynamical and Statistical Downscaling maps show the projected ",
                        "change in rainfall under RCP 4.5 and 8.5 conditions. At ",SNameF,
                        ", annual rainfall is projected to ",RF100_4.5t," by ",abs(RF100_4.5a),
                        " inches (RCP 4.5). The models do not agree on the direction of change for RCP 8.5."), fp_NM2)}
 if(a!=b && c!=d)
 {RF2100t<-ftext(paste0("These Dynamical and Statistical Downscaling maps show the projected ",
                        "change in rainfall under RCP 4.5 and 8.5 conditions. At ",SNameF,
                        ", the models do not agree on the overall direction of change for either RCP 4.5 or 8.5."), fp_NM2)}
 AnnualR
  RF_Down <-  block_list(
    fpar(ftext("Rainfall in the Year 2100" , fp_NM5)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext("Annual", fp_NM6)),
    fpar(ftext(paste0(RF100_4.5a," to ",RF100_8.5a,RFUnit,"/year"), fp_NM3)),
    fpar(ftext(AnnualR,fp_NM3)),
    # fpar(ftext("Dry Season", fp_NM6)),
    # fpar(ftext(paste0(MinRFD," to ",MaxRFD,RFUnit,"/month"), fp_NM3)),
    # fpar(ftext(paste0("(",DryR,")"),fp_NM3)),
    # fpar(ftext("Wet Season", fp_NM6)),
    # fpar(ftext(paste0(MinRFW," to ",MaxRFW,RFUnit,"/month"), fp_NM3)),
    # fpar(ftext(paste0("(",WetR,")"),fp_NM3)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(RF2100t))
  RF_Down  

  RF10085file <- paste0(R_FOLDER,"/",NameF,"/",NameF," DS_RF_8.5_v2.png")
  RF10085img <- external_img(src = RF10085file, height = 4,width = 4) 
  
  FIG_17<- block_list(
    fpar(ftext(paste0("Figure 26. Downscaled future rainfall projections (% change from present) ",
                      "at ",SNameF," by end-of-century (2100), using both Dynamical and Statistical ",
                      "downscaling. RCP 4.5 (top row) and RCP 8.5 (bottom row)."), fp_Fig)))

  

#################### Slide 25 Temperature 2040
  p<-p+1
  p29<-p
  
  # S23_TIT<- block_list(
  #   fpar(ftext("Average Air Temperature Change", FTXTT),fp_p = fp_par(text.align = "center")))
  
  S23_TIT<- block_list(
    fpar(ftext("Average Air Temperature Change", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext("Mid-Century (2040-2070)", FTXTT),fp_p = fp_par(text.align = "center")))
 
  # get AMMOUNT change (deg F)
  Down
  TAA_Thresh40_4.5 <-   Down[c(23),2]
  TAA_Thresh40_4.5
  TAA_Thresh40_8.5 <-   Down[c(24),2]
  
  TAA_RM  <- paste0(TAA_Thresh40_4.5," to ",TAA_Thresh40_8.5,"°F")
  TAA_RM
  
  # get percent change
  TAA <- CLIMA[3,14]
  TA40_4.5 <- round((TAA_Thresh40_4.5*100)/TAA,0)
  TA40_8.5 <- round((TAA_Thresh40_8.5*100)/TAA,0)
  
  TAAnnualRM  <- paste0(TA40_4.5," to ",TA40_8.5,"% change from present")
  TAAnnualRM
  
  # choose whether text says "increase" or "decrease"
  ifelse(TAA_Thresh40_4.5>=0, TA40_4.5_ch<-"increase", TA40_4.5_ch<-"decrease")
  ifelse(TAA_Thresh40_8.5>=0, TA40_8.5_ch<-"increase", TA40_8.5_ch<-"decrease")
  
  TA_Down4 <-  block_list(
    fpar(ftext("Air Temp. for Years 2040-2070" , fp_NM5)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext("Change in Air Temperature", fp_NM6)),
    fpar(ftext(TAA_RM, fp_NM3)),
    fpar(ftext(TAAnnualRM,fp_NM3)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext(paste0("These Statistical Downscaling maps show the projected change in air temperature under RCP 4.5 and 8.5 ",
                      "conditions. At ",SNameF,", air temperature is projected to ",TA40_4.5_ch," by ",abs(TAA_Thresh40_4.5),"°F (RCP 4.5), ",
                      "or ",TA40_8.5_ch," by ",abs(TAA_Thresh40_8.5),"°F (RCP 8.5) by mid-century."), fp_NM2)))
  TA_Down4
  
  TA40file <- paste0(R_FOLDER,"/",NameF,"/",NameF," StDs_Temp2040.png")
  TA40img <- external_img(src = TA40file) 
  
  FIG_20 <- block_list(
    fpar(ftext(paste0("Figure 27. Downscaled future temperature projections (°F change from present) ",
                      "at ",SNameF," by mid-century (2040-2070) using Statistical Downscaling. ",
                      "RCP 4.5 (left) and RCP 8.5 (right)."), fp_Fig)))
  
  #################### Slide 25 Temperature 2100
  p<-p+1
  p30<-p
  
  S22_TIT<- block_list(
    fpar(ftext("Average Air Temperature Change", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext("End-of-Century (2100)", FTXTT),fp_p = fp_par(text.align = "center")))
  
  # get ammount change (deg F)
  TAA_Thresh100_D4.5 <- Down[7,2]
  TAA_Thresh100_D8.5 <- Down[8,2]
  TAA_Thresh100_S4.5 <- Down[21,2]
  TAA_Thresh100_S8.5 <- Down[22,2]
  
  # Find largest change value for each RCP scenario
  ifelse(abs(TAA_Thresh100_D4.5)>=abs(TAA_Thresh100_S4.5), TA100_4.5c<-TAA_Thresh100_S4.5, TA100_4.5c<-TAA_Thresh100_D4.5)
  ifelse(abs(TAA_Thresh100_D8.5)>=abs(TAA_Thresh100_S8.5), TA100_8.5c<-TAA_Thresh100_D8.5, TA100_8.5c<-TAA_Thresh100_S8.5)
  
  # find minimum change
  TAA_Thresh100 <-  Down[c(7,8,21,22),2]
  TAA_100m<-min(TAA_Thresh100)
  TAA_100ma<-max(TAA_Thresh100)

  # get percent change
  TAA <- CLIMA[3,14]
  MinTA <- round((min(TAA_Thresh100*100))/TAA,0)
  MaxTA <- round((max(TAA_Thresh100*100))/TAA,0)

  TA100_4.5a <- round((TA100_4.5c*100)/TAA,0)
  TA100_8.5a <- round((TA100_8.5c*100)/TAA,0)
  
  # choose whether text says "increase" or "decrease"
  ifelse(TA100_4.5c>=0, TA100_4.5t<-"increase", TA100_4.5t<-"decrease")
  ifelse(TA100_8.5c>=0, TA100_8.5t<-"increase", TA100_8.5t<-"decrease")

  # # choose whether change value has + or - before it
  # ifelse(TA100_4.5c>=0, TA100_4.5c<-paste0("+ ",TA100_4.5c), TA100_4.5c<-TA100_4.5c)
  # ifelse(TA100_8.5c>=0, TA100_8.5c<-paste0("+ ",TA100_8.5c), TA100_8.5c<-TA100_8.5c)
  
  AnnualT  <- paste(MinTA, " to ", MaxTA, "% change from present")
  AnnualT
  # DryR  <- paste(min(TAD_Thresh100), " to ", max(TAD_Thresh100), "% Change")
  # WetR  <- paste(min(TAW_Thresh100), " to ", max(TAW_Thresh100), "% Change") 
  #  
  
  # determine if the downscaling method outputs agree on direction of change, and choose the appropriate text
  ifelse(TAA_Thresh100_D4.5>=0, a<-1, a<-2)
  ifelse(TAA_Thresh100_S4.5>=0, b<-1, b<-2)
  ifelse(TAA_Thresh100_D8.5>=0, c<-1, c<-2)
  ifelse(TAA_Thresh100_S8.5>=0, d<-1, d<-2)
  
  
  if(a==b && c==d)
  {TA2100t<-ftext(paste0("These Dynamical and Statistical Downscaling maps show the projected ",
                         "change in air temperature under RCP 4.5 and 8.5 conditions. At ",SNameF,
                         ", air temperature is projected to ",TA100_4.5t," by ",abs(TA100_4.5c),
                         "°F (RCP 4.5), or ",TA100_8.5t," by ",abs(TA100_8.5c),"°F (RCP 8.5) by end-of-century."), fp_NM2)}
  if(a!=b && c==d)
  {TA2100t<-ftext(paste0("These Dynamical and Statistical Downscaling maps show the projected ",
                         "change in rainfall under RCP 4.5 and 8.5 conditions. At ",SNameF,
                         ", the models do not agree on the direction of change for RCP 4.5. Annual ",
                         "rainfall is projected to ",TA100_8.5t," by ",abs(TA100_8.5c),"°F (RCP 8.5) by end-of-century."), fp_NM2)}
  if(a==b && c!=d)
  {TA2100t<-ftext(paste0("These Dynamical and Statistical Downscaling maps show the projected ",
                         "change in rainfall under RCP 4.5 and 8.5 conditions. At ",SNameF,
                         ", annual rainfall is projected to ",TA100_4.5t," by ",abs(TA100_4.5c),
                         "°F (RCP 4.5). The models do not agree on the direction of change for RCP 8.5."), fp_NM2)}
  if(a!=b && c!=d)
  {TA2100t<-ftext(paste0("These Dynamical and Statistical Downscaling maps show the projected ",
                         "change in rainfall under RCP 4.5 and 8.5 conditions. At ",SNameF,
                         ", the models do not agree on the overall direction of change for either RCP 4.5 or 8.5."), fp_NM2)}
  TA2100t
  
  TA_Down <-  block_list(
    fpar(ftext("Air Temp. in the Year 2100" , fp_NM5)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext("Change in Air Temperature", fp_NM6)),
    fpar(ftext(paste0(TA100_4.5c," to ",TA100_8.5c,"°F "), fp_NM3)),
    fpar(ftext(AnnualT, fp_NM3)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(TA2100t))
  TA_Down  
    
  RF10085file <- paste0(R_FOLDER,"/",NameF,"/",NameF," DS_RF_8.5_v2.png")
  RF10085img <- external_img(src = RF10085file, height = 4,width = 4) 
  
  TA100file <- paste0(R_FOLDER,"/",NameF,"/",NameF," DS_Temp2100.png")
  TA100img <- external_img(src = TA100file)
  
  FIG_19<- block_list(
    fpar(ftext(paste0("Figure 28. Downscaled future rainfall projections (% change from present) ",
                      "at ",SNameF," by end-of-century (2100), using both Dynamical and Statistical ",
                      "downscaling. RCP 4.5 (top row) and RCP 8.5 (bottom row)."), fp_Fig)))
  
  
  ################ Slide 27 Summary and Conclusions ###############
  p<-p+1
  p31<-p
  
  S24_TIT<- block_list(
    fpar(ftext("Part 6: CCVD Summary", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext(SNameF, FTXTT3),fp_p = fp_par(text.align = "center")))
  
  
  DecD <- (CNTDRT/10)
  
  Summary<- paste0(NameF," (",SNameF,") is located on the island of ", ISL, " at mean elevation of ", El_Mean, ELUnit," (rangF: ",El_Min, " to ", El_Max,ELUnit,"). ",  
                 "Rainfall varies over the course of the year with a maximum of ",RFx, RFUnit3," occurring in ", MONLIST[RFxm], " and a minimum of ",RFn, RFUnit3,
                 " occurring in ", MONLIST[RFnm],". On average, wet season months (Nov-Apr) receive ",SeaD,RFUnit," more rainfall than dry season months (May-Oct). ",
                 "Seasonal rainfall can vary within the unit as well, with dry season rainfall ranging from ",D_RF_R, RFUnit3," and wet season rainfall ",
                 "ranging from ",W_RF_R, RFUnit3," across the ",EL_Dif,ELUnit, " elevation gradient. ",
                 "Rainfall can also vary considerably from year-to-year with the driest years occurring during a Strong El Niño event, when on ",
                 "average, ", DifWRF,"% less rainfall is received, relative to the long-term average. ",
                 "The average temperature at ",SNameF," is ", TAA,TUnit," but temperature ranges from ",TAn, TUnit, " to ", TAx, TUnit," over the course of the year. ",
                 "Drought is a reoccurring feature in the climate system of ", SNameF, " with a total of ",
                 CNTDRT, " occurring over the record which is approximately ",DecD," per decade. A total of ",SoG, " drought events were at severe strength ",
                 "or greater and the longest drought lasted for a total of ", LongD, " consecutive months. ",
                 "Future projections of rainfall are uncertain, with end-of-century annual changes ranging from ", min(RFA_Thresh100)," to " ,max(RFA_Thresh100),
                 ". Future projections of temperature suggest an increase of ", TAA_Thresh40_4.5,TUnit, 
                 " to ", TAA_Thresh40_8.5,TUnit," by mid century (2040-2070) ",
                 "and an increase of ", TAA_100m,TUnit," to ",  TAA_100ma, TUnit," by the end of the century (2100). ")
                 
                 fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 17) 
                 fp_Summary <- fpar(ftext(Summary, fp_Tx))

 
   ################ Slide 28 Resource PAGE ############### 
                 p<-p+1
                 p32<-p
                 
                 S25_TIT<- block_list(
                   fpar(ftext("External Resources", FTXTT),fp_p = fp_par(text.align = "center")))
                 
                 fp_NM3a <- fp_text(font.size = 14, color = "darkgreen")
                 
                 fp_NM6a<- fp_text(font.size = 14, color = "black",bold = TRUE)
                 
                 fp_RES <-  block_list(
                   fpar(ftext( "For more Information"   , fp_NM5)),
                   fpar(ftext( "                    "   , fp_Fig2)),
                   fpar(ftext("US Drought Monitor"      , fp_NM6a)),
                   fpar(ftext("https://droughtmonitor.unl.edu/", fp_NM3a)),
                   fpar(ftext( "                    ", fp_Fig2)),
                   fpar(ftext("NIDIS Current Hawaii Drought Maps"      , fp_NM6a)),
                   fpar(ftext("https://www.drought.gov/states/hawaii", fp_NM3a)),
                   fpar(ftext( "                    ", fp_Fig2)),
                   # fpar(ftext("National Drought Mitigation Center (NDMC)" , fp_NM6a)),
                   # fpar(ftext("https://drought.unl.edu/ranchplan/Overview.aspx", fp_NM3a)),
                   # fpar(ftext( "     "   , fp_Fig2)),
                   fpar(ftext("State of Hawaii Drought Plan"      , fp_NM6a)),
                   fpar(ftext(" https://files.hawaii.gov/dlnr/cwrm/planning/HDP2017.pdf", fp_NM3a)),
                   fpar(ftext( "                    ", fp_Fig2)),
                   fpar(ftext( "Hawaii Climate Data Portal", fp_NM6a)),
                   fpar(ftext("https://www.hawaii.edu/climate-data-portal/", fp_NM3a)),
                   fpar(ftext( "                    ",   fp_Fig2)),
                   fpar(ftext("ENSO Current Phase and Discussion" , fp_NM6a)),
                   fpar(ftext("https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso_advisory/ensodisc.shtml", fp_NM3a)),
                   fpar(ftext( "     "   , fp_Fig2)),
                   fpar(ftext( "Pacific Drought Knowledge Exchange", fp_NM6a)),
                   fpar(ftext("http://www.soest.hawaii.edu/pdke/", fp_NM3a)),
                   fpar(ftext( "                    ", fp_Fig2)),
                   fpar(ftext( "Pacific Fire Exchange", fp_NM6a)),
                   fpar(ftext("https://www.pacificfireexchange.org/", fp_NM3a)),
                   fpar(ftext( "                    ", fp_Fig2)),
                   fpar(ftext( "Ahupuaa GIS Layer Storymap", fp_NM6a)),
                   fpar(ftext("https://storymaps.arcgis.com/stories/1ad9fcf7f2c345a58adef0997fce9b5d", fp_NM3a)),
                   fpar(ftext( "                    ", fp_Fig2)))
                 
                 #ph_hyperlink(doc,
                 #             ph_label = "Content Placeholder 2", href = "https://cran.r-project.org")
                 
  ################ Slide 29 Acknowledgements ###############
                 p<-p+1
                 p33<-p
                 
                 S26_TIT<- block_list(
                   fpar(ftext("Acknowledgements", FTXTT),fp_p = fp_par(text.align = "center")))
       
  
  Ack <- paste0( "Ryan Longman, Abby Frazier, Christian Giardina, Derek Ford, Keri Kodama, Alyssa Anderson (PDKE), Mari-Vaughn Johnson, Heather Kerkering, Patrick Grady, Elliott Parsons (PI-CASC), ",
                    "Clay Trauernicht, Melissa Kunz, Cherryle Heu (NREM, UH Manoa), Amanda Sheffield, Britt Parker, John Marra (NOAA), ",
                    "Sean Cleveland, Jared McLean (ITS, UH Manoa), Katie Kamelamela (Arizona State University), Thomas Giambelluca (WRRC, UH Manoa), ",
                    "Jim Poterma (SOEST, UH Manoa), Carolyn Auweloa (NRCS), Nicole Galase (Hawaii Cattlemen's Council), Mark Thorne (CTAHR, UH Manoa)")
  
  fp_Ack2 <- fp_text(color = "black", font.size = 15) 
  fp_Ack <- fpar(ftext(Ack, fp_Ack2 ))
  
  
  fp_CITA1 <- fp_text(bold= TRUE, color = "darkred", font.size = 16,underlined = TRUE) 
  fp_CITA2 <- fp_text(color = "black", font.size = 16) 
  fp_CITA3 <- fp_text(color = "red", font.size = 14) 
  
  CITA <- block_list(
    fpar(ftext       ("Suggested Citation", fp_CITA1)),
    fpar(ftext(paste0("Longman, R.J., Ford, D.J., Anderson, A., Frazier, A.G, Giardina, C.P (2023). ",
                      "Climate Change, Climate Variability and Drought Portfolio ", NameF,
                      " Pacific Drought knowledge Exchange CCVD Series Version ",ver,"."), fp_CITA2)))


  Draft <- block_list(
    fpar(ftext       ("For the most up-to-date version of this portfolio contact Derek Ford: fordd@eastwestcenter.org for more information.", fp_CITA3)))

   ############## Slide 30 Works Cited ###########
  p<-p+1
  p34<-p
  
  S27_TIT<- block_list(
    fpar(ftext("Works Cited", FTXTT),fp_p = fp_par(text.align = "center")))
  
  Worksfile <- paste0(I_FOLDER,"Works2.jpg")
  Worksimg <- external_img(src =Worksfile , height = 1, width = 1) 
  
  
  
  ############### Slide 31 ANNEX I ###############
  p<-p+1
  p35<-p
  
  RF100<- paste0("The 100+ year monthly rainfall dataset was drawn from two unique gridded products. We used data from Frazier et al. (2016) ",
                 "for the period 1920-1989 and Lucas et al. (2022) for the period 1990-",yr,". Given that two unique data sets and methods ",
                 "were used to make these two products we show the 1:1 Statistical relationship between the two products for a 23-year overlap ",
                 "(1990-2012) with the datasets and associated error metrics.")
  
  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 18) 
  fp_RF100 <- fpar(ftext(RF100, fp_Tx))
  
  SCAfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," 23yr_RF_Compare.png")
  SCAimg <- external_img(src = SCAfile) 
  
  FIG_A1 <- block_list(
    fpar(ftext(paste0("Figure A1: One to one comparison of 23-years (1990-2012) of monthly rainfall ",  
                      "from two unique datasets for ",SNameF,", and associated error metrics; R2, is the coefficient of determination, ",  
                      "MBE, is the mean bias error, MAE, mean absolute error."), fp_Fig)))
  

  
    ################ Slide 32 ANNEX II #############
  p<-p+1
  p36<-p

  DownA<- paste0("Two types of downscaling products were used in this analysis. Here we explain some of the nuances between the two. ",
                 "Dynamical Downscaling (Zhang et al., 2016), feeds GCM output into a regional model that can account for local topographic and ",
                 "atmospheric phenomena at much finer resolutions (e.g. 1 km). End-of-century projections (2100) encompass the period 2080-2099. ",
                 "Statistical Downscaling ",
                 "(Timm et al., 2015, Timm, 2017), develops a relationship between GCM model output and station data for a historical ",
                 "period and then uses this established relationship to make projections for two future scenarios. End-of-century projections (2100) ",
                 "encompass the period 2070-2099 (2100), Mid-century projections encompass the period 2040-2070.")

       
  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 18) 
  fp_DownA <- fpar(ftext(DownA, fp_Tx))
  
 
  
  # ################ Slide 33 ANNEX III #############
  # 
  # Roadfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Forest_Roads.png")
  # 
  # Roadimg <- external_img(src = Roadfile) 
  # 
  # Trailfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Trails.png")
  # Trailimg <- external_img(src = Trailfile) 
  # 
  # FIG_A2 <- block_list(
  #   fpar(ftext(paste0("Figure A2: Forest roads ",ISL,"."), fp_Fig)))
  # FIG_A3 <- block_list(
  #   fpar(ftext(paste0("Figure A3: Trail Inventory ",ISL,"."), fp_Fig)))
  
  
  #################### History of Drought Events table

  p<-p+1
  p37<-p
  
  DT_TIT<- block_list(
    fpar(ftext(paste0("Drought Events (1920-",yr,")"), FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext(SNameF, FTXTT3),                      fp_p = fp_par(text.align = "center")))
  
  
  TAB3 <- block_list(
    fpar(ftext(paste0("Table A1. SPI-12 drought characteristics at ",NameF, " identified in the SPI-12 timeseries. Duration is the number of months the drought persisted; ",
                      "Average Intensity is the average absolute SPI; Peak Intensity is the highest SPI value calculated during the drought ", 
                      "Magnitude is sum of absolute SPI values during the drought."), fp_Fig)))
  
  #NOTE Try and do the Auto FIt for this table to see if you can get it on the slide 
  #See Slide 10


  
  
############################################################################################################################################################  
###############
###############  SLIDE DECK 
###############

  #Look up office themes for R PPT 
  
#Slide 1
  mypowerpoint <- read_pptx() %>%
    add_slide("Title Slide","Office Theme") %>%
    #Add Title
       ph_with(TIT,ph_location_type("ctrTitle")) %>%
    #Add dates left footer NOTE#Check ph_location_type
       ph_with(value = format(Sys.Date()), location = ph_location_type(type = "dt"))%>%
    #Add subtitle 
       ph_with(SUB,ph_location_type("subTitle")) %>%
    # #add logos
    #    ph_with(value = EWCimg, location = ph_location(label = "my_name",
    #                   left = 4.38, top = 6.0, width = 1.4, height = 1.4)) %>%
    #    ph_with(value = CASCimg, location = ph_location(label = "my_name",
    #                   left = 2.1, top = 6.2, width = 1, height = 1)) %>%
    #    ph_with(value = FSimg, location = ph_location(label = "my_name",
    #                   left = 6.9, top = 6.2, width = 1, height = 1)) %>%
       ph_with(value = PDKE_L, location = ph_location(label = "my_name",
                      left = 0, top = 0.3, width = 10, height = 2))%>%
       ph_with(value = EWClimg, location = ph_location(label = "my_name",
                      left = 6.6, top = 6.6, width = 2.9, height = 0.6))%>%

    #Slide 2 
    add_slide("Title and Content","Office Theme") %>%
      #ph_with(S2_TIT, ph_location_type("title",position_left = TRUE)) %>%
      ph_with(fp_HDKE,ph_location_type("body"))%>%
      ph_with(value = p1, location = ph_location_type(type = "sldNum"))%>%
      # Logos
      # ph_with(value = EWCimg, location = ph_location(label = "my_name",
      #                 left = 1.1, top = 6.1, width = 1.4, height = 1.4)) %>%
      # ph_with(value = CASCimg, location = ph_location(label = "my_name",
      #                 left = 2.6, top = 6.3, width = 1, height = 1)) %>%
      # ph_with(value = FSimg, location = ph_location(label = "my_name",
      #                 left = 3.8, top = 6.3, width = 1, height = 1)) %>%
      # # ph_with(value = RISAimg, location = ph_location(label = "my_name",
      # #                 left = 4.2, top = 6.3, width = 1.4, height = 1)) %>%
      # ph_with(value = WRRCimg, location = ph_location(label = "my_name",
      #                 left = 5, top = 6.3, width = 1.1, height = 1)) %>%
      # ph_with(value = NOAAimg, location = ph_location(label = "my_name",
      #                 left = 6.3, top = 6.3, width = 1, height = 1)) %>%
      # ph_with(value = UHimg, location = ph_location(label = "my_name",
      #                 left = 7.5, top = 6.3, width = 1, height = 1)) %>%
     ph_with(value = PDKE_L, location = ph_location(label = "my_name",
                      left = 0, top = 0, width = 10, height = 2))%>%
     ph_with(value = "Longman et al. (2022)", location = ph_location_type(type = "dt"))%>%
    

    #Slide 3
        add_slide("Two Content","Office Theme") %>%
          ph_with(S3_TIT,         ph_location_type("title")) %>%
          # ph_with(fp_CCVD,        ph_location_type("body",position_right = FALSE)) %>%
          ph_with(value = fp_CCVD, location = ph_location(label = "my_name", left = 0.5, top = 1.1, width = 5.1, height = 6)) %>%
          ph_with(value = p2, location = ph_location_type(type = "sldNum")) %>%
          # ph_with(value = MAPimg, ph_location_type("body",position_right = TRUE)) %>%
          ph_with(value = D_AHimg, location = ph_location(label = "my_name",
                                                      left = 5.85, top = 1.75, height = 4.5, width = 3.3)) %>%
          ph_with(value = FIG_1, location = ph_location(label = "my_name",
                          left = 5.8, top = 6.05, width = 3.3, height = 0.77))%>%
    
# #Slide 3 
#     add_slide("Two Content","Office Theme") %>%
#       ph_with(S3_TIT,         ph_location_type("title")) %>%
#       ph_with(fp_CCVD,        ph_location_type("body",position_right = FALSE)) %>%
#       ph_with(value = "3", location = ph_location_type(type = "sldNum")) %>%
#       # ph_with(value = MAPimg, ph_location_type("body",position_right = TRUE)) %>%
#     ph_with(value = MAPimg, location = ph_location(label = "my_name",
#                                                   left = 5, top = 1.7, height = 4.6, width = 4.6)) %>%
#       ph_with(value = FIG_1, location = ph_location(label = "my_name",
#                       left = 5.6, top = 6.4, width = 3.81, height = 0.77))%>%
    
#Slide 4 Land Divisions
    add_slide("Title and Content","Office Theme") %>%
    ph_with(S4a_TIT,ph_location_type("title",position_left = TRUE)) %>%
    ph_with(fp_HLD1, location = ph_location(label = "my_name",
                                           left = 0.5, top = 1, width = 9, height = 2)) %>%
    ph_with(value = p3, location = ph_location_type(type = "sldNum"))%>%
    ph_with(value = D_HLDimg, location = ph_location(label = "my_name",
                                                   left = 0.39, top = 2.7, width = 9.3, height = 4.3)) %>%
    ph_with(value = MPmimg, location = ph_location(label = "my_name",
                                                  left = 0.565, top = 3.335, width = 2.83, height = 2.2)) %>%
    ph_with(value = MOmimg, location = ph_location(label = "my_name",
                                                  left = 3.6, top = 3.335, width = 2.83, height = 2.2)) %>%
    ph_with(value = AHmimg, location = ph_location(label = "my_name",
                                                   left = 6.68, top = 3.335, width = 2.83, height = 2.2)) %>%
    ph_with(value = MP, location = ph_location(label = "my_name",
                                                   left = 0.565, top = 4.6, width = 2.85, height = 3)) %>%
    ph_with(value = MO, location = ph_location(label = "my_name",
                                                   left = 3.6, top = 4.6, width = 2.85, height = 3)) %>%
    ph_with(value = AH, location = ph_location(label = "my_name",
                                                   left = 6.68, top = 4.6, width = 2.85, height = 3)) %>%
    ph_with(value = ft1, location = ph_location(label = "my_name",
                                                left = 0.5, top = 5.73, width = 3, height = 3)) %>%
  
#Slide 5 
    add_slide("Two Content","Office Theme") %>%
    ph_with(S5_TIT,       ph_location_type("title")) %>%
    ph_with(M_Elev,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(value = p4, location = ph_location_type(type = "sldNum")) %>%
    ph_with(value = Elevimg, location = ph_location(label = "my_name",
                                                   left = 5, top = 1.7, height = 4.4, width = 4.4)) %>%
    ph_with(value = FIG_2, location = ph_location(label = "my_name",
                                                  left = 5.1, top = 5.6, width = 4.2, height = 1.4))%>%
      ph_with(value = ft2, location = ph_location(label = "my_name",
                                                  left = 0.5, top = 5.73, width = 3, height = 3)) %>%

#Slide 6.1
  add_slide("Two Content","Office Theme") %>%
  ph_with(S6.1_TIT,       ph_location_type("title")) %>%
  # ph_with(M_LAND,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = M_LAND, location = ph_location(label = "my_name", left = 0.5, top = 1.1, width = 4.65, height = 3)) %>%
  ph_with(value = p5, location = ph_location_type(type = "sldNum")) %>%
  ph_with(value = LClegimg, location = ph_location(label = "my_name", left = 7.5, top = 1.4, width = 1.5, height = 1.5)) %>%
  ph_with(value = LCbarimg, location = ph_location(label = "my_name", left = 0.5, top = 3.75, width = 4.7, height = 2.8)) %>%
  ph_with(value = LCmapimg, location = ph_location(label = "my_name",
                                                    left = 5.5, top = 3.2)) %>%
  ph_with(value = FIG_3.1, location = ph_location(label = "my_name",
                                                left = 0.5, top = 5.9, width = 4.2, height = 1.4))%>%
  ph_with(value = FIG_3.2, location = ph_location(label = "my_name",
                                                  left = 5.5, top = 6, width = 3.5, height = 1.4))%>%
      ph_with(value = ft3, location = ph_location(label = "my_name",
                                                  left = 0.5, top = 5.73, width = 4, height = 3)) %>%

    
  #Slide 7.1 Aquifers
  add_slide("Two Content","Office Theme") %>%
  ph_with(S6.2_TIT,       ph_location_type("title")) %>%
  ph_with(value = AQ_T, location = ph_location(label = "my_name", left = 0.6, top = 4, width = 3.5, height = 3)) %>%
  ph_with(value = p6, location = ph_location_type(type = "sldNum")) %>%
  ph_with(value = M_AQ, location = ph_location(label = "my_name", left = 0.5, top = 0.8, width = 5.5, height = 4)) %>%
  ph_with(value = AQmapimg, location = ph_location(label = "my_name", left = 6, top = 2.3)) %>%
  ph_with(value = FIG_5.b, location = ph_location(label = "my_name",
                                                  left = 6.5, top = 5.9, width = 3.2, height = 1.3))%>%
  ph_with(value = TAB_0, location = ph_location(label = "my_name",
                                                  left = 0.5, top = 6, width = 3.5, height = 1.4))%>%
  # ph_with(value = ft3, location = ph_location(label = "my_name",
  #                                             left = 0.5, top = 5.73, width = 4, height = 3)) %>%

  #Slide 7.2 Streams
  add_slide("Two Content","Office Theme") %>%
  ph_with(S6.2_TIT,       ph_location_type("title")) %>%
  ph_with(value = M_HF, location = ph_location(label = "my_name", left = 0.5, top = 2, width = 5.4, height = 3)) %>%
  ph_with(value = p7, location = ph_location_type(type = "sldNum")) %>%
  ph_with(value = HFlmapimg, location = ph_location(label = "my_name", left = 2, top = 5.2, width = 2.8, height = 1.9)) %>%
  ph_with(value = HFmapimg, location = ph_location(label = "my_name", left = 5.8, top = 2.3)) %>%
  ph_with(value = FIG_6.b, location = ph_location(label = "my_name",
                                                  left = 6.2, top = 5.3, width = 3.2, height = 2))%>%

  #Slide 7  PART 2
    add_slide("Title and Content","Office Theme") %>%
    ph_with(P1_TIT,ph_location_type("title",position_left = TRUE)) %>%
    ph_with(fp_PART1,ph_location_type("body"))%>%
    ph_with(value = p8, location = ph_location_type(type = "sldNum"))%>%
    ph_with(value = MODimg, location = ph_location(label = "my_name",
                                                   left = 1.1, top = 4.5, width = 2.6, height = 2)) %>%
    ph_with(value = RGimg, location = ph_location(label = "my_name",
                                                    left = 3.8, top = 4.5, width = 2.4, height = 2)) %>%
    ph_with(value = RFimg, location = ph_location(label = "my_name",
                                                    left = 6.4, top = 4.5, width = 2.4, height = 2)) %>%

#Slide 8a    
  add_slide("Two Content","Office Theme") %>%
    ph_with(S8.0_TIT,           ph_location_type("title")) %>%
  
    ph_with(value = RS_T, location = ph_location(label = "my_name",
                                                   left = 0.6, top = 4.9, width = 3, height = 1.4))%>%
    ph_with(value = RSimg, location = ph_location(label = "my_name",
                                                    left = 5.3, top = 1.5, width = 4.3, height = 3.2))%>%
    ph_with(value = fp_RS1, ph_location_type("body",position_right = FALSE)) %>%
    ph_with(value = FIG_5.a, location = ph_location(label = "my_name",
                                                    left = 8, top = 4.5, width = 1.5, height = 2))%>%
    ph_with(value = FIG_6.a, location = ph_location(label = "my_name",
                                                    left = 0.5, top = 6.3, width = 7, height = 1.4))%>%
    ph_with(value = p9, location = ph_location_type(type = "sldNum")) %>%
    
#Slide 6 
  add_slide("Two Content","Office Theme") %>%
    ph_with(S6_TIT,           ph_location_type("title")) %>%
    ph_with(value = M_ANNCLIM, location = ph_location(label = "my_name",
                     left = 0.5, top = 0.05, width = 4.7, height = 5))%>%
    ph_with(value = Tt, location = ph_location(label = "my_name",
                     left = 0.8, top = 3.9, width = 1, height = 1))%>%
    ph_with(value = p10, location = ph_location_type(type = "sldNum")) %>%
    # ph_with(value = Climimg, ph_location_type("body",position_right = TRUE)) %>%
    ph_with(value = "Giambelluca et al. (2013;2014)", location = ph_location_type(type = "dt"))%>%
  
    ph_with(value = RF2img, location = ph_location(label = "my_name",
                     left = 5.4, top = 1.2, width=si, height= si))%>%
    ph_with(value = TAAimg, location = ph_location(label = "my_name",
                     left = 7.5, top = 1.2, width=si, height=si))%>%
    ph_with(value = RHimg, location = ph_location(label = "my_name",
                     left = 5.4, top = 2.7, width=si, height=si))%>%
    ph_with(value = SRimg, location = ph_location(label = "my_name",
                     left = 7.5, top = 2.7, width=si, height=si))%>%
    ph_with(value = SMimg, location = ph_location(label = "my_name",
                     left = 5.4, top = 4.2, width=si, height=si))%>%
    ph_with(value = ETimg, location = ph_location(label = "my_name",
                     left = 7.5, top = 4.2, width=si, height=si))%>%
    ph_with(value = WSimg, location = ph_location(label = "my_name",
                     left = 5.4, top = 5.7, width=si, height=si))%>%
    
  ph_with(value = FIG_5.1, location = ph_location(label = "my_name",
                                                left = 0.5, top = 6.1, width = 4.6, height = 1.4))%>%
  
  ph_with(value = FIG_6.1, location = ph_location(label = "my_name",
                                                left = 7.5, top = 5.6, width = si, height = si))%>%
  
#Slide 8
  add_slide("Two Content","Office Theme") %>%
    ph_with(S7_TIT,          ph_location_type("title")) %>%
    ph_with(value = RFBimg, location = ph_location(label = "my_name",
                                                   left = 0.3, top = 3.2, width = 6.3, height = 3.7)) %>%
    ph_with(fp_CLIMO_RF, location = ph_location(label = "my_name",
                                             left = 0.4, top = 1.1, width = 6.1, height = 3)) %>%
    ph_with(value = p11, location = ph_location_type(type = "sldNum")) %>%
    ph_with(value = RFWimg, location = ph_location(label = "my_name",
                                                   left = 6.9, top = 1.35, width = 2.55, height = 2.55)) %>%
    ph_with(value = RFDimg, location = ph_location(label = "my_name",
                                                   left = 6.9, top = 3.9, width = 2.55, height = 2.55)) %>%
    ph_with(value = "Giambelluca et al. (2013;2014)", location = ph_location_type(type = "dt"))%>%
    ph_with(value = FIG_6a, location = ph_location(label = "my_name",
                                                  left = 0.7, top = 5.9, width = 4, height = 1.4))%>%
    ph_with(value = FIG_7a, location = ph_location(label = "my_name",
                                                   left = 6.9, top = 6, width = 2.4, height = 1.4))%>%
   
#Slide 9
  add_slide("Two Content","Office Theme") %>%
    ph_with(S8a_TIT,          ph_location_type("title")) %>%
    ph_with(value = TALimg, location = ph_location(label = "my_name",
                                                   left = 0.3, top = 3.1, width = 6.3, height = 3.7)) %>%
    ph_with(fp_CLIMO_TA, location = ph_location(label = "my_name",
                                             left = 0.4, top = 1.1, width = 6.1, height = 3)) %>%
    ph_with(value = p12, location = ph_location_type(type = "sldNum")) %>%
    ph_with(value = TACimg, location = ph_location(label = "my_name",
                                                   left = 6.9, top = 1.35, width = 2.55, height = 2.55)) %>%
    ph_with(value = TAHimg, location = ph_location(label = "my_name",
                                                   left = 6.9, top = 3.9, width = 2.55, height = 2.55)) %>%
    ph_with(value = "Giambelluca et al. (2013;2014)", location = ph_location_type(type = "dt"))%>%
    ph_with(value = FIG_8a, location = ph_location(label = "my_name",
                                                   left = 0.7, top = 6.05, width = 4, height = 1.4))%>%
    ph_with(value = FIG_9a, location = ph_location(label = "my_name",
                                                   left = 6.9, top = 6, width = 2.6, height = 1.4))%>%
    
#Slide 9
add_slide("Two Content","Office Theme") %>%
  ph_with(S8_TIT, ph_location_type("title")) %>%
  ph_with(value = RFFimg, ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = p13, location = ph_location_type(type = "sldNum")) %>%
  ph_with(value = TAimg, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_5, location = ph_location(label = "my_name",
                                                left = 0.5, top = 6.08, width = 4.5, height = 1.4))%>%
  ph_with(value = FIG_6, location = ph_location(label = "my_name",
                                                left = 5.1, top = 6.08, width = 4.5, height = 1.4))%>%
  ph_with(my_pres, value = "Giambelluca et al. (2013;2014)", location = ph_location_type(type = "dt"))%>%

#Slide 10 
    
add_slide("Two Content","Office Theme") %>%
  ph_with(S9_TIT,      ph_location_type("title")) %>%
  ph_with(fp_SEA,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = p14, location = ph_location_type(type = "sldNum")) %>%
  ph_with(value = SEAimg, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_12.1, location = ph_location(label = "my_name",
                                                left = 6, top = 6.2, width = 3, height = 1.4))%>%
  ph_with(value = "Giambelluca et al. (2013)", location = ph_location_type(type = "dt"))%>%
  
  
  #Slide 10
  add_slide("Title and Content","Office Theme") %>%
    ph_with(S10_TIT,      ph_location_type("title")) %>%
    ph_with(value = p15, location = ph_location_type(type = "sldNum"))%>% 
    #ph_with(value = ACLIM_T, ph_location_type("body")) %>%
    ph_with(value = "Giambelluca et al. (2013;2014)", location = ph_location_type(type = "dt"))%>%
    ph_with(value =   TAB1, location = ph_location(label = "my_name",
                                                   left = 0.2, top = 4.45, width = 9.5, height = 4))%>%
    ph_with(value =   ACLIM_T, location = ph_location(label = "my_name",
                                                      left = 0.3, top = 1.8, width = 7, height = 3.5))%>%

#Slide 12.1
add_slide("Two Content","Office Theme") %>%
  ph_with(P2_TIT,ph_location_type("title",position_left = TRUE)) %>%
  ph_with(value = InterAn, location = ph_location(label = "my_name",
                                                    left = 0.6, top = 1.3, width = 5, height = 3)) %>%
  ph_with(value = p16, location = ph_location_type(type = "sldNum"))%>%
  ph_with(value = ENSO2img, location = ph_location(label = "my_name",
                                           left = 6, top = 1.5, width = 3.2, height = 4.4)) %>%
  ph_with(value = MEIimg, location = ph_location(label = "my_name",
                                                 left = .3, top = 4.15, width = 5.3, height = 2.7)) %>%
  ph_with(value = FIG_10.1, location = ph_location(label = "my_name",
                                                  left = 6, top = 5.9, width = 3.4, height = 1))%>%
    ph_with(value = ft4, location = ph_location(label = "my_name",
                                                left = 0.5, top = 5.73, width = 4, height = 3)) %>%
  
  #Slide 12.2
add_slide("Two Content","Office Theme") %>%
  ph_with(EN_TIT,ph_location_type("title")) %>%
  ph_with(value = MEI3, location = ph_location(label = "my_name", left = 0.4, top = .25, width = 9, height = 5)) %>%
  ph_with(value = p17, location = ph_location_type(type = "sldNum")) %>%
  ph_with(value = MEISimg, location = ph_location(label = "my_name",
                                                   left = 2, top = 4.15, width = 6, height = 2.5)) %>%
  ph_with(value = FIG_14.1, location = ph_location(label = "my_name",
                                                left = 2, top = 6.45, width = 5.5, height = 1))%>%

  #Slide 14 
  add_slide("Two Content","Office Theme") %>%
  ph_with(S13_TIT,         ph_location_type("title")) %>%
  ph_with(value = fp_Trend, location = ph_location(label = "my_name",
                                                   left = 0.5, top = 0.3, width = 4.5, height = 5.1)) %>%
  ph_with(value = Trendimg, location = ph_location(label = "my_name",
                                                   left = 5.2, top = 1.4, width = 4.5, height = 5.1))%>%
  ph_with(value = RFT_T, location = ph_location(label = "my_name",
                                                   left = 1, top = 4.5, width = 4.5, height = 5.1))%>%
  ph_with(value = FIG_9, location = ph_location(label = "my_name",
                                                left = 5.3, top = 6.1, width = 4.3, height = 1.4))%>%
 ph_with(value = TAB2.1, location = ph_location(label = "my_name",
                                                  left = 0.6, top = 5.9, width = 4, height = 1.4))%>%
  ph_with(my_pres, value = "Frazier et al. (2016); Lucas et al. (2022); See Annex I", location = ph_location_type(type = "dt"))%>%
  ph_with(value = p18, location = ph_location_type(type = "sldNum"))%>%

  #Slide 14 (NEW - Air Temp)  
  add_slide("Two Content","Office Theme") %>%
  ph_with(S14_TIT,         ph_location_type("title")) %>%

  ph_with(value = fp_AirTrend, location = ph_location(label = "my_name", left = 0.5, top = 0.05, width = 6, height = 5)) %>%
  ph_with(value = fp_AirTrend2, location = ph_location(label = "my_name", left = 0.5, top = 2.5, width = 3, height = 5)) %>%  
  ph_with(value = Trendimg2, location = ph_location(label = "my_name", left = 7.1, top = 1.7, width = 2.5, height = 1.7)) %>%
  ph_with(value = Trendimg1, location = ph_location(label = "my_name", left = 4, top = 3.5, width = 5.5, height = 3)) %>%
    
  ph_with(value = FIG_10_11, location = ph_location(label = "my_name", 
                                                left = 4.3, top = 6, width = 5, height = 1.4))%>%
  ph_with(value = p19, location = ph_location_type(type = "sldNum"))%>%
  ph_with(value = ft5, location = ph_location(label = "my_name",
                                              left = 0.5, top = 5.73, width = 4, height = 3)) %>%

#Slide 15 #Part 3
  add_slide("Title and Content","Office Theme") %>%
  ph_with(P3_TIT,ph_location_type("title",position_left = TRUE)) %>%
  ph_with(fp_Part3,ph_location_type("body"))%>%
  ph_with(value = p20, location = ph_location_type(type = "sldNum"))%>%
  ph_with(value = Droughtimg, location = ph_location(label = "my_name",
                                                 left = 1.9, top = 4.5, width = 3, height = 2)) %>%
  ph_with(value = Fireimg, location = ph_location(label = "my_name",
                                                 left = 5.45, top = 4.5, width = 3, height = 2)) %>%
    
###  Slide 16    
    add_slide("Title and Content","Office Theme") %>%
    ph_with(TIT_D,ph_location_type("title",position_left = TRUE)) %>%
    # ph_with(value = D_METimg, location = ph_location(label = "my_name",
    #                                                  left = 0.3, top = 1.8, width = 2, height = 1.3)) %>%
    # ph_with(value = D_AGimg, location = ph_location(label = "my_name",
    #                                                 left = 0.3, top = 3.6, width = 2, height = 1.3)) %>%
    # ph_with(value = D_HYDimg, location = ph_location(label = "my_name",
    #                                                  left = 0.3, top = 5.4, width = 2, height = 1.3)) %>%
    # ph_with(value =  D_ECOimg, location = ph_location(label = "my_name",
    #                                                   left = 5, top = 2.6, width = 2, height = 1.3))%>%
    # ph_with(value =  D_SOCimg, location = ph_location(label = "my_name",
    #                                                   left = 5, top = 4.4, width = 2, height = 1.3))%>%
    # ## TEXT
    # ph_with(value =  fp_D_M, location = ph_location(label = "my_name",
    #                                                 left = 2.4, top = 1.8, width = 2.6, height = 0.69))%>%
    # ph_with(value =  fp_D_A, location = ph_location(label = "my_name",
    #                                                 left = 2.4, top = 3.6, width = 2.6, height = 1.17))%>% 
    # ph_with(value =  fp_D_H, location = ph_location(label = "my_name",
    #                                                 left = 2.4, top = 5.4, width = 2.6, height = 1.17))%>% 
    # ph_with(value =  fp_D_E, location = ph_location(label = "my_name",
    #                                                 left = 7.1, top = 2.6, width = 2.6, height = 1.3))%>% 
    # ph_with(value =  fp_D_S, location = ph_location(label = "my_name",
    #                                                 left = 7.1, top = 4.3, width = 2.78, height = 1.56))%>% 
    
    ph_with(value = D_fiveimg, location = ph_location(label = "my_name",
                                                      left = 0.65, top = 2.4, width = 8.8, height = 4.9))%>%
    ph_with(value = fp_D_t, location = ph_location(label = "my_name",
                                                      left = 0.5, top = 1, width = 9.1, height = 2))%>%
    ph_with(value = "Frazier et al. (2019)", location = ph_location_type(type = "dt"))%>%
    
    ph_with(value = p21, location = ph_location_type(type = "sldNum"))%>%
    
  #Slide 17 
  add_slide("Two Content","Office Theme") %>%
  ph_with(S15_TIT,       ph_location_type("title")) %>%
  ph_with(fp_SPI,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = p22, location = ph_location_type(type = "sldNum")) %>%
  #ph_with(value = SPIimg, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_12, location = ph_location(label = "my_name",
                                                left = 5.3, top = 5, width = 4.5, height = 1.4))%>%
  ph_with(value = SPIimg, location = ph_location(label = "my_name",
                                                left = 5, top = 2.5, width = 4.42, height = 2.72))%>%
  ph_with(my_pres, value = "Frazier et al. (2016); Lucas et al. (2022)", location = ph_location_type(type = "dt"))%>%

 
  #Slide 18
  add_slide("Two Content","Office Theme") %>%
  ph_with(S16_TIT,         ph_location_type("title")) %>%
  ph_with(fp_DHist,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = p23, location = ph_location_type(type = "sldNum")) %>%
  #ph_with(value = DHistimg, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_13, location = ph_location(label = "my_name",
                             left = 5.3, top = 5, width = 4.5, height = 1.4))%>%
  ph_with(value = DHistimg, location = ph_location(label = "my_name",
                             left = 5, top = 2.5, width = 4.42, height = 2.72))%>%
  ph_with(my_pres, value = "Frazier et al. (2016); Lucas et al. (2022)", location = ph_location_type(type = "dt"))%>% 
    
  #Slide 19
  add_slide("Title and Content","Office Theme") %>%
  ph_with(SP_TIT,ph_location_type("title",position_left = TRUE)) %>%
  #ph_with(value = "Lucas et al. (2022)", location = ph_location_type(type = "dt"))%>%
  ph_with(value = p24, location = ph_location_type(type = "sldNum"))%>%
  ph_with(value = SPI3_SMimg, location = ph_location(label = "my_name",
                                                     left = 0.4, top = 1.7, width = 4.2, height = 2.58)) %>%
  ph_with(value = SPI12_SMimg , location = ph_location(label = "my_name",
                                                       left = 5.1, top = 1.7, width = 4.2, height = 2.58)) %>%
  ph_with(value = FIG_14, location = ph_location(label = "my_name",
                                                 left = 0.5, top = 4.1, width = 4.2, height = 0.5)) %>%
  ph_with(value = FIG_15, location = ph_location(label = "my_name",
                                                 left = 5.2, top = 4.1, width = 4.2, height = 0.5)) %>%
  ph_with(value = SPI3v12, location = ph_location(label = "my_name",
                                                  left = 0.2, top = 3.9, width = 9.68, height = 2.5)) %>%
  ph_with(value = SPI3v12.1, location = ph_location(label = "my_name",
                                                  left = 0.4, top = 4.6, width = 9.68, height = 2.5)) %>%
  ph_with(value = SPI3v12.2, location = ph_location(label = "my_name",
                                                  left = 0.4, top = 5.1, width = 8.7, height = 2.5)) %>%
  ph_with(value = SPI3v12.3, location = ph_location(label = "my_name",
                                                  left = 0.4, top = 5.6, width = 8.7, height = 2.5)) %>%
    
   #Slide 20
   add_slide("Two Content","Office Theme") %>%
   ph_with(S18_TIT,  ph_location_type("title")) %>%
   ph_with(fp_Fire,        ph_location_type("body",position_right = FALSE)) %>%
   ph_with(value = p25, location = ph_location_type(type = "sldNum")) %>%
   # ph_with(value = FireMimg, ph_location_type("body",position_right = TRUE)) %>%
   ph_with(value = FireMimg, location = ph_location(label = "my_name",
                                                    left = 5.2, top = 1.7, height = 4.5, width = 4.5)) %>%
   ph_with(value = FIG_16, location = ph_location(label = "my_name",
                                                 left = 5.2, top = 6.1, width = 4.3, height = 1.4))%>%
   ph_with(my_pres, value = "Trauernicht, 2019; Frazier et al. (2022)", location = ph_location_type(type = "dt"))%>% 

   #Slide 21 #Part 4
   add_slide("Title and Content","Office Theme") %>%
   ph_with(P4_TIT, ph_location_type("title",position_left = TRUE)) %>%
   ph_with(fp_Part4,ph_location_type("body"))%>%
   ph_with(fp_Part4.2, location = ph_location(label = "my_name",
                                                    left = 0.45, top = 2.9, width = 6, height = 2)) %>%
   ph_with(fp_Part4.3, location = ph_location(label = "my_name",
                                               left = 0.45, top = 4, width = 6, height = 2)) %>%
    ph_with(fp_Part4.4, location = ph_location(label = "my_name",
                                               left = 0.45, top = 4.8, width = 9, height = 3)) %>%
   ph_with(value = p26, location = ph_location_type(type = "sldNum"))%>%
   ph_with(my_pres, value = "See Annex II", location = ph_location_type(type = "dt"))%>% 
   ph_with(value = Downimg, location = ph_location(label = "my_name",
                                                     left = 6.5, top = 3.4, width = 3, height = 2)) %>%
  # ph_with(value = M_Down, location = ph_location(label = "my_name",
  #                                                    left = 0.5, top = 4.5, width = 6.3, height = 2)) %>%
  
  #Slide 22
  add_slide("Two Content","Office Theme") %>%
  ph_with(S21_TIT,         ph_location_type("title")) %>%
  # ph_with(RF_Down4,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(RF_Down4, location = ph_location(label = "my_name",
                                                 left = 0.5, top = 1.7, width = 4, height = 5))%>%
  ph_with(value = p27, location = ph_location_type(type = "sldNum")) %>%
  #ph_with(value = RF40img, ph_location_type("body",position_right = TRUE))  %>%
  ph_with(value = FIG_18, location = ph_location(label = "my_name",
                                                 left = 5.2, top = 5.5, width = 4, height = 1.4))%>%
  ph_with(my_pres, value = "Timm et al., 2015; See Annex II", location = ph_location_type(type = "dt"))%>% 
  ph_with(value = RF40img, location = ph_location(label = "my_name",
                                                  left = 4.6, top = 2, width = 5.17, height = 3.51)) %>%
    
  #Slide 23
  add_slide("Two Content","Office Theme") %>%
  ph_with(S20_TIT,   ph_location_type("title")) %>%
  # ph_with(RF_Down,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(RF_Down, location = ph_location(label = "my_name",
                                           left = 0.5, top = 1.7, width = 4, height = 5))%>%
  ph_with(value = p28, location = ph_location_type(type = "sldNum")) %>%
  #ph_with(value = RF10085img, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_17, location = ph_location(label = "my_name",
                                                 left = 5.2, top = 5.5, width = 4, height = 1.4))%>%
  ph_with(my_pres, value = "Timm et al., 2015; Zhang et al., 2016; See Annex II", location = ph_location_type(type = "dt"))%>% 
  ph_with(value = RF10085img, location = ph_location(label = "my_name",
                                                   left = 4.6, top = 2, width = 5.3, height = 3.6)) %>%
    
    #Slide 24
    add_slide("Two Content","Office Theme") %>%
    ph_with(S23_TIT,         ph_location_type("title")) %>%
    # ph_with(TA_Down4,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(TA_Down4, location = ph_location(label = "my_name",
                                            left = 0.5, top = 1.7, width = 4, height = 5))%>%
    ph_with(value = p29, location = ph_location_type(type = "sldNum")) %>%
    #ph_with(value = TA40img, ph_location_type("body",position_right = TRUE)) %>%
    ph_with(value = FIG_20, location = ph_location(label = "my_name",
                                                   left = 5.2, top = 5.5, width = 4, height = 1.4))%>%
    ph_with(my_pres, value = "Timm 2017; See Annex II", location = ph_location_type(type = "dt"))%>% 
    ph_with(value = TA40img, location = ph_location(label = "my_name",
                                                    left = 4.6, top = 2, width = 5.17, height = 3.51)) %>%
    
    #Slide 25
    add_slide("Two Content","Office Theme") %>%
    ph_with(S22_TIT,       ph_location_type("title")) %>%
    # ph_with(TA_Down,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(TA_Down, location = ph_location(label = "my_name",
                                            left = 0.5, top = 1.7, width = 4, height = 5))%>%
    ph_with(value = p30, location = ph_location_type(type = "sldNum")) %>%
    #ph_with(value = TA100img, ph_location_type("body",position_right = TRUE)) %>%
    ph_with(value = FIG_19, location = ph_location(label = "my_name",
                                                   left = 5.2, top = 5.5, width = 4, height = 1.4))%>%
    ph_with(my_pres, value = "Timm 2017; Zhang et al., 2016; See Annex II", location = ph_location_type(type = "dt"))%>% 
    ph_with(value = TA100img, location = ph_location(label = "my_name",
                                                       left = 4.6, top = 2, width = 5.3, height = 3.6)) %>%
    
    #Slide 26 #Summary
    add_slide("Title and Content","Office Theme") %>%
    ph_with(S24_TIT,   ph_location_type("title",position_left = TRUE)) %>%
    ph_with(fp_Summary,ph_location_type("body"))%>%
    ph_with(value = p31, location = ph_location_type(type = "sldNum"))%>%

    #Slide 27 #External Resources
    add_slide("Title and Content","Office Theme") %>%
    ph_with(S25_TIT, ph_location_type("title",position_left = TRUE)) %>%
    # ph_with(fp_RES, ph_location_type("body"))%>%
    ph_with(value =  fp_RES , location = ph_location(label = "my_name",
                                                 left = 0.5, top = 0.4, width = 10, height = 8))%>%
    ph_with(value = p32, location = ph_location_type(type = "sldNum"))%>%
    # ph_with(value =  DPlanimg , location = ph_location(label = "my_name",
    #                                              left = 7.3, top = 1.9, width = 2.3, height = 2.83))%>%

    #Slide 28 #Acknowledgements
    add_slide("Title and Content","Office Theme") %>%
    ph_with(S26_TIT, ph_location_type("title",position_left = TRUE)) %>%
    ph_with(fp_Ack, location = ph_location(label = "myname",
                                                 left = 0.5, top = 0.08, width = 9, height = 5))%>%
    ph_with(value = p33, location = ph_location_type(type = "sldNum"))%>%
    ph_with(value = EWCimg, location = ph_location(label = "my_name",
                                                   left = 0.09, top = 3.53, width = 1.2, height = 1.2)) %>%
    ph_with(value = CASCimg, location = ph_location(label = "my_name",
                                                    left = 1.3, top = 3.7, width = 0.9, height = 0.9)) %>%
    ph_with(value = FSimg, location = ph_location(label = "my_name",
                                                  left = 2.25, top = 3.7, width = 0.95, height = 0.9)) %>%
    ph_with(value = PIDPimg, location = ph_location(label = "my_name",
                                                    left = 3.22, top = 3.7, width = 1.4, height = 0.8)) %>%
    ph_with(value = NIDISimg, location = ph_location(label = "my_name",
                                                     left = 4.64, top = 3.7, width = 0.8, height = 0.8)) %>%
    ph_with(value = WRRCimg, location = ph_location(label = "my_name",
                                                    left = 5.59, top = 3.7, width = 0.9, height = 0.9)) %>%
    ph_with(value = NOAAimg, location = ph_location(label = "my_name",
                                                    left = 6.6, top = 3.7, width = 0.9, height = 0.9)) %>%
    ph_with(value = UHimg, location = ph_location(label = "my_name",
                                                  left = 7.69, top = 3.7, width = 0.9, height = 0.9)) %>%
    ph_with(value = CLimg, location = ph_location(label = "my_name",
                                                  left = 8.75, top = 3.7, width = 0.9, height = 0.9)) %>%
    ph_with(value = CITA, location = ph_location(label = "my_name",
                                                   left = 0.85, top = 5.1, width = 8, height = 1.4))%>%
    ph_with(value = Draft, location = ph_location(label = "my_name",
                                                 left = 0.85, top = 6.6, width = 8, height = 0.75))%>%
    ph_with(value = PDKE_S, location = ph_location(label = "my_name",
                                                left = 0, top = 0, width = 1.5, height = 1.5))%>%
  
    #Slide 29  ########## References 
    
    add_slide("Title Only","Office Theme") %>%
    ph_with(S27_TIT,  ph_location_type("title",position_left = TRUE)) %>%
    ph_with(value = Worksimg, location = ph_location(label = "my_name",
                                                  left = 1, top = 1.5, width = 7.9, height = 5.3)) %>%
    ph_with(value = p34, location = ph_location_type(type = "sldNum"))%>%
    
    #Slide 30  ############## Annex I
    
    add_slide("Two Content","Office Theme") %>%
    ph_with("Annex I: 100+ Year Rainfall ",ph_location_type("title")) %>%
    ph_with(fp_RF100,        ph_location_type("body",position_right = FALSE)) %>%
    #ph_with(value = p33, location = ph_location_type(type = "sldNum")) %>%
    ph_with(value = SCAimg, ph_location_type("body",position_right = TRUE)) %>%
    ph_with(value = FIG_A1, location = ph_location(label = "my_name",
                                                   left = 5, top = 6.3, width = 4.8, height = 1.4))%>%
    ph_with(my_pres, value = "Frazier et al., 2016; Lucas et al., (2022)", location = ph_location_type(type = "dt"))%>% 
    
    
    #Slide 31  ############## Annex II
    
    add_slide("Title and Content","Office Theme") %>%
    ph_with("Annex II: Climate Downscaling in Hawaii ",ph_location_type("title")) %>%
    ph_with(fp_DownA,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(value = p36, location = ph_location_type(type = "sldNum")) %>%
    ph_with(value = Downimg, location = ph_location(label = "my_name",
                                                    left = 3.5, top = 5, width = 3, height = 2)) %>%
  
    #Slide 32 ############### Annex III

    add_slide("Title and Content","Office Theme") %>%
    ph_with("Annex III: Drought Events (1920 - 2022)",ph_location_type("title")) %>%
    ph_with(value = p37, location = ph_location_type(type = "sldNum"))%>% 
    ph_with(value = DrotT_T1, ph_location_type("body")) %>%
    ph_with(value = TAB3, location = ph_location(label = "my_name",
                                                 left = 7.85, top = 2.5, width = 2, height = 4)) %>%
    
    {if(Drot_ct>12) add_slide(.,"Title and Content", "Office Theme") %>% 
        ph_with("Annex III: Drought Events (1920 - 2022)",ph_location_type("title")) %>%
        ph_with(value = (p37+1), location = ph_location_type(type = "sldNum"))%>%
        ph_with(value = DrotT_T2, ph_location_type("body"))%>%
        ph_with(value = TAB3,
                location = ph_location(label = "my_name",
                                       left = 7.85, top = 2.5, width = 2, height = 4)) else .} %>%
    
    {if(Drot_ct>24) add_slide(.,"Title and Content", "Office Theme") %>% 
        ph_with("Annex III: Drought Events (1920 - 2022)",ph_location_type("title")) %>%
        ph_with(value = (p37+2), location = ph_location_type(type = "sldNum"))%>%
        ph_with(value = DrotT_T3, ph_location_type("body"))%>%
        ph_with(value = TAB3,
                location = ph_location(label = "my_name",
                                       left = 7.85, top = 2.5, width = 2, height = 4)) else .} %>%

  #Slide 33 (or whatever the last slide is)  ############## PDKE 
  
  add_slide("Title and Content","Office Theme") %>%
  ph_with(value = PDKE_S, location = ph_location(label = "my_name",
                                                 left = 2, top = 1, width = 6, height = 6))%>%

  
  print(mypowerpoint, target = paste0(P_FOLDER,NameF,"_CCVD_Portfolio_v",ver,".pptx"))

  
  
  # ### Save to powerpoint
  # 
  # office_shot <- function( file, wd = getwd() ){
  #   cmd_ <- sprintf(
  #     "C:/Program Files/LibreOffice/program/soffice.exe",
  #     wd, file )
  #   system(cmd_)
  #   
  #   pdf_file <- gsub("\\.(docx|pptx)$", ".pdf", basename(file))
  #   pdf_file
  # }
  # office_shot(file = paste0(P_FOLDER,NameF,"_CCVD_Portfolio_v",ver,".pptx"))
  
  
  
#   #################################################################################################################
#   #Slides not used
#   
#     #Slide 29  ############## Annex I
#     
#     add_slide("Two Content","Office Theme") %>%
#     ph_with("Annex III: Roads and Trails ",ph_location_type("title")) %>%
#     ph_with(Roadimg,        ph_location_type("body",position_right = FALSE)) %>%
#     #ph_with(value = "29", location = ph_location_type(type = "sldNum")) %>%
#     ph_with(value = Trailimg, ph_location_type("body",position_right = TRUE)) %>%
#     ph_with(value = FIG_A2, location = ph_location(label = "my_name",
#                                                    left = 1.5, top = 6.1, width = 4, height = 1.4))%>%
#     ph_with(value = FIG_A3, location = ph_location(label = "my_name",
#                                                    left = 5.7, top = 6.1, width = 4, height = 1.4))%>%
#  
#   ###############   Slide xx  Relative Humidity
#   
#   S10_TIT<- block_list(
#     fpar(ftext("Average Monthly Relative Humidity", FTXTT),fp_p = fp_par(text.align = "center")))
# 
# RHx <- max(RH)
# RHxm <- match(RHx,RH)
# RHn <- min(RH)
# RHnm <- match(RHn,RH)
# RHD <- RHx - RHn
# 
# RHA <- paste0("Relative humidity (RH) is a measure of the water vapor content in the air relative to how much the air can hold (which is a function of temperature).",
#               " RH can play a critical role in plant growth.  ", 
#               "There is a ",RHD,"% difference in RH over the year at ", NameF, ", where highest monthly RH occurs in ", MONLIST[RHxm], 
#               " (",RHx, "%) and the lowest occurs in ", MONLIST[RHnm]," (",RHn,"%).")
# 
# fp_RH <- fpar(ftext(RHA, fp_Tx))
# 
# FIG_7 <- block_list(
#   fpar(ftext(paste0("Figure 7. Average monthly relative humidity ",NameF,
#                     " with area average shown in heading of each plot."), fp_Fig)))
# 
# #RH Figure
# RHfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," RH12.png")
# RHimg <- external_img(src = RHfile, height = 3,width = 3) 
# 
# 
# AUT <- block_list(
#   fpar(ftext(" ", fp_NM),fp_p = fp_par(text.align = "center")),
#   fpar(ftext("Authors", fp_NM),fp_p = fp_par(text.align = "center")),
#   fpar(ftext(" ", fp_NMS),fp_p = fp_par(text.align = "center")),
#   fpar(ftext("Ryan J. Longman, East-West Center", fp_NM2),fp_p = fp_par(text.align = "center")),
#   fpar(ftext(" ", fp_NMS),fp_p = fp_par(text.align = "center")),
#   fpar(ftext("Abby G. Frazier, East-West Center", fp_NM2),fp_p = fp_par(text.align = "center")),
#   fpar(ftext(" ", fp_NMS),fp_p = fp_par(text.align = "center")),
#   fpar(ftext("Christian P. Giardina, USDA Forest Service ", fp_NM2),fp_p = fp_par(text.align = "center")))
# 
# 
# 
# 
# 
