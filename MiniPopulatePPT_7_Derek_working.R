
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

setwd("E:/PDKE/CCVD/MINI_Phase2/")
#Path = "C:/Users/BUNNY1/Documents/Work From Home/PAPERS/Beef_Production/Specific Ranch/"
I_FOLDER <- "E:/PDKE/CCVD/IMAGE/"        # Folder with images 
R_FOLDER <- "E:/PDKE/CCVD/MINI_Phase2/"  # File with your site specific folders 
P_FOLDER <- "E:/PDKE/CCVD/MINI_PPT/"     # Output folder 
RANL <- list.files(R_FOLDER)
NF <- length(RANL)

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
RFUnit3 = " inches"
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

#media
D_AGfileM <- paste0(I_FOLDER,"Media1_SPI.MP4")
D_AGimgM <- external_img(src = D_AGfileM , height = 1, width = 1)



######
RANL
# LOOP 
f<-2





#for(f in 1:NF) {
#for(f in 79:82) {

#create path to folder 
RAN_F <- paste0(R_FOLDER,RANL[f],"/")
# list the output from Code 1 
RANL2 <- list.files(RAN_F)
#Unit Name
NameF <- basename(RAN_F)

## Read in Mean Climate File 
## Read in CSV Files 
CLIM <- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF,"MEAN Climate.csv"),sep=",")

  #Short Name
  SNameF <- CLIM[1,15]
  
  #Font styles 
  #Suggetions Type_Color_Size  ex B_DR_40
  
  fp_BR <- fp_text(bold = TRUE, color = "darkred", font.size = 40)
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
  fp_Fig <- fp_text(font.size = 12, color = "black")
  fp_Fig2 <- fp_text(font.size = 10, color = "black")
  fp_Fig3 <- fp_text(font.size = 14, color = "black")
  fp_Fig4<- fp_text(bold = TRUE, color = "darkred",font.size = 15,underlined = TRUE)
  
  FTXTT <-  fp_text(bold = TRUE,color = "darkblue", font.size = 40)
  FTXTT2 <- fp_text(italic = TRUE, color = "black", font.size = 40)
  FTXTT3 <- fp_text(italic = TRUE, color = "darkred", font.size = 40)
  FTXTT4 <- fp_text(color = "black", font.size = 38)
  
  ############# Slide 1 TITLE 
    
  #Slide one Title 
  #Standardize S1_TIT
  TIT <- fpar(ftext("Climate Change, Climate Variability, & Drought Portfolio", fp_BR3),fp_p = fp_par(text.align = "center") )
  SUB <- fpar(ftext(NameF, fp_BR),fp_p = fp_par(text.align = "center") )
 

 
################ Slide 2 PDKE
  
  #S2_TIT<- block_list(
    #fpar(ftext("Pacific Drought Knowledge Exchange", FTXTT),fp_p = fp_par(text.align = "center")))
  
  HDKE<- paste0("Climate change, climate variability, and drought (CCVD) will exert a growing impact on Hawaii's landscapes, watersheds, ",
              "and nearshore areas in the future. Similar impacts will be felt across much of the Pacific as well. While managers are tasked with utilizing the best available science, they are often ",
              "unaware of what datasets or products are available as there is no centralized, drought-focused information clearinghouse ",
              "or mechanism to engage with scientists in research, planning, product development, and knowledge co-production. The Pacific ",
              "Drought Knowledge Exchange (PDKE) funded by the Pacific Islands Climate Adaptation Science Center, focuses on ",
              "facilitating knowledge exchange between the research community and resource managers and user groups, thereby expanding the utility ",
              "of climate and drought-related products for resource managers.") 
  
  #Font Type consider moving up to top 
fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 20) 
#May or may not be used
fp_TxB <- fp_text(italic = TRUE, color = "black", font.size = 20) 

# Text for slide as a variable 
fp_HDKE <- fpar(ftext(HDKE, fp_TxB))


###############   Slide 3 CCVD
  
  #not sure why I use "block list" check this. 

  S3_TIT<- block_list(
  fpar(ftext("CCVD Portfolio", FTXTT),fp_p = fp_par(text.align = "center")))

  CCVD<- paste0("The climate change, climate variability, and drought (CCVD) portfolio is a comprehensive synthesis ",
              "of climate and drought information developed specifically for ", NameF, " (",SNameF,"). This portfolio is designed to ",
              "provide both research scientists and land managers with relevant climate and drought information ",
              "needed to inform land management and guide future research and extension.") 

  fp_CCVD <- fpar(ftext(CCVD, fp_Tx))

#Map Figure 
MAPfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Map.png")
MAPimg <- external_img(src = MAPfile, height = 1,width = 1) 


FIG_1 <- block_list(fpar(ftext(paste0("Figure 1. Map of ",SNameF,"."), fp_Fig)))



################ Slide 4 PART 1

P1_TIT<- block_list(
       fpar(ftext("Part 1: Climate Characteristics", FTXTT)),
       fpar(ftext(SNameF, FTXTT3),fp_p = fp_par(text.align = "center")))

PART1<- paste0("In developing this Portfolio, we relied on several gridded climate products available for the State of Hawaii. Both mean annual and mean monthly estimates of rainfall were ",
              "are obtained from the Rainfall Atlas of Hawaii (http://rainfall.geography.hawaii.edu/). Gridded estimates of other climate variables are ",
              "obtained from the Climate of Hawaii (http://climate.geography.hawaii.edu/). ",
              "We retrieved all the data points that fell within the boundaries of ",SNameF," from our 250 m resolution state-wide maps to support the presented analyses.")
    
fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 20) 
fp_PART1 <- fpar(ftext(PART1, fp_Tx))

#NameF_I_Tx <- fp_text(italic = TRUE, color = "darkblue", font.size = 20) 
#NameF_I <- fp_text(italic = TRUE, color = "darkblue", font.size = 20) 



###############   Slide 5 Elevation

S5_TIT<- block_list(
  fpar(ftext("Elevation", FTXTT),fp_p = fp_par(text.align = "center")))

#Pulling from table that was read in
El_Mean <- round(CLIM[1,2],0)
El_Max<-  round(CLIM[2,2],0)
El_Min <-  round(CLIM[3,2],0)
#calculat Elevation Difference
EL_Dif <- El_Max - El_Min
ISL<- CLIM[1,14]


ELEV_T<- paste0(SNameF," is located on the Island of ", ISL, " and covers a vertical elevation range of ", EL_Dif, ELUnit, ". In Hawaii, climate gradients can change ",
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
  fpar(ftext(paste0("Figure 2. Elevation for the Island of ",ISL," with ",SNameF," outlined in black. ",
                    "The maps shown in the following slides will be for the ",SNameF," area only."), fp_Fig)))

#MRead in elevation figure 
    Elevfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," ELMap.png")
    Elevimg <- external_img(src = Elevfile, height = 3,width = 3) 

###############   Slide 6 Mean Climate 
    
    S6_TIT<- block_list(
      fpar(ftext("Average Annual Climate Characteristics", FTXTT),fp_p = fp_par(text.align = "center")))
    
    #GET RANGES
    RFR   <- paste("  Spatial Range: ", CLIM[3,3], " to ", CLIM[2,3], RFUnit) 
    TavgR <- paste("  Spatial Range: ", CLIM[3,4], " to ", CLIM[2,4], TUnit2) 
    TmaxR <- paste("  Spatial Range: ", CLIM[3,5], " to ", CLIM[2,5], TUnit2) 
    TminR <- paste("  Spatial Range: ", CLIM[3,6], " to ", CLIM[2,6], TUnit2) 
    RHR   <- paste("  Spatial Range: ", CLIM[3,7], " to ", CLIM[2,7], "%") 
    SMR   <- paste("  Spatial Range: ", CLIM[3,8], " to ", CLIM[2,8], " Ratio") 
    KDR   <- paste("  Spatial Range: ", CLIM[3,9], " to ", CLIM[2,9], "W/m2") 
    ETR   <- paste("  Spatial Range: ", CLIM[3,10], " to ", CLIM[2,10], RFUnit) 
    CFR   <- paste("  Spatial Range: ", CLIM[3,11], " to ", CLIM[2,11], " Ratio") 
    

    fp_NM31 <- fp_text(font.size = 17, color = "darkgreen")
    fp_NM33 <- fp_text(font.size = 15, color = "black")
    
    
M_Climate <-  block_list(
  fpar(ftext("RF- Rainfall", fp_NM31)),
  fpar(ftext(RFR , fp_NM33)),
  fpar(ftext("RH- Relative humidity", fp_NM31)),
  fpar(ftext(RHR , fp_NM33)),
  fpar(ftext("SM- Soil moisture", fp_NM31)),
  fpar(ftext(SMR , fp_NM33)),
  fpar(ftext("Mean TA- Average temperature", fp_NM31)),
  fpar(ftext(TavgR , fp_NM33)),
  fpar(ftext("Min TA- Average minimum temperature", fp_NM31)),
  fpar(ftext(TminR  , fp_NM33)),
  fpar(ftext("Max TA- Average maximum temperature", fp_NM31)),
  fpar(ftext(TmaxR , fp_NM33)),
  fpar(ftext("CF- Cloud cover", fp_NM31)),
  fpar(ftext(CFR , fp_NM33)),
  fpar(ftext("S- Solar radiation", fp_NM31)),
  fpar(ftext(KDR , fp_NM33)),
  fpar(ftext("ET- Evapotranspiration", fp_NM31)),
  fpar(ftext(ETR , fp_NM33)))
  
FIG_3 <- block_list(
  fpar(ftext(paste0("Figure 3. Mean annual climate of ",SNameF,
                    " with area average shown in heading of each plot."), fp_Fig)))

#Mean CLIM Figure 
Climfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climate.png")
Climimg <- external_img(src = Climfile, height = 3,width = 3) 




###############   Slide 7 Climograph 

S7_TIT<- block_list(
  fpar(ftext("Average Monthly Rainfall", FTXTT),fp_p = fp_par(text.align = "center")),
  fpar(ftext("and Temperature", FTXTT),fp_p = fp_par(text.align = "center")))

CLIMA <- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Annual Climate.csv"),sep=",")

RFT <- round(CLIMA[1,2:13],0)
MinT<- round(CLIMA[2,2:13],0)
TAT <- round(CLIMA[3,2:13],0) # Mean Temperature 
MaxT <- round(CLIMA[4,2:13],0)
RH <- round(CLIMA[5,2:13],0)

MONLIST <- c("January","February","March","April","May","June","July","August","September","October","November","December")

RFx <- max(RFT)
RFxm <- match(RFx,RFT)
RFn <- min(RFT)
RFnm <- match(RFn,RFT)

TAx <- max(TAT)
TAxm <- match(TAx,TAT)
TAn <- min(TAT)
TAnm <- match(TAn,TAT)
Tdif <- abs(TAn-TAx)

CLIMO <- paste0("Average monthly rainfall and temperature patterns vary over the course of the year at ", SNameF, ". The highest monthly rainfall is received in ", MONLIST[RFxm], 
                " (",RFx, RFUnit,") and the lowest monthly rainfall is received in ", MONLIST[RFnm]," (",RFn, RFUnit,"). There is a ",Tdif,TUnit ," annual variation in temperature with the ",
                "warmest month occurring in ", MONLIST[TAxm]," (",TAx, TUnit2,") and the coolest month occurring in ", MONLIST[TAnm]," (",TAn, TUnit2,").")

fp_CLIMO <- fpar(ftext(CLIMO, fp_Tx))

FIG_4 <- block_list(
  fpar(ftext(paste0("Figure 4. Mean monthly rainfall and temperature at ",SNameF,"."), fp_Fig)))
  
#Mean CLIM Figure 
Climofile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Climograph.png")
Climoimg <- external_img(src = Climofile, height = 3,width = 3) 



###############   Slide 8 TA And RF

S8_TIT<- block_list(
  fpar(ftext("Average Monthly Climate", FTXTT),fp_p = fp_par(text.align = "center")))

TAfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," TA12.png")
TAimg <- external_img(src = TAfile , height = 3,width = 3) 

RFFfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," RF12.png")
RFFimg <- external_img(src = RFFfile , height = 3,width = 3) 

FIG_5 <- block_list(
  fpar(ftext(paste0("Figure 5. Mean monthly rainfall  ",SNameF,
                    " with area average shown in heading of each plot."), fp_Fig)))

FIG_6 <- block_list(
  fpar(ftext(paste0("Figure 6. Mean monthly temperature at ",SNameF,
                    " with area average shown in heading of each plot."), fp_Fig)))


###############   Slide 9 Seasonal Rainfall 

S9_TIT<- block_list(
  fpar(ftext("Average Seasonal Rainfall", FTXTT),fp_p = fp_par(text.align = "center")))
 

#Season 
DSeaMRF  <- round(CLIM[1,12],0)
DSeaMRFx <- round(CLIM[2,12],0)
DSeaMRFn <- round(CLIM[3,12],0)
WSeaMRF  <- round(CLIM[1,13],0)
WSeaMRFx <- round(CLIM[2,13],0)
WSeaMRFn <- round(CLIM[3,13],0)

#permonth 
DSeaMRF_M<- round((CLIM[1,12]/6),1)
DSeaMRF_Mx<- round((CLIM[2,12]/6),1)
DSeaMRF_Mn<- round((CLIM[3,12]/6),1)
WSeaMRF_M <- round((CLIM[1,13]/6),1)
WSeaMRF_Mx <- round((CLIM[2,13]/6),1)
WSeaMRF_Mn <- round((CLIM[3,13]/6),1)

#Per Month Range
D_RF_R_M  <- paste0(DSeaMRF_Mn, " to ", DSeaMRF_Mx) 
W_RF_R_M  <- paste0(WSeaMRF_Mn, " to ", WSeaMRF_Mx)

#Season Range 
D_RF_R  <- paste0(DSeaMRFn, " to ", DSeaMRFx) 
W_RF_R  <- paste0(WSeaMRFn, " to ", WSeaMRFx) 

WetM <- WSeaMRF_M
DryM <- DSeaMRF_M
  
SeaD <- round((WetM - DryM),1)

SEA <- paste0("Hawaii has two distinct 6-month seasons of rainfall. The Dry season runs from October to May and the Wet season runs from November to April. ",
                "At ",SNameF,", dry season rainfall averages ",DSeaMRF_M, RFUnit3," per month (", DSeaMRF, RFUnit3," total) and wet season rainfall averages ",
                 WSeaMRF_M, RFUnit3," per month (", WSeaMRF, RFUnit3," total). ",
                "Across the entire unit, rainfall ranges from ", D_RF_R_M, RFUnit3, " per month (",D_RF_R," total) during the dry season and ",W_RF_R_M, " during the wet season (",W_RF_R," total).")

fp_Tx19 <- fp_text(italic = TRUE, color = "black", font.size = 19) 
fp_SEA <- fpar(ftext(SEA, fp_Tx19))

FIG_7 <- block_list(
  fpar(ftext(paste0("Figure 7. Dry and Wet season total rainfall at ",SNameF,"."), fp_Fig)))

#Mean CLIM Figure 
SEAfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," SeaRF.png")
SEAimg <- external_img(src = SEAfile, height = 3,width = 3) 



###############   Slide 10  Average Climate Table

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
  fpar(ftext(paste0("Table 1. Average monthly climate variables characteristics at ",NameF, ". Where, RF is rainfall; Min TA is average minimum air temperature",
                    " Mean TA is average air temperature; Max TA is average maximum air temperature; RH is relative humidity; CF is cloud frequency; ET is ",
                    "evapotranspiration; SM is soil moisture; S is shortwave downward radiation: ANN, is annual total for rainfall and annual average for all other variables."), fp_Fig)))


################ Slide 11 PART 2

  P2_TIT<- block_list(
  fpar(ftext("Part 2: Inter-Annual Rainfall", FTXTT),fp_p = fp_par(text.align = "center")),
  fpar(ftext(SNameF, FTXTT3),fp_p = fp_par(text.align = "center")))

  InterAn <- paste0("Rainfall in Hawaii can vary greatly from year-to-year due to natural modes of climate variability such as the El Ni?o-Southern Oscillation ",
                  "(ENSO). ENSO can be explained as an interaction between the atmosphere and the ocean in the tropical Pacific that results in a somewhat periodic ",
                  "variation between below-normal and above-normal sea surface temperatures. In Hawaii, wet season rainfall is typically low during the ",
                  "warm (El Ni?o) phase of ENSO and high during the neutral and Cool (La Ni?a) Phases. This pattern is reversed in the dry season although it is not as pronounced as in the wet season. ",
                  "ENSO is the dominant mode of climate variability in Hawaii.")
 
            
fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 17) 
fp_InterAn <- fpar(ftext(InterAn, fp_Tx))

ENSOfile <- paste0(I_FOLDER,"ENSO_New.png")
ENSO2img <- external_img(src = ENSOfile, height = 1,width = 1) 


################ Slide 12 ENSO vs Rainfall 

MEIRF<- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," MEI_A.csv"),sep=",")

EN_TIT<- block_list(
  fpar(ftext("Wet Season Rainfall", FTXTT),fp_p = fp_par(text.align = "center")),
  fpar(ftext("vs El Ni?o/La Ni?a", FTXTT),fp_p = fp_par(text.align = "center")))

SEL_RFW <- round(MEIRF[1,2],0)
#  WetM # Wet-Season Avg

DifWRF <- round(((WetM - SEL_RFW) / WetM) * 100,0)
#DifWRF <- 100-DifWRF

MEIW <- paste0("93-years of monthly wet season rainfall (1920-2012) are compared with the Multivariate ENSO Index (MEI) to determine how rainfall is influenced ",
                  "by five different ENSO phases. During the strong El Ni?o phase, average monthly wet season rainfall (",SEL_RFW, RFUnit,"/month) is ", DifWRF,"% drier than the long-term average (",WSeaMRF_M, RFUnit,"/month).")

fp_TxW <- fp_text(italic = TRUE, color = "black", font.size = 18) 
fp_MEIW  <- fpar(ftext(MEIW, fp_TxW))

ENSOKfile <- paste0(I_FOLDER,"ENSOKEY.png")
ENSOKimg <- external_img(src = ENSOKfile, height = 5,width = 2) 


MEIWfile <- paste0(R_FOLDER,"/",NameF,"/",NameF,"MEI_WET.png")
MEIWimg <- external_img(src = MEIWfile, height = 2.5,width = 2.5) 

FIG_8 <- block_list(
  fpar(ftext(paste0("Figure 8. Boxplot of wet season monthly rainfall ",
                    "grouped by ENSO category. ",SNameF,"."), fp_Fig)))

TAB2 <- block_list(
  fpar(ftext("Table 2. ENSO phases and abbreviations corresponding to Figure 8.", fp_Fig)))



#################### Slide 13 ###############################

S13_TIT<- block_list(
  fpar(ftext("Long-Term Trends in Rainfall", FTXTT),fp_p = fp_par(text.align = "center")))



Trend <- paste0("Linear trends in annual, wet season, and dry season rainfall are calculated over a 100-year record at ",SNameF," for six different periods in the record. ",
                "Each trend period has a unique start year, but they all end in 2019. When the p-value is less than 0.05, the trend is determined to be statistically significant.")

fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 20) 
fp_Trend <- fpar(ftext(Trend, fp_Tx))

Trendfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," RF_Trend.png")
Trendimg <- external_img(src = Trendfile, height = 2,width = 2) 

FIG_9 <- block_list(
  fpar(ftext(paste0("Figure 9. 100-year (1920-2019) rainfall time series at " ,SNameF," with ",
                    "linear trends calculated over six unique periods in the record; Trend is ", 
                    "the slope (", RFUnit3," per year); R2 is the strength of the trend; ", 
                     "p is a measure of the statistical significance."), fp_Fig2)))



################ Slide 14 PART 3 Drought and Fire Occurrence ###############
#Today, 
#increased wildland fire activity, 
#and damage to terrestrial and aquatic habitats - all of which can contribute to substantial economic losses# 

P3_TIT<- block_list(
  fpar(ftext("Part 3: Drought and Fire History", FTXTT)),
  fpar(ftext(SNameF, FTXTT3),fp_p = fp_par(text.align = "center")))

Part3 <- paste0("Drought is a prominent feature of the climate system in Hawaii and can cause severe impacts across multiple sectors. ",
                  "Droughts in Hawaii often result in reduced crop yields, loss of livestock, drying of streams and reservoirs, depletion ",
                   "of groundwater, and increased wildland fire activity. These impacts can cause substantial economic ",
                  "losses as well as long-term damage to terrestrial and aquatic habitats.")

fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 20) 
fp_Part3 <- fpar(ftext(Part3, fp_Tx))



####################  SLID 15 FIVE TYPES OF DROUght######################################


TIT_D <- fpar(ftext("5-Types of Drought", FTXTT),fp_p = fp_par(text.align = "center") )
fp_D_M <-  block_list(
  fpar(ftext( "Meteorological Drought.", fp_Fig4)),
  fpar(ftext( "Defined as a lack of rainfall" , fp_Fig3)))

fp_D_A <-  block_list(
  fpar(ftext( "Agriculture Drought", fp_Fig4)),
  fpar(ftext( "Refers to a period of declining soil moisture and subsequent crop failure", fp_Fig3)))

fp_D_H <-  block_list( 
  fpar(ftext( "Hydrological Drought", fp_Fig4)),
  fpar(ftext( "Expressed as decreased streamflow and sub-surface water storage", fp_Fig3)))

fp_D_E<-  block_list(
  fpar(ftext( "Ecological Drought", fp_Fig4)),
  fpar(ftext( "Includes any impacts to ecosystems including an increase in wildfire occurrence" , fp_Fig3)))

fp_D_S<-  block_list(   
  fpar(ftext( "Socio-Economic Drought", fp_Fig4)),
  fpar(ftext( "Includes impacts to social and economic systems, including increased costs or revenue losses, and impacts on public health and safety", fp_Fig3)))




#################### Slide 16 ###############################

S15_TIT<- block_list(
  fpar(ftext("Identifying Droughts Using the", FTXTT),fp_p = fp_par(text.align = "center")),
  fpar(ftext("Standard Precipitation Index", FTXTT),fp_p = fp_par(text.align = "center")))


SPI <- paste0("The Standardized Precipitation Index (SPI) is one of the most widely used drought indices. SPI compares rainfall with ",
              "its multi-year average, and because droughts are generally defined relative to the local normal, this standardized index allows wet and dry climates ",
              "to be represented on a common scale. Here, Here, 100-years of monthly rainfall are used to used to calculate SPI-12, which compares how ",
              "a 12-month period compares with all 12-month periods in the record. SPI-12 is a good measure of sustained droughts that affect hydrological processes ",
              "at ",SNameF, ".")

  fp_Txa <- fp_text(italic = TRUE, color = "black", font.size = 17) 
  fp_SPI <- fpar(ftext(SPI, fp_Txa))

  SPIfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," SPI.png")
  SPIimg <- external_img(src = SPIfile, height = 2,width = 2) 

  FIG_10 <- block_list(
   fpar(ftext(paste0("Figure 10. 100-year (1920-2019) SPI-12 time series at ", 
                      SNameF ," positive SPI (blue) indicate wet periods, negative ", 
                     "SPI (red) indicate dry periods."), fp_Fig)))


  
  #################### Slide 17 ###############################
  
  S16_TIT<- block_list(
    fpar(ftext("A 100-Year History of Drought", FTXTT),fp_p = fp_par(text.align = "center")))
   

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
  
  DROT1
  DROT2
  DROT3
  
  DrotT_T1 <- flextable(DROT1)
  DrotT_T1 <- bg(DrotT_T1, bg = "orange", part = "header")
  DrotT_T1 <- autofit(DrotT_T1)
  DrotT_T1
  
  DrotT_T2 <- flextable(DROT2)
  DrotT_T2 <- bg(DrotT_T2, bg = "orange", part = "header")
  DrotT_T2 <- autofit(DrotT_T2)
  DrotT_T2
  
  DrotT_T3 <- flextable(DROT3)
  DrotT_T3 <- bg(DrotT_T3, bg = "orange", part = "header")
  DrotT_T3 <- autofit(DrotT_T3)
  DrotT_T3
  
  # # make conditional statement to get last slide number (will be used on last slide)
  # if(Drot_ct>12) {ls <- 33}
  # if(Drot_ct>24) {ls <- 34}
  
  DHist1 <- paste0("Negative SPI values (dry periods) are inverted to show a complete drought timeseries at ", SNameF, ". Dashed lines and corresponding color coding indicates instances ",
                 "of Moderate (SPI > 1), Severe (SPI > 1.5), and Extreme (SPI > 2) drought.") 
                 
  DHist2 <- paste0("A total of ",CNTDRT," Droughts were observed over the 100-year record ",
                   "with total of ", SoG, " drought events of severe strength or greater. The longest drought lasted for a total of ",LongD," months (see Annex III).")
 
  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 18)
  
  fp_DHist<- block_list(
             fpar(ftext(DHist1, fp_Tx)),
             fpar(ftext(" ", fp_Tx)),
             fpar(ftext(DHist2, fp_Tx)))
  
  DHistfile <- paste0(R_FOLDER,"/",NameF,"/",NameF,"Drought_History.png")
  DHistimg <- external_img(src = DHistfile, height = 2,width = 4) 

  FIG_11 <- block_list(
    fpar(ftext(paste0("Figure 11. 100-year (1920-2019) SPI time series (reversed axis) at ", SNameF ,
                       ". Dashed lines show, moderate (yellow), severe (red), and extreme (dark red), drought thresholds."), fp_Fig)))
  
  
  #################### Slide 18 ###############################
  
  
  DT_TIT<- block_list(
    fpar(ftext("Drought Events (1920-2019)", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext(SNameF, FTXTT3),                      fp_p = fp_par(text.align = "center")))


  TAB3 <- block_list(
    fpar(ftext(paste0("Table A1. SPI-12 drought characteristics at ",NameF, " identified in the SPI-12 timeseries. Duration is the number of months the drought persisted; ",
               "Average Intensity is the average absolute SPI; Peak Intensity is the highest SPI value calculated during the drought ", 
               "Magnitude is sum of absolute SPI values during the drought."), fp_Fig)))
   
 #NOTE Try and do the Auto FIt for this table to see if you can get it on the slide 
  #See Slide 10
  
  
  #################### Slide 19 ############################### 
  
  SPICNT<- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Drought Count.csv"),sep=",")
  
  SP_TIT<- block_list(
    fpar(ftext("Short-term vs Long-term Droughts", FTXTT),fp_p = fp_par(text.align = "center")))
  
  SPI_3_SMfile <- paste0(R_FOLDER,"/",NameF,"/",NameF,"Drought_HistoryS_3.png")
  SPI3_SMimg <- external_img(src = SPI_3_SMfile, height = 2,width = 4)
  
  SPI_12_SMfile <- paste0(R_FOLDER,"/",NameF,"/",NameF,"Drought_HistoryS_12.png")
  SPI12_SMimg <- external_img(src = SPI_12_SMfile, height = 2,width = 4)
  
  FIG_12 <- block_list(
    fpar(ftext(paste0("Figure 12. 30-year (1990-2019) SPI-3 time series (reversed axis) at ", NameF ,"."), fp_Fig)))
  FIG_13 <- block_list(
    fpar(ftext(paste0("Figure 13. 30-year (1990-2019) SPI-12 time series (reversed axis) at ", NameF ,"."), fp_Fig)))
  SPI3v12<- block_list(
    fpar(ftext(paste0("The SPI-3 provides a comparison of the precipitation over a specific 3-month period and reflects short- and medium-term moisture ", 
                       "conditions over the 30-year period (1990-2019). A total of ", SPICNT[1,3], " droughts were observed at ",
                        SNameF , ". Over this same time period, only ", SPICNT[1,4] ," droughts were identified when looking at the SPI-12 timeseries. It is important to compare the 3-month SPI with longer time scales. ",
                       "A relatively normal 3-month period could occur in the middle of a longer-term drought that would only be visible at longer time scales. Looking at longer time scales ",
                       "would prevent a misinterpretation that a drought might be over."), fp_Txa)))

  #################### Slide 20 ###############################
  
  S18_TIT<- block_list(
    fpar(ftext(paste("Fire Occurrence",ISL), FTXTT),fp_p = fp_par(text.align = "center")))

  
  FIRE <- paste0("Ecological drought in Hawaii is often drives an increase in wildfire occurrence. In Hawaii, wildfires are most extensive in dry ",
                 "and mesic non-native grasslands and shrublands. During drought events, wildfire risk in grasslands increases rapidly. ",
                 "Changes in land use that shift agricultural land to non-native cover of fire-prone grasses and shrubs ",
                 "combined with recurring incidences of drought are expected to increase the risk of future wildfire in Hawaii.")
                  
  
  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 20) 
  fp_Fire <- fpar(ftext(FIRE, fp_Tx))
  
  FireMfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Fire.png")
  FireMimg <- external_img(src = FireMfile, height = 4,width = 4) 
  
  
  FIG_14 <- block_list(
    fpar(ftext(paste0 ("Figure 14. The map shows wildfires that have occurred on the island of ", ISL, " between 1999 and 2019.") , fp_Fig)))
  
  
  ################ Slide 21 PART 4 Downscaling ###############
  
  P4_TIT<- block_list(
    fpar(ftext("Part 4: Future Climate", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext(SNameF, FTXTT3),fp_p = fp_par(text.align = "center")))
  
  Part4<- paste0("To simulate future rainfall and temperature, Global Climate Models are used. These can simulate future conditions under ",
                 "different scenarios for how much carbon dioxide we emit into the air. Two common scenarios are used: RCP 4.5 which assumes ",
                 "we reduce our carbon emissions, and RCP 8.5, a high emissions scenario. The outputs from global models are too coarse to ",
                 "accurately capture changes over the complex terrain of Hawaii. Therefore, we use an additional step called Downscaling to ",
                 "relate the global-scale information down to the local island scale. ")

  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 18) 
  fp_Part4 <- fpar(ftext(Part4, fp_Tx))

  
  
  M_Down <-  block_list(
    fpar(ftext("In Hawaii, two types of downscaled " , fp_NM3)),
    fpar(ftext("projections are available.", fp_NM3)),
    fpar(ftext("Dynamical Downscaling (End of Century)" , fp_NM3)),
    fpar(ftext("Statistical Downscaling (Mid & End of Century)" , fp_NM3)),
    fpar(ftext("Results for both types of downscaling ", fp_NM3)),
    fpar(ftext("and both scenarios will be shown here.", fp_NM3)))
   
  
  #################### Slide 22 ###############################
  
  S20_TIT<- block_list(
    fpar(ftext("Average Rainfall Change 2100", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext("Year 2100", FTXTT),fp_p = fp_par(text.align = "center")))
  
  
  Down <- read.csv(paste0(R_FOLDER,"/",NameF,"/",NameF," Downscaling.csv"),sep=",")

  RFA_Thresh100 <-  Down[c(1,4,9,12),2]
  RFD_Thresh100 <-  Down[c(2,5,10,13),2]
  RFW_Thresh100 <-  Down[c(3,6,11,14),2]
  RFA_Thresh40 <-   Down[c(15,18),2]
  RFD_Thresh40 <-   Down[c(16,19),2]
  RFW_Thresh40 <-   Down[c(17,20),2]
  
  #Annual RF
  RFA <- CLIMA[1,14]
  MinRF <- round(RFA * (min(RFA_Thresh100*0.01)),0)
  MaxRF <- round(RFA * (max(RFA_Thresh100*0.01)),0)

  #Dry Season 
  #DryM
  MinRFD <- round(DryM * (min(RFD_Thresh100*0.01)),0)
  MaxRFD <- round(DryM * (max(RFD_Thresh100*0.01)),0)

  #WetM
  MinRFW <- round(WetM * (min(RFW_Thresh100*0.01)),0)
  MaxRFW <- round(WetM * (max(RFW_Thresh100*0.01)),0)

  
 AnnualR  <- paste(min(RFA_Thresh100), " to ", max(RFA_Thresh100), "% Change")
 DryR  <- paste(min(RFD_Thresh100), " to ", max(RFD_Thresh100), "% Change")
 WetR  <- paste(min(RFW_Thresh100), " to ", max(RFW_Thresh100), "% Change") 
  
  
  RF_Down <-  block_list(
    fpar(ftext("Rainfall in the Year 2100" , fp_NM5)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext("Annual", fp_NM6)),
    fpar(ftext(paste0(MinRF," to ",MaxRF,RFUnit,"/year"), fp_NM3)),
    fpar(ftext(paste0("(",AnnualR,")"),fp_NM3)),
    fpar(ftext("Dry Season", fp_NM6)),
    fpar(ftext(paste0(MinRFD," to ",MaxRFD,RFUnit,"/month"), fp_NM3)),
    fpar(ftext(paste0("(",DryR,")"),fp_NM3)),
    fpar(ftext("Wet Season", fp_NM6)),
    fpar(ftext(paste0(MinRFW," to ",MaxRFW,RFUnit,"/month"), fp_NM3)),
    fpar(ftext(paste0("(",WetR,")"),fp_NM3)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext(paste0("The range in projections include estimates for both low emissions ",
                       "(RCP 4.5) and high emissions     (RCP 8.5) scenarios, and for both ",
                       "Dynamical and Statistical Downscaling approaches."), fp_NM2)))
    
  RF10085file <- paste0(R_FOLDER,"/",NameF,"/",NameF," DS_RF_8.5.png")
  RF10085img <- external_img(src = RF10085file, height = 4,width = 4) 
  
  FIG_15<- block_list(
    fpar(ftext(paste0("Figure 15. Downscaled future rainfall projections (% Change; ", 
                      "(2100) at ",SNameF,", Dynamical Downscaling (DyDS), ",
                      "Statistical Downscaling (StDs), for annual (ANN), dry season (DRY) ", 
                      "and wet season (WET)."), fp_Fig)))

  
##################### Slide 23 #################################################
  
  S21_TIT<- block_list(
    fpar(ftext("Average Rainfall Change", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext("2040-2070", FTXTT),fp_p = fp_par(text.align = "center")))
  
  
  #Annual RF
  RFA <- CLIMA[1,14]
  MinRF4 <- round(RFA * (min(RFA_Thresh40*0.01)),0)
  MaxRF4 <- round(RFA * (max(RFA_Thresh40*0.01)),0)

  #Dry Season 
  #DryM
  MinRFD4 <- round(DryM * (min(RFD_Thresh40*0.01)),0)
  MaxRFD4 <- round(DryM * (max(RFD_Thresh40*0.01)),0)

  #WetSeason
  #WetM
  MinRFW4 <- round(WetM * (min(RFW_Thresh40*0.01)),0)
  MaxRFW4 <- round(WetM * (max(RFW_Thresh40*0.01)),0)

  AnnualRM  <- paste0(min(RFA_Thresh40), " to ", max(RFA_Thresh40), "% Change")
  DryRM  <- paste0(min(RFD_Thresh40), " to ", max(RFD_Thresh40), "% Change")
  WetRM  <- paste0(min(RFW_Thresh40), " to ", max(RFW_Thresh40), "% Change") 
  
  
  RF_Down4 <-  block_list(
    fpar(ftext("Rainfall for Years 2040-2070" , fp_NM5)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext("Annual", fp_NM6)),
    fpar(ftext(paste0(MinRF4," to ",MaxRF4,RFUnit,"/year"), fp_NM3)),
    fpar(ftext(paste0("(",AnnualRM,")"),fp_NM3)),
    fpar(ftext("Dry Season", fp_NM6)),
    fpar(ftext(paste0(MinRFD4," to ",MaxRFD4,RFUnit,"/month"), fp_NM3)),
    fpar(ftext(paste0("(",DryRM,")"),fp_NM3)),
    fpar(ftext("Wet Season", fp_NM6)),
    fpar(ftext(paste0(MinRFW4," to ",MaxRFW4,RFUnit,"/month"), fp_NM3)),
    fpar(ftext(paste0("(",WetRM,")"),fp_NM3)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext(paste0("The range in projections include estimates for both low emissions ", 
                      "(RCP 4.5) and high emissions     (RCP 8.5) scenarios, for the ",
                      "Statistical Downscaling approach."), fp_NM2)))
  
  #RF40file <- paste0(R_FOLDER,"/",NameF,"/",NameF," StDsRF2040.png")
  #RF40img <- external_img(src = RF40file, height = 6,width = 4) 
  
  RF40file <- paste0(R_FOLDER,"/",NameF,"/",NameF," StDsRF2040.png")
  RF40img <- external_img(src = RF40file) 
  
  FIG_16<- block_list(
           fpar(ftext(paste0("Figure 16. Downscaled future rainfall projections (% Change; ",
                             "2040-2070) at ",SNameF,", for the Statistical Downscaling (StDs) ",
                             "approach, for annual (ANN), dry season (DRY), and wet season ",
                             "(WET) for RCP 4.5 (left) and RCP 8.5 (right)."), fp_Fig)))

  
  
#################### Slide 24 Temperature 2100
  
  S22_TIT<- block_list(
    fpar(ftext("End-of-Century Change in Temperature", FTXTT),fp_p = fp_par(text.align = "center")))
   
  
  TAA_Thresh100 <-  Down[c(7,8,21,22),2]
  
  #Annual RF
  TAA <- CLIMA[3,14]
  MinTA <- round(TAA + min(TAA_Thresh100),1)
  MaxTA <- round(TAA + max(TAA_Thresh100),1)
  AvgTa <- round(mean(TAA_Thresh100),1)
  
  
  Tchg2100  <- paste0("(",min(TAA_Thresh100), TUnit, " to ", max(TAA_Thresh100), TUnit," Change)") 
  

  TA_Down <-  block_list(
    fpar(ftext( "                    ", fp_Fig2)),
    fpar(ftext("Mean Temperature Now" , fp_NM6)),
    fpar(ftext(paste0(TAA,TUnit), fp_NM3)),
    fpar(ftext( "                    ", fp_Fig2)),
    fpar(ftext("Mean Temperature 2100" , fp_NM6)),
    fpar(ftext(paste0(MinTA,TUnit," to ", MaxTA,TUnit), fp_NM3)),
    fpar(ftext(Tchg2100, fp_NM3)),
    fpar(ftext("                    ", fp_Fig2)),
    fpar(ftext("                    ", fp_Fig2)),
    fpar(ftext(paste0("The range in projections include estimates for both low emission ", 
              "(RCP4.5) and high emissions (RCP8.5) scenarios, for both the ", 
              "Dynamical and Statistical Downscaling approaches."), fp_NM2)))
  
  
  
  TA100file <- paste0(R_FOLDER,"/",NameF,"/",NameF," DS_Temp2100.png")
  TA100img <- external_img(src = TA100file) 
  
  FIG_17 <- block_list(
    fpar(ftext(paste0("Figure 17. Downscaled projected change in mean temperature ",  
                      "(Year 2100) at ",SNameF,", Dynamical Downscaling (DyDs), ",  
                      "Statistical Downscaling (StDs)."), fp_Fig)))
  
  
#################### Slide 25 Temperature 2040
  
  S23_TIT<- block_list(
    fpar(ftext("Mid-Century Change in Temperature", FTXTT),fp_p = fp_par(text.align = "center")))

  TAA_Thresh40 <-  Down[c(23,24),2]
  MinTA40 <- round(TAA + min(TAA_Thresh40),1)
  MaxTA40 <- round(TAA + max(TAA_Thresh40),1)
  #Annual RF
  TAA <- CLIMA[3,14]
  
  Tchg2040  <- paste0("(",min(TAA_Thresh40),TUnit, " to ", max(TAA_Thresh40), TUnit," Change)") 
 
  AvgTa4 <- round(mean(TAA_Thresh40),1)
  
  TA_Down4 <-  block_list(
    fpar(ftext( "                    ", fp_Fig2)),
    fpar(ftext("Mean Temperature Now" , fp_NM6)),
    fpar(ftext(paste0(TAA,TUnit), fp_NM3)),
    fpar(ftext( "                    ", fp_Fig2)),
    fpar(ftext("Mean Temperature 2040-2070" , fp_NM6)),
    fpar(ftext(paste0(MinTA40,TUnit," to ", MaxTA40,TUnit), fp_NM3)),
    fpar(ftext(Tchg2040, fp_NM3)),
    fpar(ftext( "                    ", fp_Fig2)),
    fpar(ftext(        "                    ", fp_Fig2)),
    fpar(ftext(paste0("The range in projections include estimates for both low emission ", 
              "(RCP 4.5) and high emissions (RCP 8.5) scenarios, for the ", 
               "Statistical Downscaling approach."), fp_NM2)))

  
  TA40file <- paste0(R_FOLDER,"/",NameF,"/",NameF," StDs_Temp2040.png")
  TA40img <- external_img(src = TA40file) 
  
  FIG_18 <- block_list(
    fpar(ftext(paste0("Figure 18. Downscaled Future temperature projections ",
                      "(2040-2070) at ",SNameF," for the Statistical ", 
                      "Downscaling (StDs) approach."), fp_Fig)))
  
  
  ################ Slide 26 Summary and Conclusions ###############
  
  
  S24_TIT<- block_list(
    fpar(ftext("Part 5: CCVD Summary", FTXTT),fp_p = fp_par(text.align = "center")),
    fpar(ftext(SNameF, FTXTT3),fp_p = fp_par(text.align = "center")))
  
  
  DecD <- (CNTDRT/10)
  
  Summary<- paste0(NameF," (",SNameF,") is located on the island of ", ISL, " at mean elevation of ", El_Mean, ELUnit," (range: ",El_Min, " to ", El_Max,ELUnit,"). ",  
                 "Rainfall varies over the course of the year with a maximum of ",RFx, RFUnit3,", occurring in ", MONLIST[RFxm], " and a minimum of ",RFn, RFUnit3,
                 " occurring in ", MONLIST[RFnm],". On average, wet season months (Nov-Apr) receive ",SeaD,RFUnit," of more rainfall than dry season months (May-Oct). ",
                 "Seasonal rainfall can vary within the unit as well, with dry season rainfall ranging from ",D_RF_R, RFUnit3," and wet season rainfall ",
                 "ranging from ",W_RF_R, RFUnit3," across the ",EL_Dif,ELUnit, " elevation gradient. ",
                 "Rainfall can also vary considerably from year-to-year with the driest years occurring during a Strong El Ni?o event, when on ",
                 "average, ", DifWRF,"% less rainfall is received, relative to the long-term average. ",
                 "The average temperature at ",SNameF," is ", TAA,TUnit," but temperature ranges from ",TAn, TUnit, " to ", TAx, TUnit," over the course of the year. ",
                 "Drought is a reoccurring feature in the climate system of ", SNameF, " with a total of ",
                 CNTDRT, " occurring over the record which is approximately ",DecD," per decade. A total of ",SoG, " drought events were at severe strength ",
                 "or greater and the longest drought lasted for a total of ", LongD, " consecutive months. ",
                 "Future projections of rainfall are uncertain, with end-of-century annual changes ranging from ", min(RFA_Thresh100)," to " ,max(RFA_Thresh100), "% and more pronounced changes occurring during the wet and dry seasons. ",
                 "Future projections of temperature suggest an increase of ", min(TAA_Thresh40), TUnit, 
                 " to ", max(TAA_Thresh40),TUnit," by mid century (2040-2070) ",
                 "and an increase of ", min(TAA_Thresh100),TUnit," to ",  max(TAA_Thresh100), TUnit," by the end of the century (2100). ")
                 
                 fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 17) 
                 fp_Summary <- fpar(ftext(Summary, fp_Tx))

 
   ################ Slide 27 Resource PAGE ############### 
                 
                 #Note See if you can make clickable Links 
                 
                 
                 S25_TIT<- block_list(
                   fpar(ftext("External Resources", FTXTT),fp_p = fp_par(text.align = "center")))
                 
                 fp_RES <-  block_list(
                   fpar(ftext( "For more Information"   , fp_NM5)),
                   fpar(ftext( "                    "   , fp_Fig2)),
                   fpar(ftext("US Drought Monitor"      , fp_NM6)),
                   fpar(ftext("https://droughtmonitor.unl.edu/", fp_NM3)),
                   fpar(ftext( "                    ", fp_Fig2)),
                   fpar(ftext("State of Hawaii Drought Plan"      , fp_NM6)),
                   fpar(ftext(" https://files.hawaii.gov/dlnr/cwrm/planning/HDP2017.pdf", fp_NM3)),
                   fpar(ftext( "                    ", fp_Fig2)),
                   fpar(ftext("Rainfall Atlas of Hawaii" , fp_NM6)),
                   fpar(ftext("http://rainfall.geography.hawaii.edu/", fp_NM3)),
                   fpar(ftext( "                    ", fp_Fig2)),
                   fpar(ftext( "Climate Of Hawaii", fp_NM6)),
                   fpar(ftext("http://climate.geography.hawaii.edu/", fp_NM3)),
                   fpar(ftext( "                    ",   fp_Fig2)),
                   fpar(ftext("ENSO Current Phase and Discussion" , fp_NM6)),
                   fpar(ftext("https://www.cpc.ncep.noaa.gov/products/analysis_monitoring/enso_advisory/ensodisc.shtml", fp_NM3)),
                   fpar(ftext( "     "   , fp_Fig2)),
                   fpar(ftext( "Pacific Fire Exchange", fp_NM6)),
                   fpar(ftext("https://www.pacificfireexchange.org/", fp_NM3)),
                   fpar(ftext( "                    ", fp_Fig2)))
                 
                 #ph_hyperlink(doc,
                 #             ph_label = "Content Placeholder 2", href = "https://cran.r-project.org")
                 
  ################ Slide 28 Acknowledgements ###############
                
                 
                 S26_TIT<- block_list(
                   fpar(ftext("Acknowledgements", FTXTT),fp_p = fp_par(text.align = "center")))
       
  
  Ack <- paste0( "Mari-Vaughn Johnson, Heather Kerkering, Darcy Yogi (PI-CASC), Victoria Keener,",
                 " Laura Brewington (East-West Center Pacific-RISA), ",
                    "Melissa Kunz, Clay Trauernicht (NREM, UH Manoa), Katie Kamelamela (USDA, Forest Service), ",
                    "Thomas Giambelluca (WRRC, UH Manoa), and John Marra (NOAA).")
  
  fp_Ack2 <- fp_text(color = "black", font.size = 18) 
  fp_Ack <- fpar(ftext(Ack, fp_Ack2 ))
  
  
  fp_CITA1 <- fp_text(bold= TRUE, color = "darkred", font.size = 18,underlined = TRUE) 
  fp_CITA2 <- fp_text(color = "black", font.size = 18) 
  fp_CITA3 <- fp_text(color = "red", font.size = 16) 
  
  CITA <- block_list(
    fpar(ftext       ("Suggested Citation", fp_CITA1)),
    fpar(ftext(paste0("Longman, R.J., Ford, D.J., Frazier, A.G, Giardina, C.P (2021). ",
                      "Climate Change, Climate Variability and Drought Portfolio ", NameF,
                      " Pacific Drought knowledge Exchange CCVD Series Version 3.0, XXpp."), fp_CITA2)))


  Draft <- block_list(
    fpar(ftext       ("For the most up-to-date version of this portfolio contact Ryan Longman: rlongman@Hawai'i.edu for more information.", fp_CITA3)))

   ############## Slide 29 Works Cited ###########
  
  S27_TIT<- block_list(
    fpar(ftext("Works Cited", FTXTT),fp_p = fp_par(text.align = "center")))
  
  Worksfile <- paste0(I_FOLDER,"Works.png")
  Worksimg <- external_img(src =Worksfile , height = 1, width = 1) 
  
  
  
  ############### Slide 30 ANNEX I ###############
  
  
  RF100<- paste0("The 100-year monthly rainfall dataset was drawn from two unique gridded products. We used data from Frazier et al. (2016) ",
                 "for the period 1920-1989 and Lucas et al. (In Review) for the period 1990-2019. Given that two unique data sets and methods ",
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
  

  
    ################ Slide 31 ANNEX II #############
  

  DownA<- paste0("Two types of downscaling products were used in this analysis. Here we explain some of the nuances between the two. ",
                 "Dynamical Downscaling (Zhang et al., 2016), feeds GCM output into a regional model that can account for local topographic and ",
                 "atmospheric phenomena at much finer resolutions (e.g. 1 km). End-of-century projections (2100) encompass the period 2080-2099. ",
                 "Statistical Downscaling ",
                 "(Elison Timm et al., 2015, Elison Timm, 2017), develops a relationship between GCM model output and station data for a historical ",
                 "period and then uses this established relationship to make projections for two future scenarios. End-of-century projections (2100) ",
                 "encompass the period 2070-2099 (2100), Mid-century projections encompass the period 2040-2070.")

       
  fp_Tx <- fp_text(italic = TRUE, color = "black", font.size = 18) 
  fp_DownA <- fpar(ftext(DownA, fp_Tx))
  
 
  
  ################ Slide 32 ANNEX III #############
  
  Roadfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Forest_Roads.png")
  
  Roadimg <- external_img(src = Roadfile) 
  
  Trailfile <- paste0(R_FOLDER,"/",NameF,"/",NameF," Trails.png")
  Trailimg <- external_img(src = Trailfile) 
  
  FIG_A2 <- block_list(
    fpar(ftext(paste0("Figure A2: Forest roads ",ISL,"."), fp_Fig)))
  FIG_A3 <- block_list(
    fpar(ftext(paste0("Figure A3: Trail Inventory ",ISL,"."), fp_Fig)))
  
  
  
  
  
  
  
############################################################################################################################################################  
###############
###############  SLIDE DECK 
###############

  #Look up office themes for R PPT 
  # Do a review of existing software to see if there are any updates now packages for this 
  
#Slide 1
  mypowerpoint <- read_pptx() %>%
    add_slide("Title Slide","Office Theme") %>%
    #Add Title
       ph_with(TIT,ph_location_type("ctrTitle")) %>%
    #Add dates left footer NOTE#Check ph_location_type
       ph_with(value = format(Sys.Date()), location = ph_location_type(type = "dt"))%>%
    #Add subtitle 
       ph_with(SUB,ph_location_type("subTitle")) %>%
    #add logos 
       ph_with(value = EWCimg, location = ph_location(label = "my_name",
                      left = 4.38, top = 6.0, width = 1.4, height = 1.4)) %>%
       ph_with(value = CASCimg, location = ph_location(label = "my_name",
                      left = 2.1, top = 6.2, width = 1, height = 1)) %>%
       ph_with(value = FSimg, location = ph_location(label = "my_name",
                      left = 6.9, top = 6.2, width = 1, height = 1)) %>%
       ph_with(value = PDKE_L, location = ph_location(label = "my_name",
                      left = 0, top = 0.3, width = 10, height = 2))%>%

  
    #Slide 2 
    add_slide("Title and Content","Office Theme") %>%
      #ph_with(S2_TIT, ph_location_type("title",position_left = TRUE)) %>%
      ph_with(fp_HDKE,ph_location_type("body"))%>%
      ph_with(value = "2", location = ph_location_type(type = "sldNum"))%>%
      ph_with(value = EWCimg, location = ph_location(label = "my_name",
                      left = 1.1, top = 6.1, width = 1.4, height = 1.4)) %>%
      ph_with(value = CASCimg, location = ph_location(label = "my_name",
                      left = 2.6, top = 6.3, width = 1, height = 1)) %>%
      ph_with(value = FSimg, location = ph_location(label = "my_name",
                      left = 3.8, top = 6.3, width = 1, height = 1)) %>%
      # ph_with(value = RISAimg, location = ph_location(label = "my_name",
      #                 left = 4.2, top = 6.3, width = 1.4, height = 1)) %>%
      ph_with(value = WRRCimg, location = ph_location(label = "my_name",
                      left = 5, top = 6.3, width = 1.1, height = 1)) %>%
      ph_with(value = NOAAimg, location = ph_location(label = "my_name",
                      left = 6.3, top = 6.3, width = 1, height = 1)) %>%
      ph_with(value = UHimg, location = ph_location(label = "my_name",
                      left = 7.5, top = 6.3, width = 1, height = 1)) %>%
     ph_with(value = PDKE_L, location = ph_location(label = "my_name",
                      left = 0, top = 0, width = 10, height = 2))%>%
    
      
#Slide 3 
    add_slide("Two Content","Office Theme") %>%
      ph_with(S3_TIT,         ph_location_type("title")) %>%
      ph_with(fp_CCVD,        ph_location_type("body",position_right = FALSE)) %>%
      ph_with(value = "3", location = ph_location_type(type = "sldNum")) %>%
      ph_with(value = MAPimg, ph_location_type("body",position_right = TRUE)) %>%
      ph_with(value = "Digital elevation model NAD84", location = ph_location_type(type = "dt"))%>%
      ph_with(value = FIG_1, location = ph_location(label = "my_name",
                      left = 5.6, top = 6.4, width = 3.81, height = 0.77))%>%
  
    
#Slide 4  PART 1  
  add_slide("Title and Content","Office Theme") %>%
  ph_with(P1_TIT,ph_location_type("title",position_left = TRUE)) %>%
  ph_with(fp_PART1,ph_location_type("body"))%>%
  ph_with(value = "4", location = ph_location_type(type = "sldNum"))%>%
  ph_with(value = MODimg, location = ph_location(label = "my_name",
                                                 left = 1.1, top = 5, width = 2.6, height = 2)) %>%
  ph_with(value = RGimg, location = ph_location(label = "my_name",
                                                  left = 3.8, top = 5, width = 2.4, height = 2)) %>%
  ph_with(value = RFimg, location = ph_location(label = "my_name",
                                                  left = 6.4, top = 5, width = 2.4, height = 2)) %>%

    
#Slide 5 
    add_slide("Two Content","Office Theme") %>%
    ph_with(S5_TIT,       ph_location_type("title")) %>%
    ph_with(M_Elev,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(value = "5", location = ph_location_type(type = "sldNum")) %>%
    ph_with(value = Elevimg, ph_location_type("body",position_right = TRUE)) %>%
    ph_with(value = FIG_2, location = ph_location(label = "my_name",
                                                  left = 5.1, top = 5.6, width = 4.2, height = 1.4))%>%
#Slide 6 
  add_slide("Two Content","Office Theme") %>%
    ph_with(S6_TIT,           ph_location_type("title")) %>%
    ph_with(M_Climate,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(value = "6", location = ph_location_type(type = "sldNum")) %>%
    ph_with(value = Climimg, ph_location_type("body",position_right = TRUE)) %>%
    ph_with(value = "Giambelluca et al. (2013;2014)", location = ph_location_type(type = "dt"))%>%
    ph_with(value = FIG_3, location = ph_location(label = "my_name",
                     left = 5.1, top = 6, width = 4.4, height = 1.4))%>%
  
#Slide 7 
  add_slide("Two Content","Office Theme") %>%
    ph_with(S7_TIT,          ph_location_type("title")) %>%
    ph_with(fp_CLIMO,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(value = "7", location = ph_location_type(type = "sldNum")) %>%
    ph_with(value = Climoimg, ph_location_type("body",position_right = TRUE)) %>%
    ph_with(value = "Giambelluca et al. (2013;2014)", location = ph_location_type(type = "dt"))%>%
    ph_with(value = FIG_4, location = ph_location(label = "my_name",
                                                  left = 5.2, top = 6.2, width = 4, height = 1.4))%>%

#Slide 8 
add_slide("Two Content","Office Theme") %>%
  ph_with(S8_TIT, ph_location_type("title")) %>%
  ph_with(value = RFFimg, ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = "8", location = ph_location_type(type = "sldNum")) %>%
  ph_with(value = TAimg, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_5, location = ph_location(label = "my_name",
                                                left = 0.5, top = 6.08, width = 4.5, height = 1.4))%>%
  ph_with(value = FIG_6, location = ph_location(label = "my_name",
                                                left = 5.1, top = 6.08, width = 4.5, height = 1.4))%>%
  ph_with(my_pres, value = "Giambelluca et al. (2013;2014)", location = ph_location_type(type = "dt"))%>%

#Slide 9 
add_slide("Two Content","Office Theme") %>%
  ph_with(S9_TIT,      ph_location_type("title")) %>%
  ph_with(fp_SEA,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = "9", location = ph_location_type(type = "sldNum")) %>%
  ph_with(value = SEAimg, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_7, location = ph_location(label = "my_name",
                                                left = 6, top = 6.2, width = 3, height = 1.4))%>%
  ph_with(value = "Giambelluca et al. (2013)", location = ph_location_type(type = "dt"))%>%
    
    
#####################################################################################    
  
  #Slide 10
  add_slide("Title and Content","Office Theme") %>%
    ph_with(S10_TIT,      ph_location_type("title")) %>%
    ph_with(value = "10", location = ph_location_type(type = "sldNum"))%>% 
    #ph_with(value = ACLIM_T, ph_location_type("body")) %>%
    ph_with(value = "Giambelluca et al. (2013;2014)", location = ph_location_type(type = "dt"))%>%
    ph_with(value =   TAB1, location = ph_location(label = "my_name",
                                                   left = 0.2, top = 4.45, width = 9.5, height = 4))%>%
    ph_with(value =   ACLIM_T, location = ph_location(label = "my_name",
                                                      left = 0.3, top = 1.8, width = 7, height = 3.5))%>%

    
#Slide 11   
add_slide("Two Content","Office Theme") %>%
  ph_with(P2_TIT,ph_location_type("title",position_left = TRUE)) %>%
  ph_with(fp_InterAn,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = "11", location = ph_location_type(type = "sldNum"))%>%
  ph_with(value = "https://www.climate.gov/enso", location = ph_location_type(type = "dt"))%>%
  ph_with(value = ENSO2img, location = ph_location(label = "my_name",
                                           left = 5.5, top = 2, width = 3, height = 4)) %>%
  
    
    
  
  #Slide 12 
add_slide("Two Content","Office Theme") %>%
  ph_with(EN_TIT,ph_location_type("title")) %>%
  ph_with(fp_MEIW,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = "12", location = ph_location_type(type = "sldNum")) %>%
  ph_with(value = MEIWimg, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_8, location = ph_location(label = "my_name",
                                                left = 5.2, top = 6.05, width = 4.1, height = 1))%>%
  ph_with(my_pres, value = "Frazier et al. (2016); Wolter and Timlin (2011)", location = ph_location_type(type = "dt"))%>%
  ph_with(value = ENSOKimg, location = ph_location(label = "my_name",
                                                  left = 1.5, top = 5.5, width = 2.1, height = 1.1)) %>%
  ph_with(value = TAB2, location = ph_location(label = "my_name",
                                                  left = 0.8, top = 4.5, width = 4, height = 1.4))%>%

  #Slide 13 
  add_slide("Two Content","Office Theme") %>%
  ph_with(S13_TIT,         ph_location_type("title")) %>%
  ph_with(fp_Trend,        ph_location_type("body",position_right = FALSE)) %>%
  #ph_with(value = "13", location = ph_location_type(type = "sldNum")) %>%
  ph_with(value = Trendimg, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_9, location = ph_location(label = "my_name",
                                                left = 5.3, top = 6.4, width = 4.3, height = 1.4))%>%
  ph_with(my_pres, value = "Frazier et al. (2016); Lucas et al. (In Prep); See Annex I", location = ph_location_type(type = "dt"))%>%
 
  
#Slide 14 #Part 3
  add_slide("Title and Content","Office Theme") %>%
  ph_with(P3_TIT,ph_location_type("title",position_left = TRUE)) %>%
  ph_with(fp_Part3,ph_location_type("body"))%>%
  ph_with(value = "14", location = ph_location_type(type = "sldNum"))%>%
  ph_with(value = Droughtimg, location = ph_location(label = "my_name",
                                                 left = 2.5, top = 4.5, width = 2, height = 2)) %>%
  ph_with(value = Fireimg, location = ph_location(label = "my_name",
                                                 left = 5.5, top = 4.5, width = 2, height = 2)) %>%
    
###  Slide 15    
    
    add_slide("Title and Content","Office Theme") %>%
    ph_with(TIT_D,ph_location_type("title",position_left = TRUE)) %>%
    ph_with(value = "Frazier et al. (2019)", location = ph_location_type(type = "dt"))%>%
    ph_with(value = "15", location = ph_location_type(type = "sldNum"))%>%
    ph_with(value = D_METimg, location = ph_location(label = "my_name",
                                                     left = 0.3, top = 1.8, width = 2, height = 1.3)) %>%
    ph_with(value = D_AGimg, location = ph_location(label = "my_name",
                                                    left = 0.3, top = 3.6, width = 2, height = 1.3)) %>%
    ph_with(value = D_HYDimg, location = ph_location(label = "my_name",
                                                     left = 0.3, top = 5.4, width = 2, height = 1.3)) %>%
    ph_with(value =  D_ECOimg, location = ph_location(label = "my_name",
                                                      left = 5, top = 2.6, width = 2, height = 1.3))%>%
    ph_with(value =  D_SOCimg, location = ph_location(label = "my_name",
                                                      left = 5, top = 4.4, width = 2, height = 1.3))%>%
    ## TEXT
    ph_with(value =  fp_D_M, location = ph_location(label = "my_name",
                                                    left = 2.4, top = 1.8, width = 2.6, height = 0.69))%>%
    ph_with(value =  fp_D_A, location = ph_location(label = "my_name",
                                                    left = 2.4, top = 3.6, width = 2.6, height = 1.17))%>% 
    ph_with(value =  fp_D_H, location = ph_location(label = "my_name",
                                                    left = 2.4, top = 5.4, width = 2.6, height = 1.17))%>% 
    ph_with(value =  fp_D_E, location = ph_location(label = "my_name",
                                                    left = 7.1, top = 2.6, width = 2.6, height = 1.3))%>% 
    ph_with(value =  fp_D_S, location = ph_location(label = "my_name",
                                                    left = 7.1, top = 4.3, width = 2.78, height = 1.56))%>% 
  
  
  
  #Slide 16 
  add_slide("Two Content","Office Theme") %>%
  ph_with(S15_TIT,       ph_location_type("title")) %>%
  ph_with(fp_SPI,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = "16", location = ph_location_type(type = "sldNum")) %>%
  #ph_with(value = SPIimg, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_10, location = ph_location(label = "my_name",
                                                left = 5.3, top = 5, width = 4.5, height = 1.4))%>%
  ph_with(value = SPIimg, location = ph_location(label = "my_name",
                                                left = 5, top = 2.5, width = 4.42, height = 2.72))%>%
  ph_with(my_pres, value = "Frazier et al. (2016); Lucas et al. (In Prep)", location = ph_location_type(type = "dt"))%>%

 
  #Slide 17

  add_slide("Two Content","Office Theme") %>%
  ph_with(S16_TIT,         ph_location_type("title")) %>%
  ph_with(fp_DHist,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = "17", location = ph_location_type(type = "sldNum")) %>%
  #ph_with(value = DHistimg, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_11, location = ph_location(label = "my_name",
                             left = 5.3, top = 5, width = 4.5, height = 1.4))%>%
  ph_with(value = DHistimg, location = ph_location(label = "my_name",
                             left = 5, top = 2.5, width = 4.42, height = 2.72))%>%
  ph_with(my_pres, value = "Frazier et al. (2016); Lucas et al. (In Prep)", location = ph_location_type(type = "dt"))%>% 
 
    
    
  #   #Slide 18 #Drought Table (Derek moved to annex)
  # 
  #   
  # mypowerpoint <- read_pptx() %>%
  # 
  #   add_slide("Title and Content","Office Theme") %>%
  #   ph_with(DT_TIT,ph_location_type("title",position_left = TRUE)) %>%
  #   ph_with(value = "18", location = ph_location_type(type = "sldNum"))%>% 
  #   ph_with(value = DrotT_T1, ph_location_type("body")) %>%
  #   ph_with(value = TAB3, location = ph_location(label = "my_name",
  #                                                left = 7.7, top = 2.5, width = 2.2, height = 3.5)) %>%
  # 
  #   {if(Drot_ct>12) add_slide(.,"Title and Content", "Office Theme") %>% 
  #              ph_with(DT_TIT,ph_location_type("title",position_left = TRUE)) %>%
  #              ph_with(value = "19", location = ph_location_type(type = "sldNum"))%>%
  #              ph_with(value = DrotT_T2, ph_location_type("body")) %>%
  #              ph_with(value = TAB3,
  #                      location = ph_location(label = "my_name",
  #                                             left = 7.7, top = 2.5, width = 2.2, height = 3.5)) else .} %>%
  # 
  # {if(Drot_ct>24) add_slide(.,"Title and Content", "Office Theme") %>% 
  #     ph_with(DT_TIT,ph_location_type("title",position_left = TRUE)) %>%
  #     ph_with(value = "20", location = ph_location_type(type = "sldNum"))%>%
  #     ph_with(value = DrotT_T3, ph_location_type("body")) %>%
  #     ph_with(value = TAB3,
  #             location = ph_location(label = "my_name",
  #                                    left = 7.7, top = 2.5, width = 2.2, height = 3.5)) else .}
  # 
  # 
  #   print(mypowerpoint , target = paste0(P_FOLDER,NameF,"_CCVD_Portfolio_V2_Working",RFUnit2,"_TEST.pptx"))
  # 
  # 
  # 
    
    #Slide 18
    add_slide("Title and Content","Office Theme") %>%
    ph_with(SP_TIT,ph_location_type("title",position_left = TRUE)) %>%
    #ph_with(value = "Lucas et al. (In Review)", location = ph_location_type(type = "dt"))%>%
    ph_with(value = "18", location = ph_location_type(type = "sldNum"))%>%
    ph_with(value = SPI3_SMimg, location = ph_location(label = "my_name",
                                                       left = 0.4, top = 1.8, width = 4.25, height = 2.62)) %>%
    ph_with(value = SPI12_SMimg , location = ph_location(label = "my_name",
                                                         left = 5.1, top = 1.8, width = 4.25, height = 2.62)) %>%
    ph_with(value = FIG_12, location = ph_location(label = "my_name",
                                                   left = 0.5, top = 4.35, width = 4, height = 0.5)) %>%
    ph_with(value = FIG_13, location = ph_location(label = "my_name",
                                                   left = 5.2, top = 4.35, width = 4, height = 0.5)) %>%
    ph_with(value = SPI3v12, location = ph_location(label = "my_name",
                                                    left = 0.2, top = 4.85, width = 9.68, height = 2.5)) %>%
    
    
    
   #Slide 19
   add_slide("Two Content","Office Theme") %>%
   ph_with(S18_TIT,  ph_location_type("title")) %>%
   ph_with(fp_Fire,        ph_location_type("body",position_right = FALSE)) %>%
   ph_with(value = "19", location = ph_location_type(type = "sldNum")) %>%
   ph_with(value = FireMimg, ph_location_type("body",position_right = TRUE)) %>%
   ph_with(value = FIG_14, location = ph_location(label = "my_name",
                                                 left = 5.3, top = 6.1, width = 4.3, height = 1.4))%>%
   ph_with(my_pres, value = "Trauernicht, 2019; Frazier et al. (In Review)", location = ph_location_type(type = "dt"))%>% 
  

   #Slide 20 #Part 4
   add_slide("Title and Content","Office Theme") %>%
   ph_with(P4_TIT, ph_location_type("title",position_left = TRUE)) %>%
   ph_with(fp_Part4,ph_location_type("body"))%>%
   ph_with(value = "20", location = ph_location_type(type = "sldNum"))%>%
   ph_with(my_pres, value = "See Annex II", location = ph_location_type(type = "dt"))%>% 
   ph_with(value = Downimg, location = ph_location(label = "my_name",
                                                     left = 6, top = 4.5, width = 3, height = 2)) %>%
  ph_with(value = M_Down, location = ph_location(label = "my_name",
                                                     left = 0.5, top = 4.5, width = 6.3, height = 2)) %>%
  
  
  #Slide 21
  add_slide("Two Content","Office Theme") %>%
  ph_with(S20_TIT,   ph_location_type("title")) %>%
  ph_with(RF_Down,        ph_location_type("body",position_right = FALSE)) %>%
  ph_with(value = "21", location = ph_location_type(type = "sldNum")) %>%
  #ph_with(value = RF10085img, ph_location_type("body",position_right = TRUE)) %>%
  ph_with(value = FIG_15, location = ph_location(label = "my_name",
                                                 left = 5.2, top = 5.5, width = 4, height = 1.4))%>%
  ph_with(my_pres, value = "Elison Timm et al., 2015; Zhang et al., 2016; See Annex II", location = ph_location_type(type = "dt"))%>% 
  ph_with(value = RF10085img, location = ph_location(label = "my_name",
                                                   left = 4.6, top = 2, width = 5.17, height = 3.51)) %>%
  
    #Slide 22
    add_slide("Two Content","Office Theme") %>%
    ph_with(S21_TIT,         ph_location_type("title")) %>%
    ph_with(RF_Down4,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(value = "22", location = ph_location_type(type = "sldNum")) %>%
    #ph_with(value = RF40img, ph_location_type("body",position_right = TRUE))  %>%
    ph_with(value = FIG_16, location = ph_location(label = "my_name",
                                                   left = 5.2, top = 5.5, width = 4, height = 1.4))%>%
    ph_with(my_pres, value = "Elison Timm et al., 2015; See Annex II", location = ph_location_type(type = "dt"))%>% 
    ph_with(value = RF40img, location = ph_location(label = "my_name",
                                                       left = 4.6, top = 2, width = 5.17, height = 3.51)) %>%
    
    
    
    #Slide 23
    add_slide("Two Content","Office Theme") %>%
    ph_with(S22_TIT,       ph_location_type("title")) %>%
    ph_with(TA_Down,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(value = "23", location = ph_location_type(type = "sldNum")) %>%
    #ph_with(value = TA100img, ph_location_type("body",position_right = TRUE)) %>%
    ph_with(value = FIG_17, location = ph_location(label = "my_name",
                                                   left = 5.2, top = 5.5, width = 4, height = 1.4))%>%
    ph_with(my_pres, value = "Elison Timm 2017; Zhang et al., 2016; See Annex II", location = ph_location_type(type = "dt"))%>% 
    ph_with(value = TA100img, location = ph_location(label = "my_name",
                                                       left = 4.6, top = 2, width = 5.17, height = 3.51)) %>%
    
    #Slide 24
    add_slide("Two Content","Office Theme") %>%
    ph_with(S23_TIT,         ph_location_type("title")) %>%
    ph_with(TA_Down4,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(value = "24", location = ph_location_type(type = "sldNum")) %>%
    #ph_with(value = TA40img, ph_location_type("body",position_right = TRUE)) %>%
    ph_with(value = FIG_18, location = ph_location(label = "my_name",
                                                   left = 5.2, top = 5.5, width = 4, height = 1.4))%>%
    ph_with(my_pres, value = "Elison Timm 2017; See Annex II", location = ph_location_type(type = "dt"))%>% 
    ph_with(value = TA40img, location = ph_location(label = "my_name",
                                                     left = 4.6, top = 2, width = 5.17, height = 3.51)) %>%
    
    
    #Slide 25 #Summary
    add_slide("Title and Content","Office Theme") %>%
    ph_with(S24_TIT,   ph_location_type("title",position_left = TRUE)) %>%
    ph_with(fp_Summary,ph_location_type("body"))%>%
    #ph_with(value = "25", location = ph_location_type(type = "sldNum"))%>%
 
  
    #Slide 26#Respirces 
    add_slide("Title and Content","Office Theme") %>%
    ph_with(S25_TIT, ph_location_type("title",position_left = TRUE)) %>%
    ph_with(fp_RES, ph_location_type("body"))%>%
    ph_with(value = "26", location = ph_location_type(type = "sldNum"))%>%
    ph_with(value =  DPlanimg , location = ph_location(label = "my_name",
                                                 left = 7.3, top = 1.9, width = 2.3, height = 2.83))%>%
  
 
    #Slide 27 #Acknowledgements
    add_slide("Title and Content","Office Theme") %>%
    ph_with(S26_TIT, ph_location_type("title",position_left = TRUE)) %>%
    ph_with(fp_Ack, ph_location_type("body"))%>%
    ph_with(value = "27", location = ph_location_type(type = "sldNum"))%>%
    ph_with(value = EWCimg, location = ph_location(label = "my_name",
                                                   left = 0.85, top = 3.5, width = 1.4, height = 1.4)) %>%
    ph_with(value = CASCimg, location = ph_location(label = "my_name",
                                                    left = 2.25, top = 3.7, width = 1, height = 1)) %>%
    ph_with(value = FSimg, location = ph_location(label = "my_name",
                                                  left = 3.35, top = 3.7, width = 1, height = 1)) %>%
  
    # Replace RISA with PIDP
    ph_with(value = RISAimg, location = ph_location(label = "my_name",
                                                    left = 4.45, top = 3.7, width = 1.4, height = 1)) %>%
    ph_with(value = WRRCimg, location = ph_location(label = "my_name",
                                                    left = 5.9, top = 3.7, width = 1.1, height = 1)) %>%
    ph_with(value = NOAAimg, location = ph_location(label = "my_name",
                                                    left = 7.05, top = 3.7, width = 1, height = 1)) %>%
    ph_with(value = UHimg, location = ph_location(label = "my_name",
                                                  left = 8.15, top = 3.7, width = 1, height = 1)) %>%
    ph_with(value = CITA, location = ph_location(label = "my_name",
                                                   left = 0.85, top = 5.1, width = 8, height = 1.4))%>%
    ph_with(value = Draft, location = ph_location(label = "my_name",
                                                 left = 0.85, top = 6.6, width = 8, height = 0.75))%>%
    ph_with(value = PDKE_S, location = ph_location(label = "my_name",
                                                left = 0, top = 0, width = 1.5, height = 1.5))%>%
 
  
  
    #Slide 28  ########## References 

    add_slide("Title Only","Office Theme") %>%
    ph_with(S27_TIT,  ph_location_type("title",position_left = TRUE)) %>%
    ph_with(value = Worksimg, location = ph_location(label = "my_name",
                                                  left = 1.5, top = 1.5, width = 7, height = 6)) %>%
    ph_with(value = "28", location = ph_location_type(type = "sldNum"))%>%
    
    
    
    #Slide 29  ############## Annex I
    
    add_slide("Two Content","Office Theme") %>%
    ph_with("Annex I: 100-Year Rainfall ",ph_location_type("title")) %>%
    ph_with(fp_RF100,        ph_location_type("body",position_right = FALSE)) %>%
    #ph_with(value = "29", location = ph_location_type(type = "sldNum")) %>%
    ph_with(value = SCAimg, ph_location_type("body",position_right = TRUE)) %>%
    ph_with(value = FIG_A1, location = ph_location(label = "my_name",
                                                   left = 5, top = 6.3, width = 4.8, height = 1.4))%>%
    ph_with(my_pres, value = "Frazier et al., 2016; Lucas et al., (In Review)", location = ph_location_type(type = "dt"))%>% 
    
    
    #Slide 30  ############## Annex II
    
    add_slide("Title and Content","Office Theme") %>%
    ph_with("Annex II: Climate Downscaling in Hawaii ",ph_location_type("title")) %>%
    ph_with(fp_DownA,        ph_location_type("body",position_right = FALSE)) %>%
    ph_with(value = "30", location = ph_location_type(type = "sldNum")) %>%
    ph_with(value = Downimg, location = ph_location(label = "my_name",
                                                    left = 3.5, top = 5, width = 3, height = 2)) %>%
  
    #Slide 31 ############### Annex III

    add_slide("Title and Content","Office Theme") %>%
    ph_with("Annex III: Drought Events (1920 - 2019)",ph_location_type("title")) %>%
    ph_with(value = "31", location = ph_location_type(type = "sldNum"))%>% 
    ph_with(value = DrotT_T1, ph_location_type("body")) %>%
    ph_with(value = TAB3, location = ph_location(label = "my_name",
                                                 left = 7.85, top = 2.5, width = 2, height = 4)) %>%
    
    {if(Drot_ct>12) add_slide(.,"Title and Content", "Office Theme") %>% 
        ph_with("Annex III: Drought Events (1920 - 2019)",ph_location_type("title")) %>%
        ph_with(value = "32", location = ph_location_type(type = "sldNum"))%>%
        ph_with(value = DrotT_T2, ph_location_type("body"))%>%
        ph_with(value = TAB3,
                location = ph_location(label = "my_name",
                                       left = 7.85, top = 2.5, width = 2, height = 4)) else .} %>%
    
    {if(Drot_ct>24) add_slide(.,"Title and Content", "Office Theme") %>% 
        ph_with("Annex III: Drought Events (1920 - 2019)",ph_location_type("title")) %>%
        ph_with(value = "33", location = ph_location_type(type = "sldNum"))%>%
        ph_with(value = DrotT_T3, ph_location_type("body"))%>%
        ph_with(value = TAB3,
                location = ph_location(label = "my_name",
                                       left = 7.85, top = 2.5, width = 2, height = 4)) else .} %>%

  
  #Slide 32 (or whatever the last slide is)  ############## PDKE 
  
  add_slide("Title and Content","Office Theme") %>%
  ph_with(value = PDKE_S, location = ph_location(label = "my_name",
                                                 left = 2, top = 1, width = 6, height = 6))%>%

  
  print(mypowerpoint, target = paste0(P_FOLDER,NameF,"_CCVD_Portfolio",RFUnit2,".pptx"))
  
  

    

  
  
  
  
  
  
  
  
  
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
