###### Air Temperature Graphs for CCVD Portfolio #####

library(rgdal)
library(raster)

################################################################################
##### Extract monthly air temp stats (min, max, mean) from study area and
##### write to CSV

setwd("E:/PDKE/CCVD/")

# study area
sa <- readOGR("HHL_kahikinui_maui_WGS.dbf")

# set monthly air temp maps directory
maps<-"E:/PDKE/CCVD/air_temp/statewide/data_map/"

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
    m2 <- crop(map, extent(sa))
    m3 <- mask(m2, sa)
    
    # calculate min, max, mean air temp
    table[no,]$min<-round(cellStats(m3, stat="min"), digits=2)
    table[no,]$max<-round(cellStats(m3, stat="max"), digits=2)
    table[no,]$mean<-round(cellStats(m3, stat="mean"), digits=2)
    
    no<-no+1
  }
}

table

# write to csv
setwd("E:/PDKE/CCVD/air_temp")
write.csv(table, "hawaii_airtemp.csv")


################################################################################
##### Air temperature trend (graph with trendline)

setwd("E:/PDKE/CCVD/air_temp")

dat<-read.csv("hawaii_airtemp.csv")

head(dat)
tail(dat)

RF_IN <- as.numeric(dat[,4]) 
summary(RF_IN)

dat$Day <- 1
dat$Ndate <- as.Date(with( dat, paste(Year, Month, Day,sep="-")), "%Y-%m-%d")

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

short.date_M = strftime(dat$Ndate, "%Y/%m")

Mean_M_RF = aggregate(as.numeric(dat$RF) ~ short.date_M, FUN = mean)
colnames(Mean_M_RF) <- c("Date","RF")
short.date_Y = strftime(as.Date(dat$Ndate), "%Y")

Mean_Y_RF = aggregate(as.numeric(dat$RF) ~ short.date_Y, FUN = mean)
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
WET_RF <- subset(dat, Month==c('01','02','03','04','11','12'))
head(WET_RF)
WET_RF2 <- rbind(c(NA,"12",NA,NA,NA), WET_RF)
WET_RF3 <- rbind(c(NA,"11",NA,NA,NA), WET_RF2)
WET_RF4 <- as.numeric(WET_RF3$RF)

# WET_RF<-dat[dat$Month !="5" & dat$Month !="6" & dat$Month != "7" & 
#                  dat$Month != "8" & dat$Month != "9" & dat$Month != "10",]
# WET_RF5<-as.numeric(WET_RF$RF)

#Get seasonal average 
WET_RF5 <-  as.vector(tapply(WET_RF4, gl(length(WET_RF4)/6, 6), mean,na.rm=T))
DRY_RF <-dat[dat$Month != "01" & dat$Month  != "02" & 
                  dat$Month  != "03" & dat$Month  != "04" & 
                  dat$Month  != 11 & dat$Month  != 12,]  

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
