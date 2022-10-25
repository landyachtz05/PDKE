###### Drought Bar Graphs ######

# use the "Drought History" spreadsheet created in MINI_Phase2
setwd("E:/PDKE/CCVD/MINI_Phase2/Hawaiian Homelands - Kahikinui/")

dat<-read.csv("Hawaiian Homelands - Kahikinui Drought History.csv")
dat

# make dataframe with one year per row and duration, intensity and mag columns
dat2<-data.frame(year = 1920:2019)
dat2$dur<-NA
dat2$int<-NA
dat2$mag<-NA

head(dat2)

for (x in 1:nrow(dat)) {

  t<-dat[x,]

# get start year
sy<-as.numeric(substr(t$Start,1,4))

# add values to dataframe
dat2[which(dat2$year == sy),]$dur<-t$Duration
dat2[which(dat2$year == sy),]$int<-t$P.Intensity
dat2[which(dat2$year == sy),]$mag<-t$Magnitude

}

dat2$date<-as.Date(paste0(dat2$year,"-01-01"))
dat2

##### Plots
library(ggplot2)
library(scales)

# Duration
ggplot(dat2, aes(x=date, y=dur)) +
  geom_bar(position="dodge", stat="identity", width=600, fill="darkgreen", color="black") +
  geom_smooth(method=lm, se=F, color="black") +
  scale_x_date(date_breaks="10 years", labels=date_format(format="%Y"), 
              limits=c(as.Date("1920-01-01"), as.Date("2019-01-01")), expand=c(0,0)) +
  labs(title="Drought Duration in Months", x="Year", y="Months") +
  theme_bw()

# Magnitude
ggplot(dat2, aes(x=date, y=mag)) +
  geom_bar(position="dodge", stat="identity", width=600, fill="darkred", color="black") +
  geom_smooth(method=lm, se=F, color="black") +
  scale_x_date(date_breaks="10 years", labels=date_format(format="%Y"), 
               limits=c(as.Date("1920-01-01"), as.Date("2019-01-01")), expand=c(0,0)) +
  labs(title="Drought Magnitude", x="Year", y="Magnitude (SPI)") +
  theme_bw()

  
