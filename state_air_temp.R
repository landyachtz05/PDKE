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
write.csv(table, "kahikinui_airtemp.csv")


################################################################################

##### Monthly air temp for whole time series

setwd("E:/PDKE/CCVD/air_temp")

dat<-read.csv("kahikinui_airtemp.csv")


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



################################################################################
##### Annual air temperature trend (graph with trend line)

library(ggplot2)
library(scales)
library(ggpubr)

setwd("E:/PDKE/CCVD/air_temp")

dat<-read.csv("kahikinui_airtemp.csv")

head(dat)
tail(dat)

# convert all celcius to farenheit
dat$min<-(dat$min*(9/5)) + 32
dat$max<-(dat$max*(9/5)) + 32
dat$mean<-(dat$mean * (9/5)) + 32

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

################################################################################
##### Monthly averages and deviation from those normals #####

#make copy
dat2<-dat
head(dat2)

month.avg<-aggregate(mean ~ month, dat2, mean)
month.avg

# make month name column
month.avg$month2<-forcats::fct_reorder(month.abb[month.avg$month],1:12)

# plot monthly averages
ggplot(month.avg, aes(x=month2, y=meanF, group=1)) +
  geom_line(size=1.5, color="red") +
  labs(title="Average Monthly Air Temperature",
       y="Temperature (F)", x="") +
  theme_bw()

# add avg and difference columns to dat2
dat2$avg<-0
dat2$diff<-0

# for each month-year, calculate the deviation from month average
for (x in 1:nrow(dat2)) {

  # get row
  my<-dat2[x,]
  # get month
  month<-my$month
  # get month average
  avg<-month.avg[which(month.avg$month==month),]$meanF
  # write to dat2
  dat2[x,]$avg<-avg
  # subtract average from month-year value
  diff<-my$meanF - avg
  # write to dat2
  dat2[x,]$diff<-diff
  
}

head(dat2,20)

# get date column
dat2$date<-as.Date(paste0(dat2$year,"-",dat2$month,"-01"))

# plot difference time series
ggplot(dat2, aes(x=date,y=diff)) +
  geom_line() +
  geom_smooth(span=0.2, se=F, size=1.3) +
  geom_smooth(method=lm, se=F, color="red") +
  stat_cor(method="pearson", label.x=as.Date("1993-01-01"), label.y=2.5) +
  scale_x_date(date_breaks = "4 years", labels = date_format(format="%Y")) +
  labs(title="Monthly Deviation from Normal Air Temperature",
       y="Deviation (F)", x="") +
  geom_hline(yintercept=0) +
  theme_bw()

# calculate annual average deviations
dat3<-aggregate(diff ~ year, dat2, mean)
dat3

ggplot(dat3, aes(x=year,y=diff, group=1)) +
  geom_line() 


# make date column....again
dat3$date<-as.Date(paste0(dat3$year,"-01-01"))

ggplot(dat3, aes(x=date, y=diff, group=1)) +
  geom_line(size=1.5, color="red") +
  geom_smooth(method=lm, se=F, size=1.2) +
  stat_cor(method="pearson", label.x=as.Date("1993-01-01"), label.y=0.9) +
  scale_x_date(date_breaks = "4 years", labels = date_format(format="%Y")) +
  labs(title="Annual Deviation from Normal Air Temperature",
       y="Deviation (F)", x="") +
  theme_bw()

################################################################################
#####
