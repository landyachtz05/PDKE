################################################################################
### Standardized Precipitation Index (SPI)
# https://search.r-project.org/CRAN/refmans/precintcon/html/spi.html

setwd("E:/PDKE/CCVD/MINI_Phase2/Haleakala National Park/")

# load monthly precipitation data with 4 columns (Date, Year, Month, RF)
dat<-read.csv("Haleakala National Park Monthly Rainfall_in.csv")
head(dat,20)

# calculate SPI and SPEI
library(SPEI)

spi<-spi(ts(dat$RF, freq=12, start=c(1920, 1)), scale = 12)

# write spi to dataframe
spi_val<-spi$fitted
spi_val<-data.frame(spi=as.matrix(spi_val), date=time(spi_val))
spi_val$date<-sub("\\..*","",spi_val$date)
head(spi_val, 20)

# fix date column
library(zoo)
spi_val$date2<-rep(c(1:12))
spi_val$date<-as.Date(paste0(spi_val$date,"/",spi_val$date2, "/01"), format = "%Y/%m/%d")

spi_val<-spi_val[1:2]
colnames(spi_val)<-c("spi_12","date")

head(spi_val,20)
tail(spi_val)
summary(spi_val$spi_12)

# see strongest drought months
head(spi_val[order(spi_val$spi_12),])

plot(spi, main = "Haleakala 12-Month SPI")

plot(spei, main="Haleakala 12-month SPEI")

write.csv(spi_val, "haleakala_spi_12month.csv")

################################################################################
### Create inverted (negatives only) SPI graph like fig. 11 in MKWP CCVD portfolio
# (https://www.eastwestcenter.org/sites/default/files/filemanager/Research_Program/PDKE%20Page/Mauna%20Kahalawai%20Watershed%20Partnership_CCVD_Portfolio_V2_Workingin.pdf)
setwd("E:/PDKE/")

# load SPI csv from above
spi_full<-read.csv("haleakala_spi_12month.csv")
spi_full<-spi_full[c("spi_12","date")]
head(spi_full)
summary(spi_full$spi_12)

# format date column
spi_full$date<-as.Date(spi_full$date)

# sort by date and make a new copy
spi2<-spi_full[order(as.Date(spi_full$date)),]
head(spi2,20)
tail(spi2)
summary(spi2$spi_12)

# plot SPI
library(ggplot2)
library(scales)

ggplot(spi2, aes(x=date, y=spi_12)) +
  geom_area(fill = "blue") +
  geom_line() +
  scale_x_date(date_breaks = "2 years", labels = date_format(format = "%Y")) +
  ylim(-3,3) +
  labs(x = "Year", y = "Drought Intensity", title = "Guam SPI-12") +
  geom_hline(yintercept = -1, linetype = "dashed", color = "orange", size = 1)


# Change all positive values to 0, and convert negative values to positive
spi2$spi_negs<-ifelse(spi2$spi_12>0,0,spi2$spi_12)
spi2$spi_negs<-abs(spi2$spi_negs)

summary(spi2$spi_negs)
head(spi2, 50)

# plot spi_negs
ggplot(spi2, aes(x=date, y=spi_negs)) +
  geom_area(fill = "orange") +
  geom_line() +
  scale_x_date(date_breaks = "5 years", labels = date_format(format = "%Y")) +
  ylim(0,3) +
  labs(x = "Year", y = "Drought Intensity", title = "Guam SPI-12 Drought Events") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "orange", size = 1) +
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "brown", size = 1)


################################################################################
### count number of drought events and calculate stats (duration, peak, magnitude, mean)

  # create binary column of drought (spi_negs>0) yes or no
  spi2$drought<-ifelse(spi2$spi_negs>0, 1, 0)
  head(spi2, 100)

  # make column of consecutive drought month count
  library(data.table)

  spi2$DRGHT_ct<-with(spi2, (drought == 1)*
                        ave(drought, rleid(drought == 1),
                            FUN=seq_along))

  summary(spi2$DRGHT_ct)
  
  # subset data for drought event months only
  spi3<-spi2[which(spi2$DRGHT_ct >= 1),]
  head(spi3, 20)

  ### create another empty event count column and fill
  spi3$event_ct<-0
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
    
    # fill event_Ct2 column with final event number
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
  summary(spi5$event_ct2)
  
  # get rid of intermediate event column
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
      
      SPI_I$months<-nrow(SPI_I)
      
      if(max(SPI_I$spi_negs) <= 1.5) {SPI_I$intensity<-1}
      if(max(SPI_I$spi_negs) > 1.5 && max(SPI_I$spi_negs) <= 2) {SPI_I$intensity<-2}
      if(max(SPI_I$spi_negs) > 2) {SPI_I$intensity<-3}

      SPI_I$peak<-max(SPI_I$spi_negs)
      SPI_I$mean<-mean(SPI_I$spi_negs)
      SPI_I$mag<-sum(SPI_I$spi_negs)

      # keep only these columns
      SPI_I<-SPI_I[c("date","months","intensity","peak","mean","mag")]
      
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
      
      # make empty table
      event.t<-data.frame(matrix(ncol=8, nrow=0))
      x<-c("event_ct","start", "stop", "months","intensity", "peak", "mean", "mag")
      colnames(event.t)<-x
      event.t$start<-as.POSIXct(event.t$start)
      event.t$stop<-as.POSIXct(event.t$stop)
      event.t
      
      ### loop through rows and fill table
      
      for (y in 1:nrow(spi6)) {
        
        # get row
        row<-spi6[y,]
        
        # get previous and next row
        prev<-spi6[y-1,]
        nex<-spi6[y+1,]
        
        # get event number
        number<-row$event_ct
        
        if(!is.na(row$intensity) && is.na(prev$intensity)) {event.t[number,]$start<-nex$date}
        if(!is.na(row$intensity) && is.na(prev$intensity)) {event.t[number,]$intensity<-row$intensity}
        if(!is.na(row$intensity) && is.na(nex$intensity)) {event.t[number,]$stop<-nex$date}
        if(!is.na(row$intensity)) {event.t[number,]$months<-row$months}
        if(!is.na(row$intensity)) {event.t[number,]$event_ct<-row$event_ct}
        if(!is.na(row$intensity)) {event.t[number,]$peak<-row$peak}
        if(!is.na(row$intensity)) {event.t[number,]$mean<-row$mean}
        if(!is.na(row$intensity)) {event.t[number,]$mag<-row$mag}
      }
      
      # fix date formats
      library(zoo)
      event.t$start<-as.yearmon(sub(" .*","", event.t$start), "%Y-%m-%d")
      event.t$stop<-as.yearmon(sub(" .*","", event.t$stop), "%Y-%m-%d")
      
      event.t
      
      # get count of event intensities
      table(event.t$intensity)
  
      write.csv(event.t, "haleakala_spi_drought_events_derek2.csv")
    
### END ###
 
      
     
##################################
### EGWAS ###

# Calculate site's min, max, mean SPI over entire rainfall record (1920 - 2022)

# ranch name
ranch<-"Z Bar Ranch"
setwd(paste0("E:/PDKE/CCVD/MINI_Phase2/",ranch,"/"))

spi12<-read.csv(paste0(ranch,"SPI_12.csv"))
head(spi12, 20)

min<-min(spi12$SP, na.rm=T)
mean<-mean(spi12$SP, na.rm=T)
max<-max(spi12$SP, na.rm=T)

### Make monthly min, mean, max SPI value dataset

# get month and year columns
spi12$year<-substr(spi12$DT,1,4)
spi12$month<-as.numeric(substr(spi12$DT,6,7))

# remove NA rows
spi12<-spi12[which(!is.na(spi12$SP)),]

# aggregate over months
spi12mean<-sapply(split(spi12$SP,spi12$month),mean)
spi12min<-sapply(split(spi12$SP,spi12$month),min)
spi12max<-sapply(split(spi12$SP,spi12$month),max)
spi12mean

spi12m<-rbind(spi12mean,spi12min,spi12max)
spi12m

write.csv(spi12m, paste0(ranch,"SPI_12 monthly.csv"))

# # ### seasons are wet november-march, dry april-october, so needs to start in 1921 dry season
# # spi12<-spi12[16:nrow(spi12),]
# # spi12<-spi12[which(spi12$DT>as.Date("1921-04-01")),]
# 
# spi12$season<-"na"
# 
# # make season column
# for (x in 1:nrow(spi12)) {
#   d<-spi12[x,]
#   if(d$month > 4 && d$month < 11) {spi12[x,]$season <-"dry"}
#   if(d$month > 10 | d$month < 5) {spi12[x,]$season <- "wet"}
# }
# 
# # calculate average seasonal SPI value
# dry<-mean(spi12[which(spi12$season == "dry"),]$SP)
# dry
# wet<-mean(spi12[which(spi12$season == "wet"),]$SP)
# wet
