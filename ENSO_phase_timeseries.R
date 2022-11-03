##### Calculate rainfall by ENSO phase state-wide #####

setwd("E:/PDKE/CCVD/CCVD INPUTS/")

# load ENSO phase dataset (using same ONI dataset from Guam analysis)
enso<-read.csv("enso_oni_1950_2022.csv")
head(enso)
enso$date<-as.Date(enso$date)

### plot ENSO phases over time
library(ggplot2)
library(scales)
  
  # make table for background color bars
  data_breaks<-data.frame(start=c(-2.8,-1.5,-0.5,0.5,1.5),
                          end=c(-1.5,-0.5,0.5,1.5,2.8),
                          colors=c("blue3","lightskyblue1","white","indianred1","red3"))
  
  # plot
  ggplot() +
    geom_rect(data=data_breaks, aes(ymin=start, ymax=end,
                                    xmin=as.Date("1950-01-01"), 
                                    xmax=as.Date("2022-01-01")),
              fill=data_breaks$colors, alpha=0.5) +
    geom_line(data=enso, aes(x=date, y=delta_t, group=1),size=1) +
    scale_x_date(date_breaks = "5 years", labels = date_format(format="%Y"),
                 expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    labs(title="ENSO Phases based on Sea Surface Temperature (SST)",
         y="Change in SST (°F)", x="Year") +
    geom_hline(yintercept=0) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, vjust=0.65, size=12),
          axis.text.y=element_text(size=12))

    
  
  
  
### process statewide rainfall maps into monthly values
setwd("E:/PDKE/CCVD/data/production/rainfall/legacy/month/statewide/data_map/")

  # loop through year folders and calculate average value for each month

  # list of years
  years<-as.list(list.files())
  years
  
  # create empty table for raster stats
  table<-data.frame(matrix(ncol=3,nrow=0, 
                           dimnames=list(NULL, c("year","month", "rainfall"))))
  
  # set row
  no<-1
  
  ### for each year, loop through months and extract mean rainfall value
  library(raster)
  
  for (i in years) {
  
    setwd(paste0("E:/PDKE/CCVD/data/production/rainfall/legacy/month/statewide/data_map/",i,"/"))
  
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
      
      # calculate min, max, mean air temp
      table[no,]$rainfall<-round(cellStats(map, stat="mean"), digits=2)
      
      no<-no+1
    }
  }
  
  table
  
  # add date column
  table$date<-paste0(table$year,"-",table$month,"-01")
  tail(table)
  summary(table$rainfall)

  # write to csv
  setwd("E:/PDKE/CCVD/")
  write.csv(table, "statewide_monthly_rainfall_1920_2012.csv")

# combine ENSO phase and rainfall tables (1950 - 2012)
head(table)
head(enso)

table2<-merge(enso, table, by="date")
head(table2)

# calculate average rainfall value for each month and ENSO phase
table3<-aggregate(rainfall ~ month + phase5, table2, FUN=mean)
table3

table4<-aggregate(rainfall ~ phase5, table2, FUN=mean)
table4
barplot(table4$rainfall ~ table4$phase5)


# plot
library(ggplot2)

ggplot(data=table3, aes(x=month, y=rainfall, group=phase5, color=phase5)) +
  geom_line()