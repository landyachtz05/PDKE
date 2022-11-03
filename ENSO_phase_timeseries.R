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

  
#################################################################################
###### Assign ENSO phases to season-years (1990 dry, 1990 wet, 1991 dry, 1991 wet, etc.)

# combine ENSO phase and rainfall tables (1950 - 2012)
table<-read.csv("statewide_monthly_rainfall_1920_2012.csv")

# add date column if not already there
table$date<-as.Date(paste0(table$year,"-",table$month,"-01"))
tail(table)
summary(table$rainfall)

head(table)
head(enso)

table2<-merge(enso, table, by="date")
head(table2, 20)

# create season column and fill using month
table2$season<-NA

for (x in 1:nrow(table2)) {
  
  r<-table2[x,]
  
  if(r$month < 5 || r$month > 10) {table2[x,]$season <- "wet"}
  if(r$month < 11 && r$month > 4) {table2[x,]$season <- "dry"}
}





### for each year aggregate delta_t (ONI value) by season and keep maximum

#list of years
year<-as.list(unique(table2$year))

# create empty dataframe
seasons<-setNames(data.frame(matrix(ncol = 4, nrow = 0)), 
                  c("year", "season", "delta_t", "s.phase"))
head(seasons)

# loop through years and aggregate season ONI values
for (y in year) {
  
  b<-table2[which(table2$year == y),]
  
  # find min and max ONI values for each season-year
  h<-aggregate(delta_t ~ season, b, max)
  l<-aggregate(delta_t ~ season, b, min)
  
  # for each season, determine which value is larger (absolute) and assign ENSO phases
  seas<-as.list(unique(table2$season))
  
  for (s in seas) {
    
    h2<-h[which(h$season == s),]
    l2<-l[which(l$season == s),]
    l2
    
    if(abs(h2$delta_t) > abs(l2$delta_t)) {seasons[nrow(seasons) + 1,]$delta_t <- h2$delta_t}
    if(abs(h2$delta_t) < abs(l2$delta_t)) {seasons[nrow(seasons) + 1,]$delta_t <- l2$delta_t}
    
    if(seasons[nrow(seasons),]$delta_t>1.5) {seasons[nrow(seasons),]$s.phase <- "SEL"}
    if(seasons[nrow(seasons),]$delta_t<=1.5 & seasons[nrow(seasons),]$delta_t>=0.5) {seasons[nrow(seasons),]$s.phase <- "WEL"}
    if(seasons[nrow(seasons),]$delta_t>(-0.5) & seasons[nrow(seasons),]$delta_t<0.5) {seasons[nrow(seasons),]$s.phase <- "NUT"}
    if(seasons[nrow(seasons),]$delta_t<=(-0.5) & seasons[nrow(seasons),]$delta_t>=(-1.5)) {seasons[nrow(seasons),]$s.phase <- "WLA"}
    if(seasons[nrow(seasons),]$delta_t<(-1.5)) {seasons[nrow(seasons),]$s.phase <- "SLA"}
    
    seasons[nrow(seasons),]$year<-y
    seasons[nrow(seasons),]$season<-s
  }
}

head(seasons, 10)

### add column of rainfall values (sum of season-year)
# calculate sum of each season-year rainfall
head(table2)

rain.seas<-aggregate(rainfall ~ year+season, table2, sum)

colnames(rain.seas)[3]<-"rain_sum"

head(rain.seas)

# add to seasonal phase dataset
seasons2<-merge(seasons, rain.seas, all=T)
seasons2

write.csv(seasons2, "statewide_ENSO_phase_seasonyear_rainfall.csv")


################################################################################
### barplot

seasons3<-seasons2
head(seasons3)

# remove delta_t column from seasons2 so it isn't used in the join
seasons4<-subset(seasons3, select=-c(delta_t))
head(seasons4)

# add s.phase and rain_sum columns
library(dplyr)
rain4<-left_join(table2, seasons4)

rain4$month_no<-as.factor(rain4$month_no)
head(rain4, 20)

# # If there's "NA" in the rain_sum or s.phase columns
# for (i in 1:nrow(rain4)) {
#   if(is.na(rain4[i,]$rain_sum)) {rain4[i,]$rain_sum<-rain4[i-1,]$rain_sum}
#   if(is.na(rain4[i,]$s.phase)) {rain4[i,]$s.phase<-rain4[i-1,]$s.phase}
#   
# }
# 
head(rain4,100)

### See if all months are represented by each ENSO phase
# list of enso phases
p<-unique(rain4$phase5)

p

d2<-data.frame()

for (i in p) {
  
  d<-setNames(data.frame(matrix(ncol = 2, nrow = 1)),
              c("phase","months"))
  
  d$phase<-i
  
  sub<-rain4[which(rain4$phase5 == i),]
  
  m<-as.numeric(unique(sub$month_no))
  
  s<-sum(unique(m))
  
  if(s == 78) {d$months<-c("ALL")}
  if(s != 78) {d$months<-c("NOT ALL")}
  
  d2<-rbind(d2, d)
}

d2

### count how many seasons are in each ENSO phase
# make season.year column
rain4$season.year<-paste0(rain4$year,"_",rain4$season)


rain4<-rain4[order(rain4$phase5, rain4$Year, rain4$month_no),]
rain4<-rain4[order(rain4$month_year),]

head(rain4)

# make list of phases
phase<-as.list(unique(rain4$s.phase))

# make list of seasons
s<-as.list(unique(rain4$season))

dat_all<-data.frame()

for (i in phase) {
  
  # subset for each phase
  sub<-rain4[which(rain4$s.phase == i),]
  head(sub)
  
  for(x in s) {
    
    # subset for each season
    sub2<-sub[which(sub$season == x),]
    head(sub2)
    
    # count unique season.years
    sy<-data.frame(as.numeric(length(unique(sub2$season.year))))
    sy$s.phase<-i
    sy$season<-x
    colnames(sy)<-c("sy.count","s.phase","season")
    sy
    
    # add to dataframe
    dat_all<-rbind(dat_all, sy)
  }
}

dat_all
# 
# # set count values
# sel.ct<-dat_all[4,1]
# wel.ct<-dat_all[1,1]
# nut.ct<-dat_all[5,1]
# wla.ct<-dat_all[2,1]
# sla.ct<-dat_all[3,1]


# aggregate over each s.phase and season
rain5<-aggregate(rain_sum ~ s.phase + pval5 + season, rain4, FUN=mean)
rain5

# add season.year count column
rain6<-full_join(rain5, dat_all)
rain6

# plot
ggplot(data=rain6, 
       aes(x=season, y=rain_sum, group=reorder(s.phase, pval5))) +
  geom_bar(aes(fill=s.phase), position = "dodge", stat="identity", color="black", 
           alpha=.7, width=0.7) +
  labs(title="Average Seasonal Rainfall by ENSO Phase",
       y = "Rainfall (mm)", x= "Season") +
  scale_fill_manual(values=c("darkred", "darkorange1", "darkgoldenrod1", "Dark Green", "blue3"),
                    limits=c("SEL","WEL","NUT","WLA","SLA")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2500)) +
  guides(fill=guide_legend(title="ENSO Phase")) +
  geom_text(aes(label=sy.count), position=position_dodge(width=0.7), vjust=-0.8) +
  theme_bw()










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