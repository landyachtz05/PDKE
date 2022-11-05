##### Calculate rainfall by ENSO phase state-wide #####

setwd("E:/PDKE/CCVD/CCVD INPUTS/")

# load ENSO phase dataset (using same ONI dataset from Guam analysis)
enso<-read.csv("enso_oni_1950_2022.csv")
enso

enso$date<-as.Date(enso$date)

### OOOORRRR load MEI dataset from CCVD Inputs
enso<-read.csv("MEI_Season.csv")
enso

# Add date column
enso$date<-c(1920:2012)
enso$date<-as.Date(paste0(enso$date,"-01-01"))

# add seasonal ENSO phase columns based on MEI values
enso$phase.w<-NA
enso$phase.d<-NA

# assign 5 ENSO phase based on MEI_W
enso$phase.w<-ifelse(enso$MEI_W>1.5,"SEL",0)
enso$phase.w<-ifelse(enso$MEI_W>=0.5 & enso$MEI_W<=1.5, "WEL", enso$phase.w)
enso$phase.w<-ifelse(enso$MEI_W>(-0.5) & enso$MEI_W<0.5, "NUT", enso$phase.w)
enso$phase.w<-ifelse(enso$MEI_W<=(-0.5) & enso$MEI_W>=(-1.5), "WLA", enso$phase.w)
enso$phase.w<-ifelse(enso$MEI_W<(-1.5), "SLA", enso$phase.w)
    
# assign 5 ENSO phase based on MEI_D
enso$phase.d<-ifelse(enso$MEI_D>1.5,"SEL",0)
enso$phase.d<-ifelse(enso$MEI_D>=0.5 & enso$MEI_D<=1.5, "WEL", enso$phase.d)
enso$phase.d<-ifelse(enso$MEI_D>(-0.5) & enso$MEI_D<0.5, "NUT", enso$phase.d)
enso$phase.d<-ifelse(enso$MEI_D<=(-0.5) & enso$MEI_D>=(-1.5), "WLA", enso$phase.d)
enso$phase.d<-ifelse(enso$MEI_D<(-1.5), "SLA", enso$phase.d)

head(enso, 20)

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
                                    xmin=as.Date("1920-01-01"), 
                                    xmax=as.Date("2012-01-01")),
              fill=data_breaks$colors, alpha=0.5) +
    geom_line(data=enso, aes(x=date, y=MEI_W, group=1), size=1, color="blue") +
    geom_line(data=enso, aes(x=date, y=MEI_D, group=1), size=1, color="red3") +
    scale_x_date(date_breaks = "5 years", labels = date_format(format="%Y"),
                 expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    labs(title="ENSO Phases based on Sea Surface Temperature (SST)",
         y="Change in SST (°F)", x="Year") +
    geom_hline(yintercept=0) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, vjust=0.65, size=12),
          axis.text.y=element_text(size=12))

# load Kahikinui monthly rainfall
setwd("E:/PDKE/CCVD/MINI_Phase2/Hawaiian Homelands - Kahikinui/")
table<-read.csv("Hawaiian Homelands - Kahikinui Monthly Rainfall_in.csv")
head(table)

#################################################################################
###### Assign ENSO phases to season-years (1990 dry, 1990 wet, 1991 dry, 1991 wet, etc.)

# combine ENSO phase and rainfall tables (1950 - 2012)

# add date column if not already there
table$date<-as.Date(paste0(table$Date,"/01"), format="%Y/%m/%d")
head(table)
summary(table$RF)

head(table)
head(enso)

# add month-year date column to enso and remove full date column so it isn't used
# in the join
library(stringr)
enso$Year<-str_sub(enso$date, end=-7)
enso<-subset(enso, select=-c(date))

library(plyr)
dfs<-list(table, enso)
table2<-join_all(dfs, match="all")
table2

### make season column
table2$season<-NA

for (i in 1:nrow(table2)) {
  
  s<-table2[i,]
  
  if(s$Month>4 && s$Month<11) {table2[i,]$season<-"dry"}
  if(s$Month<5 || s$Month>10) {table2[i,]$season<-"wet"}
}

# make MEI column 

# loop through years and aggregate season ONI values
for (y in year) {
  
  b<-table2[which(table2$Year.y == y),]
  
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

# fix table2 year column name
colnames(table2)[which(names(table2) == "Year.y")] <- "year"

### add column of rainfall values (sum of season-year)
# calculate sum of each season-year rainfall
head(table2)

rain.seas<-aggregate(RF ~ year+season, table2, sum)

colnames(rain.seas)[3]<-"rain_sum"

rain.seas

seasons

# add to seasonal phase dataset
seasons2<-merge(seasons, rain.seas, all=T)
seasons2

write.csv(seasons2, "Kahikinui_ENSO_phase_seasonyear_rainfall.csv")





################################################################################
################################################################################
### Boxplots of seasonal rainfall by ENSO phase
### rainy season (July-December)
library(ggplot2)

head(seasons2)
unique(seasons2$s.phase)

# remove NA seasons
seasons2<-seasons2[which(!is.na(seasons2$s.phase)),]

# # get months column
# rain2$month<-as.numeric(substring(rain2$month_year, 6))

### 3 phases
seasons2$pval<-ifelse(seasons2$s.phase=="SEL" | seasons2$s.phase == "WEL",1,0)
seasons2$pval<-ifelse(seasons2$s.phase=="NUT", 2, seasons2$pval)
seasons2$pval<-ifelse(seasons2$s.phase=="WLA" | seasons2$s.phase == "SLA",3, seasons2$pval)
summary(seasons2$pval)

### 5 phases
# assign values to ENSO phases
seasons2$pval5<-ifelse(seasons2$s.phase=="SEL",1,0)
seasons2$pval5<-ifelse(seasons2$s.phase=="WEL",2,seasons2$pval5)
seasons2$pval5<-ifelse(seasons2$s.phase=="NUT",3,seasons2$pval5)
seasons2$pval5<-ifelse(seasons2$s.phase=="WLA",4,seasons2$pval5)
seasons2$pval5<-ifelse(seasons2$s.phase=="SLA",5,seasons2$pval5)
summary(seasons2$pval5)

head(seasons2)

########################################################
### keep only wet season months
season.w<-subset(seasons2, season == "wet")
head(season.w)

### stats
# count
library(dplyr)
count<-season.w %>% count(s.phase, sort=T)
count
c.sel<-as.numeric(count[4,]$n)
c.sel<-paste0("count = ",c.sel)

c.wel<-as.numeric(count[2,]$n)
c.wel<-paste0("count = ",c.wel)

c.nut<-as.numeric(count[1,]$n)
c.nut<-paste0("count = ",c.nut)

c.wla<-as.numeric(count[3,]$n)
c.wla<-paste0("count = ",c.wla)

c.sla<-as.numeric(count[5,]$n)
c.sla<-paste0("count = ",c.sla)

# min, mean, max
stats<-summary(season.w[which(season.w$s.phase == "SEL"),])

sel.min<-stats[1,5]
sel.min<-sub(".*:","",sel.min)
sel.min<-sub(" .*","",sel.min)
sel.min<-paste0("min = ",sel.min)

sel.mean<-stats[4,5]
sel.mean<-sub(".*:","",sel.mean)
sel.mean<-sub(" .*","",sel.mean)
sel.mean<-paste0("mean = ",sel.mean)

sel.max<-stats[6,5]
sel.max<-sub(".*:","",sel.max)
sel.max<-sub(" .*","",sel.max)
sel.max<-paste0("max = ",sel.max)

stats<-summary(season.w[which(season.w$s.phase == "WEL"),])

wel.min<-stats[1,5]
wel.min<-sub(".*:","",wel.min)
wel.min<-sub(" .*","",wel.min)
wel.min<-paste0("min = ",wel.min)

wel.mean<-stats[4,5]
wel.mean<-sub(".*:","",wel.mean)
wel.mean<-sub(" .*","",wel.mean)
wel.mean<-paste0("mean = ",wel.mean)

wel.max<-stats[6,5]
wel.max<-sub(".*:","",wel.max)
wel.max<-sub(" .*","",wel.max)
wel.max<-paste0("max = ",wel.max)

stats<-summary(season.w[which(season.w$s.phase == "NUT"),])

nut.min<-stats[1,5]
nut.min<-sub(".*:","",nut.min)
nut.min<-sub(" .*","",nut.min)
nut.min<-paste0("min = ",nut.min)

nut.mean<-stats[4,5]
nut.mean<-sub(".*:","",nut.mean)
nut.mean<-sub(" .*","",nut.mean)
nut.mean<-paste0("mean = ",nut.mean)

nut.max<-stats[6,5]
nut.max<-sub(".*:","",nut.max)
nut.max<-sub(" .*","",nut.max)
nut.max<-paste0("max = ",nut.max)

stats<-summary(season.w[which(season.w$s.phase == "WLA"),])

wla.min<-stats[1,5]
wla.min<-sub(".*:","",wla.min)
wla.min<-sub(" .*","",wla.min)
wla.min<-paste0("min = ",wla.min)

wla.mean<-stats[4,5]
wla.mean<-sub(".*:","",wla.mean)
wla.mean<-sub(" .*","",wla.mean)
wla.mean<-paste0("mean = ",wla.mean)

wla.max<-stats[6,5]
wla.max<-sub(".*:","",wla.max)
wla.max<-sub(" .*","",wla.max)
wla.max<-paste0("max = ",wla.max)

stats<-summary(season.w[which(season.w$s.phase == "SLA"),])

sla.min<-stats[1,5]
sla.min<-sub(".*:","",sla.min)
sla.min<-sub(" .*","",sla.min)
sla.min<-paste0("min = ",sla.min)

sla.mean<-stats[4,5]
sla.mean<-sub(".*:","",sla.mean)
sla.mean<-sub(" .*","",sla.mean)
sla.mean<-paste0("mean = ",sla.mean)

sla.max<-stats[6,5]
sla.max<-sub(".*:","",sla.max)
sla.max<-sub(" .*","",sla.max)
sla.max<-paste0("max = ",sla.max)

# get average rainfall value across all phases
avg.w<-mean(season.w$rain_sum)
avg.w<-round(avg.w, digits=2)

# create average rainfall value label
label<-paste("average rainfall = ",avg.w, sep="")

# set ylims
ycount<-2850
ymin<-2750
ymean<-2650
ymax<-2550

ggplot(season.w, aes(x=reorder(s.phase, pval5), y=rain_sum)) +
  geom_boxplot() +
  ylim(100,2950) +
  labs(title="Seasonal rainfall by ENSO category (rainy season)",
       y = "Seasonal Rainfall (mm)", x= "ENSO catogory") +
  geom_hline(yintercept=avg.w, linetype="dashed", color="blue", size=1) +
  
  annotate("text",x="SEL",y=ycount,label=c.sel) +
  annotate("text",x="WEL",y=ycount,label=c.wel) +
  annotate("text",x="NUT",y=ycount,label=c.nut) +
  annotate("text",x="WLA",y=ycount,label=c.wla) +
  annotate("text",x="SLA",y=ycount,label=c.sla) +
  
  annotate("text", x="WEL", y=250, label=label) +
  
  annotate("text",x="SEL",y=ymin,label=sel.min) +
  annotate("text",x="SEL",y=ymean,label=sel.mean) +
  annotate("text",x="SEL",y=ymax,label=sel.max) +
  
  annotate("text",x="WEL",y=ymin,label=wel.min) +
  annotate("text",x="WEL",y=ymean,label=wel.mean) +
  annotate("text",x="WEL",y=ymax,label=wel.max) +
  
  annotate("text",x="NUT",y=ymin,label=nut.min) +
  annotate("text",x="NUT",y=ymean,label=nut.mean) +
  annotate("text",x="NUT",y=ymax,label=nut.max) +
  
  annotate("text",x="WLA",y=ymin,label=wla.min) +
  annotate("text",x="WLA",y=ymean,label=wla.mean) +
  annotate("text",x="WLA",y=ymax,label=wla.max) +
  
  annotate("text",x="SLA",y=ymin,label=sla.min) +
  annotate("text",x="SLA",y=ymean,label=sla.mean) +
  annotate("text",x="SLA",y=ymax,label=sla.max)

######################################################
### keep only dry season months
season.d<-subset(seasons2, season == "dry")
head(season.d)

season.d[which(season.d$pval5 == 1),]

# sort by mean
season.d<-season.d[with(season.d, order(-rain_sum)),]

### stats
# count
library(dplyr)
count<-season.d %>% count(s.phase, sort=T)
count
c.sel<-as.numeric(count[4,]$n)
c.sel<-paste0("count = ",c.sel)

c.wel<-as.numeric(count[2,]$n)
c.wel<-paste0("count = ",c.wel)

c.nut<-as.numeric(count[1,]$n)
c.nut<-paste0("count = ",c.nut)

c.wla<-as.numeric(count[3,]$n)
c.wla<-paste0("count = ",c.wla)

c.sla<-as.numeric(count[5,]$n)
c.sla<-paste0("count = ",c.sla)


# min, mean, max
stats<-summary(season.d[which(season.d$s.phase == "SEL"),])

sel.min<-stats[1,5]
sel.min<-sub(".*:","",sel.min)
sel.min<-sub(" .*","",sel.min)
sel.min<-paste0("min = ",sel.min)

sel.mean<-stats[4,5]
sel.mean<-sub(".*:","",sel.mean)
sel.mean<-sub(" .*","",sel.mean)
sel.mean<-paste0("mean = ",sel.mean)

sel.max<-stats[6,5]
sel.max<-sub(".*:","",sel.max)
sel.max<-sub(" .*","",sel.max)
sel.max<-paste0("max = ",sel.max)

stats<-summary(season.d[which(season.d$s.phase == "WEL"),])
stats

wel.min<-stats[1,5]
wel.min<-sub(".*: ","",wel.min)
wel.min<-sub(" .*","",wel.min)
wel.min<-paste0("min = ",wel.min)

wel.mean<-stats[4,5]
wel.mean<-sub(".*: ","",wel.mean)
wel.mean<-sub(" .*","",wel.mean)
wel.mean<-paste0("mean = ",wel.mean)

wel.max<-stats[6,5]
wel.max<-sub(".*:","",wel.max)
wel.max<-sub(" .*","",wel.max)
wel.max<-paste0("max = ",wel.max)

stats<-summary(season.d[which(season.d$s.phase == "NUT"),])

nut.min<-stats[1,5]
nut.min<-sub(".*: ","",nut.min)
nut.min<-sub(" .*","",nut.min)
nut.min<-paste0("min = ",nut.min)

nut.mean<-stats[4,5]
nut.mean<-sub(".*: ","",nut.mean)
nut.mean<-sub(" .*","",nut.mean)
nut.mean<-paste0("mean = ",nut.mean)

nut.max<-stats[6,5]
nut.max<-sub(".*:","",nut.max)
nut.max<-sub(" .*","",nut.max)
nut.max<-paste0("max = ",nut.max)

stats<-summary(season.d[which(season.d$s.phase == "WLA"),])

wla.min<-stats[1,5]
wla.min<-sub(".*: ","",wla.min)
wla.min<-sub(" .*","",wla.min)
wla.min<-paste0("min = ",wla.min)

wla.mean<-stats[4,5]
wla.mean<-sub(".*: ","",wla.mean)
wla.mean<-sub(" .*","",wla.mean)
wla.mean<-paste0("mean = ",wla.mean)

wla.max<-stats[6,5]
wla.max<-sub(".*:","",wla.max)
wla.max<-sub(" .*","",wla.max)
wla.max<-paste0("max = ",wla.max)

stats<-summary(season.d[which(season.d$s.phase == "SLA"),])

sla.min<-stats[1,5]
sla.min<-sub(".*: ","",sla.min)
sla.min<-sub(" .*","",sla.min)
sla.min<-paste0("min = ",sla.min)

sla.mean<-stats[4,5]
sla.mean<-sub(".*:","",sla.mean)
sla.mean<-sub(" .*","",sla.mean)
sla.mean<-paste0("mean = ",sla.mean)

sla.max<-stats[6,5]
sla.max<-sub(".*:","",sla.max)
sla.max<-sub(" .*","",sla.max)
sla.max<-paste0("max = ",sla.max)



# calculate mean rainfall across phases
avg.d<-mean(season.d$rain_sum)
avg.d<-round(avg.d, digits=2)

# set ylims
ycount<-2850
ymin<-2750
ymean<-2650
ymax<-2550

# create average rainfall value label
label<-paste("average rainfall = ",avg.d, sep="")


ggplot(season.d, aes(x=reorder(s.phase, pval5), y=rain_sum)) +
  geom_boxplot() +
  ylim(100, 2950) +
  labs(title="Seasonal rainfall by ENSO category (dry season)",
       y = "Seasonal rainfall (mm)", x= "ENSO catogory") +
  geom_hline(yintercept=avg.d, linetype="dashed", color="blue", size=1) +
  
  annotate("text",x="SEL",y=ycount,label=c.sel) +
  annotate("text",x="WEL",y=ycount,label=c.wel) +
  annotate("text",x="NUT",y=ycount,label=c.nut) +
  annotate("text",x="WLA",y=ycount,label=c.wla) +
  annotate("text",x="SLA",y=ycount,label=c.sla) +
  
  annotate("text",x="WLA", y=250, label=label) +
  
  annotate("text",x="SEL",y=ymin,label=sel.min) +
  annotate("text",x="SEL",y=ymean,label=sel.mean) +
  annotate("text",x="SEL",y=ymax,label=sel.max) +
  
  annotate("text",x="WEL",y=ymin,label=wel.min) +
  annotate("text",x="WEL",y=ymean,label=wel.mean) +
  annotate("text",x="WEL",y=ymax,label=wel.max) +
  
  annotate("text",x="NUT",y=ymin,label=nut.min) +
  annotate("text",x="NUT",y=ymean,label=nut.mean) +
  annotate("text",x="NUT",y=ymax,label=nut.max) +
  
  annotate("text",x="WLA",y=ymin,label=wla.min) +
  annotate("text",x="WLA",y=ymean,label=wla.mean) +
  annotate("text",x="WLA",y=ymax,label=wla.max) +
  
  annotate("text",x="SLA",y=ymin,label=sla.min) +
  annotate("text",x="SLA",y=ymean,label=sla.mean) +
  annotate("text",x="SLA",y=ymax,label=sla.max)







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

rain4$month<-as.factor(rain4$month)
head(rain4, 20)

# # If there's "NA" in the rain_sum or s.phase columns
# for (i in 1:nrow(rain4)) {
#   if(is.na(rain4[i,]$rain_sum)) {rain4[i,]$rain_sum<-rain4[i-1,]$rain_sum}
#   if(is.na(rain4[i,]$s.phase)) {rain4[i,]$s.phase<-rain4[i-1,]$s.phase}
#   
# }
# 
head(rain4,100)

# ### See if all months are represented by each ENSO phase
# # list of enso phases
# p<-unique(rain4$phase5)
# p
# 
# d2<-data.frame()
# 
# for (i in p) {
# 
#   d<-setNames(data.frame(matrix(ncol = 2, nrow = 1)),
#               c("phase","months"))
#   d$phase<-i
#   sub<-rain4[which(rain4$phase5 == i),]
#   m<-as.numeric(unique(sub$month))
#   s<-sum(unique(m))
#   if(s == 78) {d$months<-c("ALL")}
#   if(s != 78) {d$months<-c("NOT ALL")}
#   d2<-rbind(d2, d)
# }
# 
# d2

### count how many seasons are in each ENSO phase
# make season.year column
rain4$season.year<-paste0(rain4$year,"_",rain4$season)

rain4<-rain4[order(rain4$phase5, rain4$Year, rain4$month),]
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
rain5<-aggregate(rain_sum ~ s.phase + season, rain4, FUN=mean)
rain5

# add season.year count column
rain6<-full_join(rain5, dat_all)
rain6<-rain6[1:10,]
rain6

# set order of phases and seasons for plotting
rain6$s.phase<-factor(rain6$s.phase, levels=c("SEL","WEL","NUT","WLA","SLA"))
rain6$season<-factor(rain6$season, levels=c("wet","dry"))

# set ylim
ylim<-max(rain6$rain_sum)*1.2

# plot
ggplot(data=rain6, 
       aes(x=season, y=rain_sum, group=s.phase)) +
  geom_bar(aes(fill=s.phase), position = "dodge", stat="identity", color="black", 
           alpha=.7, width=0.7) +
  labs(title="Average Seasonal Rainfall by ENSO Phase",
       y = "Rainfall (mm)", x= "Season") +
  scale_fill_manual(values=c("red3", "darkgoldenrod1", "forestgreen", "royalblue1", "blue3"),
                    limits=c("SEL","WEL","NUT","WLA","SLA")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  guides(fill=guide_legend(title="ENSO Phase")) +
  geom_text(aes(label=sy.count), position=position_dodge(width=0.7), vjust=-0.8) +
  theme_bw()


# calculate average rainfall value for each month and ENSO phase
table3<-aggregate(rainfall ~ month + phase5, table2, FUN=mean)
table3

table4<-aggregate(rainfall ~ phase5, table2, FUN=mean)
table4
barplot(table4$rainfall ~ table4$phase5)



################################################################################
################################################################################
### line plot of average monthly rainfall (January to December) by ENSO Phase
head(rain4)

# aggregate over season-phase and month
rain7<-aggregate(rainfall ~ month + s.phase, rain4, FUN=mean)
rain7[order(rain7$s.phase),]
rain7

### plot in inches ###
# set order of phases for plotting
rain7$s.phase<-factor(rain7$s.phase, levels=c("SEL","WEL","WLA","SLA","NUT"))

# set order of months for plotting
rain7$month<-factor(rain7$month, levels=c(11,12,1,2,3,4,5,6,7,8,9,10))

# set ylim
ymin<-min(rain7$rainfall)*0.7
ymax<-max(rain7$rainfall)*1.3

ggplot(data=rain7, 
       aes(x=month, y=rainfall, group=(s.phase))) +
  geom_line(aes(color=s.phase), stat="identity", size=1.2) +
  # geom_smooth(method="lm", formula=y~poly(x,9), se=F) +
  labs(title="Average Monthly Rainfall by ENSO Phase",
       y = "Rainfall (in.)", x= "Month") +
  scale_color_manual(values=c("red3", "darkgoldenrod1", "forestgreen", "royalblue1", "blue3"),
                     limits=c("SEL","WEL","NUT","WLA","SLA")) +
  scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax)) +
  guides(color=guide_legend(title="ENSO Phase")) +
  geom_vline(xintercept = 6.5, color="black") +
  annotate(geom="text", x=3.5, y=ymax*0.93, label = "wet season", color="blue", size=5) +
  annotate(geom="text", x=9.5, y=ymax*0.93, label = "dry season", color = "darkred", size=5) +
  theme_bw() +
  theme(axis.title.y=element_text(size=14), 
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.title.x=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12))
