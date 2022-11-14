##### Calculate rainfall by ENSO phase state-wide #####

setwd("E:/PDKE/CCVD/CCVD INPUTS/")

### Load MEI dataset from CCVD Inputs
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
setwd("E:/PDKE/CCVD/MINI_Phase2/Parker Ranch/")
table<-read.csv("Parker Ranch Monthly Rainfall_in.csv")
head(table)
tail(table)

### make season column
table$season<-NA

for (i in 1:nrow(table)) {
  
  s<-table[i,]
  
  if(s$Month>4 && s$Month<11) {table[i,]$season<-"dry"}
  if(s$Month<5 || s$Month>10) {table[i,]$season<-"wet"}
}

### calculate season-year rainfall means
# subset for correct time period (1920 - 2012)
table2<-table[which(table$Year<2013),]

# remove first and last months to start dataset at beginning of 1920 dry season
table2 =  table2[-c(1115:1116),]
table2 =  table2[-c(1:4),]
head(table2, 20)
tail(table2, 20)

# add season count column
table2$sc<-NA

r<-seq(from=1, to=1114, by=6)
n<-1

for (i in r) {

  l<-i+5
  s<-table2[i:l,]
  table2[i:l,]$sc<-n
  n<-n+1
}

# aggregate rainfall over consecutive seasons
t2<-aggregate(RF ~ sc, table2, mean)
t2

### merge other columns back in and reduce to one row per season count
  t1<-table2[c(2,3,5,6)]
  head(t1, 20)
  
  # make empty dataframe
  t<-data.frame("Year")
  t$Month<-NA
  t$Season<-NA
  t$sc<-NA
  
  # loop through rows and only write row to new table
  # if sc is greater than the previous sc
  for (i in 1:nrow(t1)) {
  
    # current sc
    x<-t1[i,]$sc
    # previous sc
    n<-t1[i-1,]$sc
  
    # keep the first row (sc 1)
    if(i==1) {t[i,]<-t1[i,]}
    if(i==1) {n<-1}
    if(i>1 && x>n) {t[x,]<-t1[i,]}
  }
  
  head(t)
  
  library(m)
  t3<-left_join(t,t2)
  t3
  
  # remove month column, fix year colname
  t3<-t3[c(1,3,4,5)]
  colnames(t3)[which(names(t3) == "X.Year.")] <- "Year"

#################################################################################
###### Assign ENSO phases to season-years (1990 dry, 1990 wet, 1991 dry, 1991 wet, etc.)

# combine ENSO phase and rainfall tables
head(t3)
head(enso)

# add year column to enso for join
enso$Year<-substr(enso$date, 1, 4)

library(plyr)
dfs<-list(t3, enso)
table3<-join_all(dfs, match="all")
head(table3, 20)
tail(table3,20)

# remove rows with no data (past 2012)
table3<-table3[which(!is.na(table3$MEI_D)),]

# make single MEI and s.phase columns with value based on season
table3$MEI<-NA
table3$x<-NA

for (i in 1:nrow(table3)) {
  
  s<-table3[i,]
  if(s$Season == "wet") {table3[i,]$MEI<-s$MEI_W}
  if(s$Season == "dry") {table3[i,]$MEI<-s$MEI_D}
  if(s$Season == "wet") {table3[i,]$x<-s$phase.w}
  if(s$Season == "dry") {table3[i,]$x<-s$phase.d}
}

# remove unneeded MEI and phase columns
table3<-table3[c(1,2,3,4,10,11)]
head(table3, 20)

write.csv(table3, "Kahikinui_ENSO_phase_seasonyear_meanrainfall_MEI.csv")

################################################################################
################################################################################
### Boxplots of seasonal rainfall by ENSO phase
### rainy season (July-December)
library(ggplot2)

unique(table3$x)
table3<-table3[which(!is.na(table3$x)),]

### 5 phases
# assign values to ENSO phases
table3$pval5<-ifelse(table3$x=="SEL",1,0)
table3$pval5<-ifelse(table3$x=="WEL",2,table3$pval5)
table3$pval5<-ifelse(table3$x=="NUT",3,table3$pval5)
table3$pval5<-ifelse(table3$x=="WLA",4,table3$pval5)
table3$pval5<-ifelse(table3$x=="SLA",5,table3$pval5)
summary(table3$pval5)

head(table3)

########################################################
### keep only wet season months
season.w<-table3[which(table3$Season == "wet"),]
head(season.w)

### stats
# count
library(dplyr)
count<-count(season.w$x)
count
# count<-season.w %>% count(s.phase, sort=T)

c.sel<-as.numeric(count[which(count$x == "SEL"),]$freq)
c.sel<-paste0("count = ",c.sel)
c.wel<-as.numeric(count[which(count$x == "WEL"),]$freq)
c.wel<-paste0("count = ",c.wel)
c.nut<-as.numeric(count[which(count$x == "NUT"),]$freq)
c.nut<-paste0("count = ",c.nut)
c.wla<-as.numeric(count[which(count$x == "WLA"),]$freq)
c.wla<-paste0("count = ",c.wla)
c.sla<-as.numeric(count[which(count$x == "SLA"),]$freq)
c.sla<-paste0("count = ",c.sla)

# min, mean, max rainfall
stats<-summary(season.w[which(season.w$x == "SEL"),])
stats
sel.min<-stats[1,4]
sel.min<-sub(".*:","",sel.min)
# sel.min<-sub(" .*","",sel.min)
sel.min<-paste0("min = ",sel.min)

sel.mean<-stats[4,4]
sel.mean<-sub(".*:","",sel.mean)
# sel.mean<-sub(" .*","",sel.mean)
sel.mean<-paste0("mean = ",sel.mean)

sel.max<-stats[6,4]
sel.max<-sub(".*:","",sel.max)
# sel.max<-sub(" .*","",sel.max)
sel.max<-paste0("max = ",sel.max)

stats<-summary(season.w[which(season.w$x == "WEL"),])

wel.min<-stats[1,4]
wel.min<-sub(".*:","",wel.min)
# wel.min<-sub(" .*","",wel.min)
wel.min<-paste0("min =",wel.min)

wel.mean<-stats[4,4]
wel.mean<-sub(".*:","",wel.mean)
# wel.mean<-sub(" .*","",wel.mean)
wel.mean<-paste0("mean = ",wel.mean)

wel.max<-stats[6,4]
wel.max<-sub(".*:","",wel.max)
# wel.max<-sub(" .*","",wel.max)
wel.max<-paste0("max = ",wel.max)

stats<-summary(season.w[which(season.w$x == "NUT"),])

nut.min<-stats[1,4]
nut.min<-sub(".*:","",nut.min)
# nut.min<-sub(" .*","",nut.min)
nut.min<-paste0("min =",nut.min)

nut.mean<-stats[4,4]
nut.mean<-sub(".*:","",nut.mean)
# nut.mean<-sub(" .*","",nut.mean)
nut.mean<-paste0("mean = ",nut.mean)

nut.max<-stats[6,4]
nut.max<-sub(".*:","",nut.max)
# nut.max<-sub(" .*","",nut.max)
nut.max<-paste0("max = ",nut.max)

stats<-summary(season.w[which(season.w$x == "WLA"),])

wla.min<-stats[1,4]
wla.min<-sub(".*:","",wla.min)
# wla.min<-sub(" .*","",wla.min)
wla.min<-paste0("min =",wla.min)

wla.mean<-stats[4,4]
wla.mean<-sub(".*:","",wla.mean)
# wla.mean<-sub(" .*","",wla.mean)
wla.mean<-paste0("mean = ",wla.mean)

wla.max<-stats[6,4]
wla.max<-sub(".*:","",wla.max)
# wla.max<-sub(" .*","",wla.max)
wla.max<-paste0("max = ",wla.max)

stats<-summary(season.w[which(season.w$x == "SLA"),])

sla.min<-stats[1,4]
sla.min<-sub(".*:","",sla.min)
# sla.min<-sub(" .*","",sla.min)
sla.min<-paste0("min = ",sla.min)

sla.mean<-stats[4,4]
sla.mean<-sub(".*:","",sla.mean)
# sla.mean<-sub(" .*","",sla.mean)
sla.mean<-paste0("mean = ",sla.mean)

sla.max<-stats[6,4]
sla.max<-sub(".*:","",sla.max)
# sla.max<-sub(" .*","",sla.max)
sla.max<-paste0("max = ",sla.max)

# get average rainfall value across all phases
avg.w<-mean(season.w$RANCH_RF)
avg.w<-round(avg.w, digits=2)

# create average rainfall value label
label<-paste("average rainfall = ",avg.w, sep="")

# set ylims
ycount<-max(season.w$RANCH_RF*1.2)
ymin<-max(season.w$RANCH_RF*1.15)
ymean<-max(season.w$RANCH_RF*1.1)
ymax<-max(season.w$RANCH_RF*1.05)
avg<-max(season.w$RANCH_RF*1.3)

ggplot(season.w, aes(x=reorder(x, pval5), y=RANCH_RF)) +
  geom_boxplot() +
  # ylim(100,2950) +
  labs(title="Seasonal rainfall by ENSO category (rainy season)",
       y = "Seasonal Rainfall (mm)", x= "ENSO catogory") +
  geom_hline(yintercept=avg.w, linetype="dashed", color="blue", size=1) +
  
  annotate("text",x="SEL",y=ycount,label=c.sel) +
  annotate("text",x="WEL",y=ycount,label=c.wel) +
  annotate("text",x="NUT",y=ycount,label=c.nut) +
  annotate("text",x="WLA",y=ycount,label=c.wla) +
  annotate("text",x="SLA",y=ycount,label=c.sla) +
  
  annotate("text", x="WEL", y=avg, label=label) +
  
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
season.d<-table3[which(table3$Season == "dry"),]
head(season.d)

# # sort by mean
# season.d<-season.d[with(season.d, order(-rain_sum)),]

### stats
# count
library(dplyr)
count<-count(season.d$x)
count

c.sel<-as.numeric(count[2,]$freq)
c.sel<-paste0("count = ",c.sel)

c.wel<-as.numeric(count[4,]$freq)
c.wel<-paste0("count = ",c.wel)

c.nut<-as.numeric(count[1,]$freq)
c.nut<-paste0("count = ",c.nut)

c.wla<-as.numeric(count[5,]$freq)
c.wla<-paste0("count = ",c.wla)

c.sla<-as.numeric(count[3,]$freq)
c.sla<-paste0("count = ",c.sla)

# min, mean, max rainfall
stats<-summary(season.d[which(season.d$x == "SEL"),])
stats
sel.min<-stats[1,4]
sel.min<-sub(".*:","",sel.min)
# sel.min<-sub(" .*","",sel.min)
sel.min<-paste0("min = ",sel.min)

sel.mean<-stats[4,4]
sel.mean<-sub(".*:","",sel.mean)
# sel.mean<-sub(" .*","",sel.mean)
sel.mean<-paste0("mean = ",sel.mean)

sel.max<-stats[6,4]
sel.max<-sub(".*:","",sel.max)
# sel.max<-sub(" .*","",sel.max)
sel.max<-paste0("max = ",sel.max)

stats<-summary(season.d[which(season.d$x == "WEL"),])

wel.min<-stats[1,4]
wel.min<-sub(".*:","",wel.min)
# wel.min<-sub(" .*","",wel.min)
wel.min<-paste0("min =",wel.min)

wel.mean<-stats[4,4]
wel.mean<-sub(".*:","",wel.mean)
# wel.mean<-sub(" .*","",wel.mean)
wel.mean<-paste0("mean = ",wel.mean)

wel.max<-stats[6,4]
wel.max<-sub(".*:","",wel.max)
# wel.max<-sub(" .*","",wel.max)
wel.max<-paste0("max = ",wel.max)

stats<-summary(season.d[which(season.d$x == "NUT"),])

nut.min<-stats[1,4]
nut.min<-sub(".*:","",nut.min)
# nut.min<-sub(" .*","",nut.min)
nut.min<-paste0("min =",nut.min)

nut.mean<-stats[4,4]
nut.mean<-sub(".*:","",nut.mean)
# nut.mean<-sub(" .*","",nut.mean)
nut.mean<-paste0("mean = ",nut.mean)

nut.max<-stats[6,4]
nut.max<-sub(".*:","",nut.max)
# nut.max<-sub(" .*","",nut.max)
nut.max<-paste0("max = ",nut.max)

stats<-summary(season.d[which(season.d$x == "WLA"),])

wla.min<-stats[1,4]
wla.min<-sub(".*:","",wla.min)
# wla.min<-sub(" .*","",wla.min)
wla.min<-paste0("min =",wla.min)

wla.mean<-stats[4,4]
wla.mean<-sub(".*:","",wla.mean)
# wla.mean<-sub(" .*","",wla.mean)
wla.mean<-paste0("mean = ",wla.mean)

wla.max<-stats[6,4]
wla.max<-sub(".*:","",wla.max)
# wla.max<-sub(" .*","",wla.max)
wla.max<-paste0("max = ",wla.max)

stats<-summary(season.d[which(season.d$x == "SLA"),])

sla.min<-stats[1,4]
sla.min<-sub(".*:","",sla.min)
# sla.min<-sub(" .*","",sla.min)
sla.min<-paste0("min = ",sla.min)

sla.mean<-stats[4,4]
sla.mean<-sub(".*:","",sla.mean)
# sla.mean<-sub(" .*","",sla.mean)
sla.mean<-paste0("mean = ",sla.mean)

sla.max<-stats[6,4]
sla.max<-sub(".*:","",sla.max)
# sla.max<-sub(" .*","",sla.max)
sla.max<-paste0("max = ",sla.max)

# get average rainfall value across all phases
avg.d<-mean(season.d$RANCH_RF)
avg.d<-round(avg.d, digits=2)

# create average rainfall value label
label<-paste("average rainfall = ",avg.d, sep="")

# set ylims
ycount<-max(season.d$RANCH_RF*1.2)
ymin<-max(season.d$RANCH_RF*1.15)
ymean<-max(season.d$RANCH_RF*1.1)
ymax<-max(season.d$RANCH_RF*1.05)
avg<-max(season.d$RANCH_RF*1.3)

ggplot(season.d, aes(x=reorder(x, pval5), y=RANCH_RF)) +
  geom_boxplot() +
  # ylim(100,2950) +
  labs(title="Seasonal rainfall by ENSO category (dry season)",
       y = "Seasonal Rainfall (mm)", x= "ENSO catogory") +
  geom_hline(yintercept=avg.w, linetype="dashed", color="blue", size=1) +
  
  annotate("text",x="SEL",y=ycount,label=c.sel) +
  annotate("text",x="WEL",y=ycount,label=c.wel) +
  annotate("text",x="NUT",y=ycount,label=c.nut) +
  annotate("text",x="WLA",y=ycount,label=c.wla) +
  annotate("text",x="SLA",y=ycount,label=c.sla) +
  
  annotate("text", x="WEL", y=avg, label=label) +
  
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

rain4<-table3
head(rain4)

# add s.phase and rain_sum columns
# library(dplyr)
# rain4<-left_join(table2, seasons4)

# rain4$month<-as.factor(rain4$month)
# head(rain4, 20)

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
rain4$season.year<-paste0(rain4$Year,"_",rain4$Season)

rain4<-rain4[order(rain4$phase5, rain4$Year, rain4$month),]
rain4<-rain4[order(rain4$month_year),]

head(rain4)

# make list of phases
phase<-as.list(unique(rain4$x))

# make list of seasons
s<-as.list(unique(rain4$Season))

dat_all<-data.frame()

for (i in phase) {
  
  # subset for each phase
  sub<-rain4[which(rain4$x == i),]
  head(sub)
  
  for(x in s) {
    
    # subset for each season
    sub2<-sub[which(sub$Season == x),]
    head(sub2)
    
    # count unique season.years
    sy<-data.frame(as.numeric(length(unique(sub2$season.year))))
    sy$x<-i
    sy$Season<-x
    colnames(sy)<-c("sy.count","x","Season")
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
rain5<-aggregate(RANCH_RF ~ x + Season, rain4, FUN=mean)
rain5

# add season.year count column
rain6<-full_join(rain5, dat_all)
rain6<-rain6[1:10,]
rain6

# set order of phases and seasons for plotting
rain6$x<-factor(rain6$x, levels=c("SEL","WEL","NUT","WLA","SLA"))
rain6$Season<-factor(rain6$Season, levels=c("wet","dry"))

# set ylim
ylim<-max(rain6$RANCH_RF)*1.2

head(rain6)
# plot
ggplot(data=rain6, 
       aes(x=Season, y=RANCH_RF, group=x)) +
  geom_bar(aes(fill=x), position = position_dodge(width=0.7), stat="identity", color="black", 
           alpha=.7, width=.55) +
  labs(title="Average Seasonal Rainfall by ENSO Phase",
       y = "Rainfall (inches)", x= "Season") +
  scale_fill_manual(values=c("red3", "darkgoldenrod1", "forestgreen", "royalblue1", "blue3"),
                    limits=c("SEL","WEL","NUT","WLA","SLA")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  guides(fill=guide_legend(title="ENSO Phase")) +
  geom_text(aes(label=sy.count), position=position_dodge(width=0.7), vjust=-0.8) +
  theme_bw()


# calculate average rainfall value for each month and ENSO phase
table3<-aggregate(rainfall ~ month + phase5, table2, FUN=mean)
table3

table4<-aggregate(RANCH_RF ~ x, rain7, FUN=mean)
table4
barplot(table4$RANCH_RF ~ table4$x)



################################################################################
################################################################################
### line plot of average monthly rainfall (January to December) by ENSO Phase
head(table2)
colnames(table2)[which(names(table2) == "season")]<-"Season"
table2<-table2[c(1:5)]

head(table3)
table3$Year<-as.numeric(table3$Year)
table3<-table3[c(1,2,3,5,6,7)]

table5<-left_join(table2, table3)
table5

# aggregate over season-phase and month
rain7<-aggregate(RANCH_RF ~ Month + x, table5, FUN=mean)
rain7

### plot in inches ###
# set order of phases for plotting
rain7$x<-factor(rain7$x, levels=c("SEL","WEL","WLA","SLA","NUT"))

# set order of months for plotting
rain7$Month<-factor(rain7$Month, levels=c(11,12,1,2,3,4,5,6,7,8,9,10))

# set ylim
ymin<-min(rain7$RANCH_RF)*0.7
ymax<-max(rain7$RANCH_RF)*1.3

ggplot(data=rain7, 
       aes(x=Month, y=RANCH_RF, group=(x))) +
  geom_line(aes(color=x), stat="identity", size=1.2) +
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
