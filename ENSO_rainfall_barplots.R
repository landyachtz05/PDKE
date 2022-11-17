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

# load monthly rainfall
setwd("E:/PDKE/CCVD/MINI_Phase2/Parker Ranch/")
table<-read.csv("Parker Ranch Monthly Rainfall_in.csv")
colnames(table)<-c("X","Year","Month","RF")
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
  head(table2)
  tail(table2)
  
  # add season count column
  table2$sc<-NA
  
  # make sequence of numbers by 6 (1,7,13,etc.) to select the first month
  # of each season
  r<-seq(from=1, to=nrow(table2), by=6)
  n<-1

  for (i in r) {
    
    # set the last row (last month of season)
    l<-i+5
    # select all season rows
    s<-table2[i:l,]
    # assign season count number
    table2[i:l,]$sc<-n
    # next number
    n<-n+1
  }
  
  # aggregate rainfall over consecutive seasons
  t2<-aggregate(RF ~ sc, table2, sum)
  head(t2)
  
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
    
    # join in the seasonal rainfall by SC
    library(dplyr)
    t3<-left_join(t,t2)
    
    # remove month column, fix year colname (Year, Season, SC, RF)
    # t3 should now contain total season-year rainfall values (sum of each season-year)
    t3<-t3[c(1,3,4,5)]
    colnames(t3)[which(names(t3) == "X.Year.")] <- "Year"
    head(t3)

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

# remove unneeded MEI and phase columns (keep Year, Season, SC, RF, MEI, x)
table3<-table3[c(1,2,3,4,10,11)]
head(table3, 20)

# unique(table3$x)
# table3<-table3[which(!is.na(table3$x)),]

# assign values to ENSO phases
table3$pval5<-ifelse(table3$x=="SEL",1,0)
table3$pval5<-ifelse(table3$x=="WEL",2,table3$pval5)
table3$pval5<-ifelse(table3$x=="NUT",3,table3$pval5)
table3$pval5<-ifelse(table3$x=="WLA",4,table3$pval5)
table3$pval5<-ifelse(table3$x=="SLA",5,table3$pval5)
summary(table3$pval5)

head(table3)

### Make dataframe with average monthly rainfall values for each ENSO phase
  # fix year values and remove RF column
  table3$Year<-as.numeric(table3$Year)
  table4<-table3[c(1,2,3,5,6,7)]
  
  # bring in monthly rainfall values
    head(table2)
    # fix column name and remove sc column
    colnames(table2)[which(names(table2) == "season")]<-"Season"
    table2<-table2[c(1:5)]
  
  # join by year and season
  table5<-left_join(table2, table4)
  head(table5)
  
  # aggregate over month to calculate total rainfall for each ENSO phase
  rain7<-aggregate(RF ~ x, table5, FUN=mean)
  rain7

# spell out the ENSO phases
  rain7$x2<-NA
  
  for (i in 1:nrow(rain7)) {
    y<-rain7[i,]
    if(y$x == "SEL") {rain7[i,]$x2 <- "Strong El Nino"}
    if(y$x == "WEL") {rain7[i,]$x2 <- "Weak El Nino"}
    if(y$x == "NUT") {rain7[i,]$x2 <- "Neutral"}
    if(y$x == "WLA") {rain7[i,]$x2 <- "Weak La Nina"}
    if(y$x == "SLA") {rain7[i,]$x2 <- "Strong La Nina"}
  }
  
### Plot
# set order of ENSO phases for plotting
library(ggplot2)
rain7$x2<-factor(rain7$x2, levels=c("Strong El Nino","Weak El Nino","Neutral","Weak La Nina","Strong La Nina"))

# set ylim
ylim<-max(rain7$RF)*1.2

ggplot(data=rain7, 
       aes(x=x2, y=RF, group=x2)) +
  geom_bar(aes(fill=x2), position = position_dodge(width=0.65), stat="identity", color="black", 
           alpha=.7, width=.50) +
  labs(title="Average Monthly Rainfall by ENSO Phase",
       y = "Rainfall (inches)", x= "ENSO Phase") +
  scale_fill_manual(values=c("darkred","red","grey","lightskyblue1","darkblue"),
                    limits=c("Strong El Nino","Weak El Nino","Neutral","Weak La Nina","Strong La Nina")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
  guides(fill=guide_legend(title="ENSO Phase")) +
  # geom_text(aes(label=sy.count), position=position_dodge(width=0.7), vjust=-0.8) +
  theme_bw() +
  theme(axis.text.x=element_blank())


################################################################################
### Seasonal rainfall by ENSO phase barplot

rain4<-table3
head(rain4)

### count how many seasons are in each ENSO phase
# make season.year column
rain4$season.year<-paste0(rain4$Year,"_",rain4$Season)

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
head(rain4)
rain5<-aggregate(RF ~ x + Season, rain4, FUN=mean)
rain5

# add season.year count column
rain6<-full_join(rain5, dat_all)
rain6<-rain6[1:10,]
rain6

# spell out the ENSO phases
rain6$x2<-NA

for (i in 1:nrow(rain6)) {
  y<-rain6[i,]
  if(y$x == "SEL") {rain6[i,]$x2 <- "Strong El Nino"}
  if(y$x == "WEL") {rain6[i,]$x2 <- "Weak El Nino"}
  if(y$x == "NUT") {rain6[i,]$x2 <- "Neutral"}
  if(y$x == "WLA") {rain6[i,]$x2 <- "Weak La Nina"}
  if(y$x == "SLA") {rain6[i,]$x2 <- "Strong La Nina"}
}

### Plot
  # set order of seasons for plotting
  rain6$Season<-factor(rain6$Season, levels=c("wet","dry"))
  
  # set order of ENSO phases for plotting
  rain6$x2<-factor(rain6$x2, levels=c("Strong El Nino","Weak El Nino","Neutral","Weak La Nina","Strong La Nina"))
  
  # set ylim
  ylim<-max(rain6$RF)*1.2

  ggplot(data=rain6, 
         aes(x=Season, y=RF, group=x2)) +
    geom_bar(aes(fill=x2), position = position_dodge(width=0.7), stat="identity", color="black", 
             alpha=.7, width=.55) +
    labs(title="Average Seasonal Rainfall by ENSO Phase",
         y = "Rainfall (inches)", x= "Season") +
    scale_fill_manual(values=c("darkred","red","grey","lightskyblue1","darkblue"),
                      limits=c("Strong El Nino","Weak El Nino","Neutral","Weak La Nina","Strong La Nina")) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ylim)) +
    # guides(fill=guide_legend(title="ENSO Phase")) +
    geom_text(aes(label=sy.count), position=position_dodge(width=0.7), vjust=-0.8) +
    theme_bw()+
    theme(axis.text.x=element_text(size=13),
          axis.text.y=element_text(size=13),
          axis.title.x=element_text(size=13),
          axis.title.y=element_text(size=13),
          legend.position="none")


  

