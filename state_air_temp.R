###### Air Temperature Graphs for CCVD Portfolio #####

library(rgdal)
library(raster)

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


