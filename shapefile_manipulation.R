library(raster)
library(readxl)

##### Join taxonomic data to soil map unit polygon shapefile
setwd("E:/PDKE/CCVD/wss_gsmsoil_HI_[2016-10-13]/")

dat<-as.data.frame(read_excel("HI _Taxonomy_Dom_Con_02_26_2023.xlsx"))
names(dat)[names(dat) == "mukey"]<-"MUKEY"

head(dat)
dat[which(dat$MUKEY == 2371887),]

shape<-shapefile("spatial/gsmsoilmu_a_hi3_MUNAME_merge.shp")

head(shape)

m<-merge(shape, dat, by="MUKEY")

shapefile(m, "spatial/gsmsoilmu_a_hi3_MUNAME_merge2.shp")

###############################################################################

##### Extract Soil info for each EGWAS ranch
library(sf)

# load soils data
setwd("E:/PDKE/CCVD/gSSURGO_HI_2_28/")
soils<-st_read("MUPOLYGON_taxonomymerge2.shp")

# load ranch sites shapefile
setwd("E:/PDKE/CCVD/Cattle_Ranch_Individual_11_Shapefiles_Climate_Jan2023/")
sites<-st_read("EGWAS_ranches_merge.shp")

# make list of ranch sites
l<-as.list(unique(sites$Name))

# match soils crs to ranch sites
soils<-st_transform(soils, st_crs(sites))
plot(soils)

# iterate through ranch sites and extract soils data

for (i in l) {
  
  i = l[[1]]
  
  # get site polygon(s)
  s<-sites[which(sites$Name == i),]

  # clip soils using site polygon
  s2<-st_crop(soils$geometry, s$geometry)
  
}
