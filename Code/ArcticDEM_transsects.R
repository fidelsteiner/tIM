################################################################################
# Transect validation of ArcticDEM ice margin 
#
# ArcticDEM_transsects.R
#
# ReadMe: 
#
# https://github.com/fidelsteiner/tGISM
#
# Input:
# 
#
# Created:          2021/10/15
# Latest Revision:  2024/02/16
#
# Jakob F Steiner | jakob.steiner@uni-graz.at | fidelsteiner.github.io 
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

# packages (if not installed yet: install.packages('examplePackage')
library(rgdal)
library(rgeos)
library(maptools)
library(raster)
library(ggplot2)
library(dplyr)
library(ggridges)
library(RColorBrewer)
library(terra)
library(sf)
library(ncdf4)
library(plotly)
##################
# Multiple subsites for analysis of the margin
##################

data_basepath <- 'E:\\GeospatialData'                     # path for heavy input data
data_ArcticDEM <- 'E:\\ArcticDEM'                         # path for ArcticDEM data
output_basepath <- 'C:\\Work\\Research\\tGISM'            # path for code/output   

##################
# Paths
##################

path_outlines <- 'E:\\Research\\OtherResearch\\Greenland\\Code\\SectorMargins\\ManualRedRock' # Not used, remove
path_figs <- paste(output_basepath,'\\Figures',sep = "")
pathDEM_Pleiades <- paste(data_basepath,'DEMData\\ValidationDEMs\\RedRock_Pleiades\\Pleiades\\DEM\\2017-08-15_RedRockCliff\\PGO\\2017-08-15_RedRockCliff',sep = "")
pathDEM_Arctic <- data_ArcticDEM
data_Bedmachine <- 'E:\\GeospatialData\\Bedmachine\\5000000451838\\160281892' # path Bedmachine
data_lakes <- 'E:\\GeospatialData\\Lakes\\Greenland' # path Lakes
folderSectors <- paste(output_basepath,'\\Sectors',sep = "")


projec <-'+proj=utm +datum=WGS84'

projec_utm <-'+proj=utm +zone=19 +north +datum=WGS84'
projec_27N <-'+proj=utm +zone=27N +datum=WGS84'

# Validation transsects
validationTS <- 'C:\\Work\\Research\\tGISM\\MarginValidation\\ArcticDEMValidation\\Transects'
SectorCount <- list.files(validationTS)

for(k in 1:length(SectorCount)){
SecFolder <- strsplit(list.files(validationTS)[k],'S')
secNum <- as.numeric(SecFolder[[1]][1])


DEMfiles <- unlist(read.table(paste('C:\\Work\\Research\\tGISM\\Sectors\\',secNum,'Sector\\',secNum,'ArcticDEMTiles.txt',sep='')))

allTs <- list.files(paste(validationTS,'\\',secNum,'Sector',sep=''),pattern = '*.shp', full.names=T)
popOgr <- ogrInfo(paste(output_basepath,'\\Sectors\\SectorBuffers_nuna\\buffertempmarginSector',secNum,'.shp',sep = ""))
margin_buf <- readOGR(dsn=paste(output_basepath,'\\Sectors\\SectorBuffers_nuna',sep = ""),layer = popOgr$layer)
margin_con <- SpatialPolygons(margin_buf@polygons,proj4string=margin_buf@proj4string)

for(i in 1 :length(allTs)){
popOgr <- ogrInfo(allTs[i])
trans_long <- readOGR(dsn=allTs[i],layer = popOgr$layer)
trans_long <- spTransform(trans_long,CRSobj = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
trans_short <- gIntersection(trans_long,margin_con)

DEMmerged <- raster(res = 2, crs = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
DEMveclist <- vector()
for(dk in 1:length(DEMfiles)){
  demdummy <- raster(DEMfiles[dk])
  if(gIntersects(trans_long,as(extent(demdummy), "SpatialPolygons"))){
    #Sys.sleep()
    DEMveclist <- c(DEMveclist, dk)
  }
}
for(dk in 1:length(DEMveclist)){
  if(dk==1){DEMmerged <- crop(raster(DEMfiles[DEMveclist[dk]]),extent(trans_long))
  }else{DEMmerged <- merge(DEMmerged,crop(raster(DEMfiles[DEMveclist[dk]]),extent(trans_long)))}
  
}



transall <- extract(DEMmerged,trans_long,cellnumbers = T)[[1]]
if(secNum==46){
  if(i==3|i==4){
    transall[,1] <- rev(transall[,1])
    transall[,2] <- rev(transall[,2])
  }
}
if(secNum==19){
  if(i==2){
    transall[,1] <- rev(transall[,1])
    transall[,2] <- rev(transall[,2])
  }
}
if(secNum==61){
  if(i==4){
    transall[,1] <- rev(transall[,1])
    transall[,2] <- rev(transall[,2])
  }
}
if(secNum==131){
  if(i==1|i==2){
    transall[,1] <- rev(transall[,1])
    transall[,2] <- rev(transall[,2])
  }
}
if(secNum==137){
  if(i==1|i==2){
    transall[,1] <- rev(transall[,1])
    transall[,2] <- rev(transall[,2])
  }
}
if(secNum==301|secNum==304){
    transall[,1] <- rev(transall[,1])
    transall[,2] <- rev(transall[,2])
}
transall_coords <- xyFromCell(DEMmerged, transall[,1]) 

df <- data.frame(id = c(1), x_origin =transall_coords[1,1], y_origin = transall_coords[1,2], x_dest =transall_coords[length(transall_coords[,1]),1], y_dest=transall_coords[length(transall_coords[,1]),2]) 
origin <- st_as_sf(df, coords = c("x_origin", "y_origin"), crs='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' )
dest <- st_as_sf(df, coords = c("x_dest", "y_dest"), crs='+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs' )
transect_distances <- c(0,st_distance(origin, dest,  by_element = TRUE ))
transect_steps <- seq(1, transect_distances[2], transect_distances[2]/length(transall[,2]))

transmargin <- extract(DEMmerged,trans_short,cellnumbers = T)
idMatch <- match(coordinates(transmargin)[,1],coordinates(transall)[,1])
coords <- cbind(coordinates(DEMmerged)[coordinates(transmargin)[,1],1],coordinates(DEMmerged)[coordinates(transmargin)[,1],2])
transmarg <- unlist(extract(DEMmerged,trans_long)) * NA
transmarg[idMatch] <- unlist(extract(DEMmerged,trans_short))

png(file=paste('C:\\Work\\Research\\tGISM\\MarginValidation\\ArcticDEMValidation\\Figures\\',secNum,'_',i,'.png',sep=''), res = 300,width=1800,height=800) 
par(mar = c(4, 4, 1, 1))
plot(transect_steps,transall[,2],type='l',lwd=2,xlab ='distance [m]',ylab = 'Elevation [m a.s.l.]',cex = 2,xaxs ='i',xaxt='n',yaxt='n')
points(transect_steps,transmarg, type='l', col='red',lwd=3)
axis(1, at = seq(0,max(transect_steps),100))
axis(2, at = seq(round(min(transall[,2])),round(max(transall[,2])),20))
abline(v=seq(0,max(transect_steps),100),h=seq(round(min(transall[,2])),round(max(transall[,2])),20),
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.5)      # Grid line width
dev.off()
}
}