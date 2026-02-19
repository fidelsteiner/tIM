################################################################################
# Transect validation of ArcticDEM ice margin 
#
# ArcticDEM_transsects.R
#
# ReadMe: Extracts DEM data for polylines placed across the margin to investigate match of the 
# actual ice margin with ice masks. Requires polylines placed for the study provided in the repository.
# To perform validation at additional sites code needs to be adapted.
# List of ArcticDEM tile names associated to sectors (provided as txt files in the repository) need to be available
#
# https://github.com/fidelsteiner/tIM
#
# Input:
# (1) numbered lines saved as shp files placed across the actual ice margin
# (2) actual ice masks
# (3) list of ArcticDEM tiles associated to each Sector number
# (4) ArcticDEM tiles (not provided in repositors, need to be procured from https://arcticdem.apps.pgc.umn.edu/)
# 
#
# Created:          2021/10/15
# Latest Revision:  2026/02/13
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
# Paths (to be adapted for local machine)
##################

validationTS <- 'C:\\Work\\Research\\tGISM\\MarginValidation\\ArcticDEMValidation\\Transects' # path with actual transect files stored in folders named by 'XXSector'
buffLocations <- 'C:\\Work\\Research\\tGISM\\Sectors\\SectorBuffers_nuna'
pathDEM_Arctic <- 'E:\\ArcticDEM'         # path for ArcticDEM data
output_basepath <- 'E:\\tGISM'            # path for any code/output   
path_figs <- paste(output_basepath,'\\Figures',sep = "")

##################
# Processing of transects
##################

SectorCount <- list.files(validationTS)

for(k in 1:length(SectorCount)){
SecFolder <- strsplit(list.files(validationTS)[k],'S')
secNum <- as.numeric(SecFolder[[1]][1])

DEMfiles <- unlist(read.table(paste(output_basepath,'\\Sectors\\',secNum,'Sector\\',secNum,'ArcticDEMTiles.txt',sep='')))

allTs <- list.files(paste(validationTS,'\\',secNum,'Sector',sep=''),pattern = '*.shp', full.names=T)
popOgr <- ogrInfo(paste(buffLocations,'\\buffertempmarginSector',secNum,'.shp',sep = ""))
margin_buf <- readOGR(dsn=paste(buffLocations,sep = ""),layer = popOgr$layer)
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

png(file=paste(path_figs,'\\',secNum,'_',i,'.png',sep=''), res = 300,width=1800,height=800) 
par(mar = c(4, 4, 1, 1))
plot(transect_steps,transall[,2],type='l',lwd=2,xlab ='distance [m]',ylab = 'Elevation [m a.s.l.]',cex = 2,xaxs ='i',xaxt='n',yaxt='n',xlim=c(20,max(transect_steps)-100),ylim=c(floor(min(transall[,2])/10)*10,floor(min(transall[,2])/10)*10 + 60))
points(transect_steps,transmarg, type='l', col='red',lwd=3)
axis(1, at = seq(100,max(transect_steps),100))
axis(2, at = seq(floor(min(transall[,2])/10)*10,floor(min(transall[,2])/10)*10 + 60,20))
abline(v=seq(0,max(transect_steps),100),h=seq(round(min(transall[,2])),round(max(transall[,2])),20),
     lty = 2,      # Grid line type
     col = "gray", # Grid line color
     lwd = 0.5)      # Grid line width
dev.off()

if(secNum==311){
  if(i==2){
    png(file=paste(path_figs,'\\',secNum,'_',i,'.png',sep=''), res = 300,width=1800,height=800) 
    par(mar = c(4, 4, 1, 1))
    plot(transect_steps,transall[,2],type='l',lwd=2,xlab ='distance [m]',ylab = 'Elevation [m a.s.l.]',cex = 2,xaxs ='i',xaxt='n',yaxt='n',xlim=c(20,max(transect_steps)-100),ylim=c(360,360 + 60))
    points(transect_steps,transmarg, type='l', col='red',lwd=3)
    axis(1, at = seq(100,max(transect_steps),100))
    axis(2, at = seq(360,360 + 60,20))
    abline(v=seq(0,max(transect_steps),100),h=seq(round(min(transall[,2])),round(max(transall[,2])),20),
           lty = 2,      # Grid line type
           col = "gray", # Grid line color
           lwd = 0.5)      # Grid line width
    dev.off()
  }
}

if(secNum==217){
  if(i==2){
    png(file=paste(path_figs,'\\',secNum,'_',i,'.png',sep=''), res = 300,width=1800,height=800) 
    par(mar = c(4, 4, 1, 1))
    plot(transect_steps,transall[,2],type='l',lwd=2,xlab ='distance [m]',ylab = 'Elevation [m a.s.l.]',cex = 2,xaxs ='i',xaxt='n',yaxt='n',xlim=c(300,1400),ylim=c(330,500))
    points(transect_steps,transmarg, type='l', col='red',lwd=3)
    axis(1, at = seq(300,1400,100))
    axis(2, at = seq(330,510,20))
    abline(v=seq(0,max(transect_steps),100),h=seq(round(min(transall[,2])),round(max(transall[,2])),20),
           lty = 2,      # Grid line type
           col = "gray", # Grid line color
           lwd = 0.5)      # Grid line width
    dev.off()
  }
}
}
}