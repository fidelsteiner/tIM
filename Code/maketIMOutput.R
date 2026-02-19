################################################################################
# Extract Ice Margin Morphology for separate segments for individual analysis
#
# maketIMOutput.R
#
# ReadMe: Code used to process outputs from the analysis of data from the 
# terrestrial ice margin processing for inclusion in the publication.
# Make sure that folder structures for Sectors stay as in the repository to make sure
# the code works, or alternatively adapt path strings in functions.
# 
# https://github.com/fidelsteiner/tIM
#
# Input:
# 
# Requires all sector outputs from the study stored under \SubbasinGeospatialData\Output
#
# Created:          2025/04/15
# Latest Revision:  2026/02/12
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

# Replace below with respective paths for data
output_basepath <- 'C:\\Work\\Research\\tGISM'            # path for code/final output  
sector_basepath <- 'E:\\tGISM\\Sectors\\'                 # path where all Sector data are stored in folders with name 'XSector'

# numbers of the individual subbabsins used to store sector data
regionSW <- c(137,136,135,134,133,132,131,130,71,70,69,40,22,21,20,19,18,17,16,15,14,13,12,11)
regionSE <- c(63,64,65,67,68,73,80,83,88,92,93,94,106,107,108,110,111,113,114,115,117,118,119,120,121,125,126,128,186,187,188,189,190,191,192,193,199,200,201,202,209,210,211,212,213,231,233,237,243)
regionCE <- c(60,61,62,74,75,76,95,96,100,102,122,203,204,248,249,250,251,252,255)
regionCW <- c(3,4,5,6,7,8,9,10,72,81,82,85,218,219,220,221,222)
regionNE <- c(58,59,87,129,138,139,140,141,142,144,145,217,244,245,246,247)
regionNO <- c(39,50,51,52,53,54,55,56,57,77,84,109,146,176,215,216,223,224)
regionNW <- c(1,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,41,42,43,44,45,46,47,48,49,78,86,103,104,105,127,147,148,159,160,161,162,163,164,165,166,167,168,169,170,172,174,175,180,185,195,196,197,207,208,214,225,238,239,240)
regionPGICNO <- c(302,303,304)
regionPGICNE <- c(305,306,307)
regionPGICCE <- c(308,309,310)
regionPGICSE <- c(311,312,313)
regionPGICSW <- c(314,315)

# process all lake margins of GrIS and summarize in regional files

lakeMarProc <- function(regionName,saveName){
  k<-0
  for(i in regionName){
    if(file.exists(paste(sector_basepath,i,'Sector\\',i,'lakemargin.shp',sep = ""))){
      k<-k+1
      lakeMar <- ogrInfo(paste(sector_basepath,i,'Sector\\',i,'lakemargin.shp',sep = ""))
      lakeMar <-   readOGR(dsn=paste(sector_basepath,i,'Sector\\',i,'lakemargin.shp',sep = ""),layer = lakeMar$layer)
      lakeMar <- spTransform(lakeMar, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
      
      if(k==1){mergedProd<-lakeMar}else{
        mergedProd <- rbind(mergedProd,lakeMar)
      }}
  }
  rgdal::writeOGR(mergedProd,paste(output_basepath,'\\',saveName,'.shp',sep = ""),layer = saveName,driver = 'ESRI Shapefile',overwrite_layer = T)
}

lakeMarProc(regionCW,'1RegionalLakMargin')
lakeMarProc(regionNW,'2RegionalLakMargin')
lakeMarProc(regionNO,'3RegionalLakMargin')
lakeMarProc(regionNE,'4RegionalLakMargin')
lakeMarProc(regionCE,'5RegionalLakMargin')
lakeMarProc(regionSE,'6RegionalLakMargin')
lakeMarProc(regionSW,'7RegionalLakMargin')

# process all marine margins of GrIS and summarize in regional files

marineMarProc <- function(regionName,saveName){
  k<-0
  for(i in regionName){
    if(file.exists(paste(sector_basepath,i,'Sector\\',i,'marinemargin.shp',sep = ""))){
      k<-k+1
      marMar <- ogrInfo(paste(sector_basepath,i,'Sector\\',i,'marinemargin.shp',sep = ""))
      marMar <-   readOGR(dsn=paste(sector_basepath,i,'Sector\\',i,'marinemargin.shp',sep = ""),layer = marMar$layer)
      marMar <- spTransform(lakeMar, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
      
      if(k==1){mergedProd<-marMar}else{
        mergedProd <- rbind(mergedProd,marMar)
      }}
  }
  rgdal::writeOGR(mergedProd,paste(output_basepath,'\\',saveName,'.shp',sep = ""),layer = saveName,driver = 'ESRI Shapefile',overwrite_layer = T)
}

marineMarProc(regionCW,'1RegionalMarMargin')
marineMarProc(regionNW,'2RegionalMarMargin')
marineMarProc(regionNO,'3RegionalMarMargin')
marineMarProc(regionNE,'4RegionalMarMargin')
marineMarProc(regionCE,'5RegionalMarMargin')
marineMarProc(regionSE,'6RegionalMarMargin')
marineMarProc(regionSW,'7RegionalMarMargin')

### Analysis of actual calving fronts vs marine margin identified in the study
# See section on uncertainties in the manuscript

# IDs of sectors to be investigated (with significant marine section)
region_marine <- c(70,69,22,14,63,80,88,106,107,108,110,121,126,189,190,191,202,213,237,243,
                   60,61,74,95,96,248,3,4,6,10,82,85,218,220,58,59,139,244,245,246,39,50,53,56,77,84,215,216,223,
                   23,25,26,27,28,29,30,31,32,33,34,35,36,37,38,41,42,43,44,45,47,48,78,86,103,127,147,148,196,197,207,208,214,238,239)

marlengths <- matrix(nrow=length(region_marine),ncol=2)
i<-0
for(k in region_marine){
  i<-i+1
  marMar <- ogrInfo(paste(output_basepath,'\\Sectors\\MarineSectors_Sensitivity\\',k,'marinemargin.shp',sep = ""))
  marMar <-   readOGR(dsn=paste(output_basepath,'\\Sectors\\MarineSectors_Sensitivity\\',k,'marinemargin.shp',sep = ""),layer = marMar$layer)
  marMar <- spTransform(marMar, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  marlengths[i,1] <- gLength(marMar)/1000 # Actual lengths
  
  marMar_calv <- ogrInfo(paste(output_basepath,'\\Sectors\\MarineSectors_Sensitivity\\',k,'marinemargin_calving.shp',sep = ""))
  marMar_calv <-   readOGR(dsn=paste(output_basepath,'\\Sectors\\MarineSectors_Sensitivity\\',k,'marinemargin_calving.shp',sep = ""),layer = marMar_calv$layer)
  marMar_calv <- spTransform(marMar_calv, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  marlengths[i,2] <- gLength(marMar_calv)/1000 # just calving lengths
}

### Analysis of connectivity of individual morphologies
# Note that in the way this is coded, disaggregated polylines need to be 
# dissolved in QGis/ArcGis after the first loop to enable the second loop to execute
# See Discussion in the manuscript

# GrIS
lengthMatrx <- matrix(NA,ncol=4*7,nrow=10000)

# GrIS Loop 1
for(sb in 1:7){
cliffmargin <- ogrInfo(paste(output_basepath,'\\Sectors\\',sb,'RegionalCliffmargin.shp',sep=''))
cliffmargin <- readOGR(dsn=paste(output_basepath,'\\Sectors\\',sb,'RegionalCliffmargin.shp',sep=''),layer = cliffmargin$layer)

cliffmargin_dis <- disaggregate(cliffmargin)

rgdal::writeOGR(cliffmargin_dis,paste(output_basepath,'\\Sectors\\',sb,'cliffmargin_dis.shp',sep = ""),layer = 'cliffmargin_dis',driver = 'ESRI Shapefile',overwrite_layer = T)

cliffsteepmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalsteepCliffmargin.shp',sep=''))
cliffsteepmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalsteepCliffmargin.shp',sep=''),layer = cliffsteepmargin$layer)

cliffsteepmargin_dis <- disaggregate(cliffsteepmargin)

rgdal::writeOGR(cliffsteepmargin_dis,paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'cliffsteepmargin_dis.shp',sep = ""),layer = 'cliffsteepmargin_dis',driver = 'ESRI Shapefile',overwrite_layer = T)

steepmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalSteepmargin.shp',sep=''))
steepmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalSteepmargin.shp',sep=''),layer = steepmargin$layer)

steepmargin_dis <- disaggregate(steepmargin)

rgdal::writeOGR(steepmargin_dis,paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'steepmargin_dis.shp',sep = ""),layer = 'steepmargin_dis',driver = 'ESRI Shapefile',overwrite_layer = T)

shallowmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalShallowmargin.shp',sep=''))
shallowmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalShallowmargin.shp',sep=''),layer = shallowmargin$layer)

shallowmargin_dis <- disaggregate(shallowmargin)

rgdal::writeOGR(shallowmargin_dis,paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'shallowmargin_dis.shp',sep = ""),layer = 'shallowmargin_dis',driver = 'ESRI Shapefile',overwrite_layer = T)
}

# GrIS Loop 2
kb <- 1
for(sb in 1:7){
  cliffmargin <- ogrInfo(paste(output_basepath,'\\Sectors\\',sb,'cliffmargin_dissolved.shp',sep=''))
  cliffmargin <- readOGR(dsn=paste(output_basepath,'\\Sectors\\',sb,'cliffmargin_dissolved.shp',sep=''),layer = cliffmargin$layer)
  
  cliffMarginL <- gLength(cliffmargin,byid=T)
  cliffMarginL[cliffMarginL<30]<-NA
  
  cliffsteepmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'cliffsteepmargin_dissolved.shp',sep=''))
  cliffsteepmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'cliffsteepmargin_dissolved.shp',sep=''),layer = cliffsteepmargin$layer)
  
  cliffsteepMarginL <- gLength(cliffsteepmargin,byid=T)
  cliffsteepMarginL[cliffsteepMarginL<30]<-NA
  
  steepmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'steepmargin_dissolved.shp',sep=''))
  steepmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'steepmargin_dissolved.shp',sep=''),layer = steepmargin$layer)
  
  steepMarginL <- gLength(steepmargin,byid=T)
  steepMarginL[steepMarginL<30]<-NA
  
  shallowmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'shallowmargin_dissolved.shp',sep=''))
  shallowmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'shallowmargin_dissolved.shp',sep=''),layer = shallowmargin$layer)
  
  shallowMarginL <- gLength(shallowmargin,byid=T)
  shallowMarginL[shallowMarginL<30]<-NA
  
  plot(density(shallowMarginL[-which(is.na(shallowMarginL))]),col='black',ylim=c(0,0.002),xlim=c(0,5000))
  lines(density(steepMarginL[-which(is.na(steepMarginL))]),col='grey')
  lines(density(cliffMarginL[-which(is.na(cliffMarginL))]),col='red')
  lines(density(cliffsteepMarginL[-which(is.na(cliffsteepMarginL))]),col='orange')
  
  median(shallowMarginL,na.rm=T)
  median(steepMarginL,na.rm=T)
  median(cliffMarginL,na.rm=T)
  median(cliffsteepMarginL,na.rm=T)
  
  lengthMatrx[1:length(shallowMarginL),kb] <- shallowMarginL
  lengthMatrx[1:length(steepMarginL),kb+1] <- steepMarginL
  lengthMatrx[1:length(cliffMarginL),kb+2] <- cliffMarginL
  lengthMatrx[1:length(cliffsteepMarginL),kb+3] <- cliffsteepMarginL
  kb <- kb + 4

}

# PGIC
lengthMatrx_PGIC <- matrix(NA,ncol=4*16,nrow=10000)

# PGIC Loop 1
for(sb in 8:23){
  cliffmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalCliffmargin.shp',sep=''))
  cliffmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalCliffmargin.shp',sep=''),layer = cliffmargin$layer)
  
  cliffmargin_dis <- disaggregate(cliffmargin)
  
  rgdal::writeOGR(cliffmargin_dis,paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'cliffmargin_dis.shp',sep = ""),layer = 'cliffmargin_dis',driver = 'ESRI Shapefile',overwrite_layer = T)
  
  
  cliffsteepmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalsteepCliffmargin.shp',sep=''))
  cliffsteepmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalsteepCliffmargin.shp',sep=''),layer = cliffsteepmargin$layer)
  
  cliffsteepmargin_dis <- disaggregate(cliffsteepmargin)
  
  rgdal::writeOGR(cliffsteepmargin_dis,paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'cliffsteepmargin_dis.shp',sep = ""),layer = 'cliffsteepmargin_dis',driver = 'ESRI Shapefile',overwrite_layer = T)
  
  
  steepmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalSteepmargin.shp',sep=''))
  steepmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalSteepmargin.shp',sep=''),layer = steepmargin$layer)
  
  steepmargin_dis <- disaggregate(steepmargin)
  
  rgdal::writeOGR(steepmargin_dis,paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'steepmargin_dis.shp',sep = ""),layer = 'steepmargin_dis',driver = 'ESRI Shapefile',overwrite_layer = T)
 
  shallowmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalShallowmargin.shp',sep=''))
  shallowmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'RegionalShallowmargin.shp',sep=''),layer = shallowmargin$layer)
  
  shallowmargin_dis <- disaggregate(shallowmargin)
  
  rgdal::writeOGR(shallowmargin_dis,paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'shallowmargin_dis.shp',sep = ""),layer = 'shallowmargin_dis',driver = 'ESRI Shapefile',overwrite_layer = T)
  
  
  
}

# PGIC Loop 2
kb <- 1
for(sb in 8:23){
  
  cliffmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'cliffmargin_dissolved.shp',sep=''))
  cliffmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'cliffmargin_dissolved.shp',sep=''),layer = cliffmargin$layer)
  
  cliffMarginL <- gLength(cliffmargin,byid=T)
  cliffMarginL[cliffMarginL<30]<-NA
  
  cliffsteepmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'cliffsteepmargin_dissolved.shp',sep=''))
  cliffsteepmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'cliffsteepmargin_dissolved.shp',sep=''),layer = cliffsteepmargin$layer)
  
  cliffsteepMarginL <- gLength(cliffsteepmargin,byid=T)
  cliffsteepMarginL[cliffsteepMarginL<30]<-NA
  
  steepmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'steepmargin_dissolved.shp',sep=''))
  steepmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'steepmargin_dissolved.shp',sep=''),layer = steepmargin$layer)
  
  steepMarginL <- gLength(steepmargin,byid=T)
  steepMarginL[steepMarginL<30]<-NA
  
  shallowmargin <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'shallowmargin_dissolved.shp',sep=''))
  shallowmargin <- readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',sb,'shallowmargin_dissolved.shp',sep=''),layer = shallowmargin$layer)
  
  shallowMarginL <- gLength(shallowmargin,byid=T)
  shallowMarginL[shallowMarginL<30]<-NA
  
  median(shallowMarginL,na.rm=T)
  median(steepMarginL,na.rm=T)
  median(cliffMarginL,na.rm=T)
  median(cliffsteepMarginL,na.rm=T)
  
  lengthMatrx_PGIC[1:length(shallowMarginL),kb] <- shallowMarginL
  lengthMatrx_PGIC[1:length(steepMarginL),kb+1] <- steepMarginL
  lengthMatrx_PGIC[1:length(cliffMarginL),kb+2] <- cliffMarginL
  lengthMatrx_PGIC[1:length(cliffsteepMarginL),kb+3] <- cliffsteepMarginL
  kb <- kb + 4
}

shallowMarginL_GrIS <- c(lengthMatrx[,seq(1,25,4)])
steepMarginL_GrIS <- c(lengthMatrx[,seq(2,26,4)])
cliffMarginL_GrIS <- c(lengthMatrx[,seq(3,27,4)])
cliffsteepMarginL_GrIS <- c(lengthMatrx[,seq(4,28,4)])

shallowMarginL_PGIC <- c(lengthMatrx_PGIC[,seq(1,61,4)])
steepMarginL_PGIC <- c(lengthMatrx_PGIC[,seq(2,62,4)])
cliffMarginL_PGIC <- c(lengthMatrx_PGIC[,seq(3,63,4)])
cliffsteepMarginL_PGIC <- c(lengthMatrx_PGIC[,seq(4,64,4)])

png(file=paste(output_basepath,'\\Figures\\MorphologySegmentLengths.png',sep = ""), res = 300,width=2800,height=1800,bg = "transparent")
par(bty = 'n') 

plot(density(shallowMarginL_PGIC[-which(is.na(shallowMarginL_PGIC))]),col="#319f28ff",ylim=c(0,0.0022),xlim=c(0,2500),lty=2,xaxt = 'n',ylab ='Density [%]',xlab ='margin length [km]',yaxt='n',main='',yaxs="i",cex.lab=1.5)
axis(1, at=seq(0,2500,500), labels=seq(0,2.5,0.5),cex.axis=1.5)
axis(2, at=seq(0,0.0020,0.0005), labels=seq(0,0.2,0.05),cex.axis=1.5)

lines(density(steepMarginL_PGIC[-which(is.na(steepMarginL_PGIC))]),col="#fc8101ff",lty=2)
lines(density(c(cliffMarginL_PGIC[-which(is.na(cliffMarginL_PGIC))],cliffsteepMarginL_PGIC[-which(is.na(cliffsteepMarginL_PGIC))])),col="#ed0c0cff",lty=2)
#lines(density(cliffsteepMarginL_PGIC[-which(is.na(cliffsteepMarginL_PGIC))]),col="#ed4c0cff",lty=2,lwd=0.5)

lines(density(shallowMarginL_GrIS[-which(is.na(shallowMarginL_GrIS))]),lwd=2,col="#319f28ff",ylim=c(0,0.002),xlim=c(0,2500))
lines(density(steepMarginL_GrIS[-which(is.na(steepMarginL_GrIS))]),lwd=2,col="#fc8101ff")
lines(density(c(cliffMarginL_GrIS[-which(is.na(cliffMarginL_GrIS))],cliffsteepMarginL_GrIS[-which(is.na(cliffsteepMarginL_GrIS))])),lwd=2,col="#ed0c0cff")
#lines(density(cliffsteepMarginL_GrIS[-which(is.na(cliffsteepMarginL_GrIS))]),lwd=0.5,col="#ed4c0cff")

colBar_morph <- c("#319f28ff","#fc8101ff","#ed0c0cff","#ed4c0cff")
seqYloc <- seq(0.00200,0.0005,-0.00008)
seqLabsize <- c(1.5,1.5,1.5,1)
points(c(median(lengthMatrx[,1],na.rm=T),median(lengthMatrx[,2],na.rm=T),median(c(lengthMatrx[,3],lengthMatrx[,4]),na.rm=T)),rep( seqYloc[1],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx[,5],na.rm=T),median(lengthMatrx[,6],na.rm=T),median(c(lengthMatrx[,7],lengthMatrx[,8]),na.rm=T)),rep( seqYloc[2],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx[,9],na.rm=T),median(lengthMatrx[,10],na.rm=T),median(c(lengthMatrx[,11],lengthMatrx[,12]),na.rm=T)),rep( seqYloc[3],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx[,13],na.rm=T),median(lengthMatrx[,14],na.rm=T),median(c(lengthMatrx[,15],lengthMatrx[,16]),na.rm=T)),rep( seqYloc[4],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx[,17],na.rm=T),median(lengthMatrx[,18],na.rm=T),median(c(lengthMatrx[,19],lengthMatrx[,20]),na.rm=T)),rep( seqYloc[5],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx[,21],na.rm=T),median(lengthMatrx[,22],na.rm=T),median(c(lengthMatrx[,23],lengthMatrx[,24]),na.rm=T)),rep( seqYloc[6],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx[,25],na.rm=T),median(lengthMatrx[,26],na.rm=T),median(c(lengthMatrx[,27],lengthMatrx[,28]),na.rm=T)),rep( seqYloc[7],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])

points(c(median(lengthMatrx_PGIC[,1],na.rm=T),median(lengthMatrx_PGIC[,2],na.rm=T),median(c(lengthMatrx_PGIC[,3],lengthMatrx_PGIC[,3]),na.rm=T)),rep( seqYloc[10],3),col=colBar_morph[1:3],cex = seqLabsize[1:3])
points(c(median(lengthMatrx_PGIC[,5],na.rm=T),median(lengthMatrx_PGIC[,6],na.rm=T),median(c(lengthMatrx_PGIC[,7],lengthMatrx_PGIC[,8]),na.rm=T)),rep( seqYloc[11],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx_PGIC[,seq(9,17,4)],na.rm=T),median(lengthMatrx_PGIC[,seq(10,18,4)],na.rm=T),median(c(lengthMatrx_PGIC[,seq(11,19,4)],lengthMatrx_PGIC[,seq(12,20,4)]),na.rm=T)),rep( seqYloc[12],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx_PGIC[,seq(21,29,4)],na.rm=T),median(lengthMatrx_PGIC[,seq(22,30,4)],na.rm=T),median(c(lengthMatrx_PGIC[,seq(23,31,4)],lengthMatrx_PGIC[,seq(24,32,4)]),na.rm=T)),rep( seqYloc[13],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx_PGIC[,seq(33,41,4)],na.rm=T),median(lengthMatrx_PGIC[,seq(34,42,4)],na.rm=T),median(c(lengthMatrx_PGIC[,seq(35,43,4)],lengthMatrx_PGIC[,seq(36,44,4)]),na.rm=T)),rep( seqYloc[14],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx_PGIC[,seq(45,53,4)],na.rm=T),median(lengthMatrx_PGIC[,seq(46,54,4)],na.rm=T),median(c(lengthMatrx_PGIC[,seq(47,55,4)],lengthMatrx_PGIC[,seq(48,56,4)]),na.rm=T)),rep( seqYloc[15],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
points(c(median(lengthMatrx_PGIC[,seq(57,61,4)],na.rm=T),median(lengthMatrx_PGIC[,seq(58,62,4)],na.rm=T),median(c(lengthMatrx_PGIC[,seq(59,63,4)],lengthMatrx_PGIC[,seq(60,64,4)]),na.rm=T)),rep( seqYloc[16],3),col=colBar_morph[1:3],cex=seqLabsize[1:3])
dev.off() 