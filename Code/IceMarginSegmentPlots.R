################################################################################
# Plot margin length statistics
#
# IceMarginSegmentPlots.R
#
# ReadMe: 
#
# Uses sector margins across all subbasins in Greenland to produce overview plot parts and saves the combined segments across
# 7 subregions forthe GrIS and 15 subregions for PGIC
#
# https://github.com/fidelsteiner/tIM
#
# Input:
# Requires results from the margin analysis (see IceMarginExtractor.R)
#
# Created:          2021/10/15
# Latest Revision:  2026/02/16
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
library(truncnorm)

##################
# Paths
##################

output_basepath <- 'C:\\Work\\Research\\tGISM'            # path for code/output   
folderSectors <- paste('E:\\tGISM\\Sectors',sep = "")     # path where all spatial sector data are stored
marginStats <- 'C:\\Work\\Research\\tGISM\\Sectors\\marginstatistics.csv' # file for all margin statistics

##################
# Region numbering
##################
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

PGICCW <- 300
PGICNW <- 301
PGICNOW <- c(302)
PGICNOC <- c(303)
PGICNOE <- c(304)
PGICNEN <- c(305)
PGICNEC <- c(306)
PGICNES <- c(307)
PGICCEN <- c(308)
PGICCEC <- c(309)
PGICCES <- c(310)
PGICSEN <- c(311)
PGICSEC <- c(312)
PGICSES <- c(313)
PGICSWS <- c(314)
PGICSWN <- c(315)

##################
# Read general margin statistics
##################

margStats <- read.csv(marginStats)


colorMargins <- c('#ffa8fcff','#edf0ffff','#a8fff9ff','#abffa7ff','#fffda6ff','#ffbbadff','#cfa8bbff',
                  '#ffa8fcff','#edf0ffff','#a8fff9ff','#a8fff9ff','#a8fff9ff','#abffa7ff','#abffa7ff','#abffa7ff',
                  '#fffda6ff','#fffda6ff','#fffda6ff','#ffbbadff','#ffbbadff','#ffbbadff','#cfa8bbff','#cfa8bbff',
                  '#abffa7ff','#fffda6ff','#ffbbadff','#cfa8bbff')
nameMargins <- c('CW','NW','NO','NE', 'CE', 'SE', 'SW', 
                 'pCW','pNW','pNOW','pNOC','pNOE',
                 'pNEN','pNEC','pNES','pCEN','pCEC','pCES','pSEN','pSEC',
                 'pSES','pSWS','pSWN',
                 'pNO','pNE','pCE','pSE','pSW')

##################
# Remove erroneous margin segments based on manually identified locations
##################

offmarginVec <- vector()

marginErr <- function(regionDat,offmarginVec){
  for(i in 1:length(regionDat)){
    offmargin <- read.csv(paste(folderSectors,'\\',regionDat[i],'Sector\\',regionDat[i],'_offmargin.csv',sep=''))
    if(length(offmargin$ID)>0){
      grid <- ogrInfo(paste(folderSectors,'\\',regionDat[i],'Sector\\',regionDat[i],'Grid.shp',sep=''))
      grid <- readOGR(dsn=paste(folderSectors,'\\',regionDat[i],'Sector\\',sep=''),layer = grid$layer)
      
      wrongGrids <- grid[match(offmargin$ID,grid$id),]
      if(is.null(grid$type)){grid$type <- grid$id*NA}
      grid$type[match(offmargin$ID,grid$id)] <- 6
      terrMargin <- ogrInfo(paste(folderSectors,'\\',regionDat[i],'Sector\\',regionDat[i],'terrestrialmargin.shp',sep = ""))
      terrMargin <- readOGR(dsn=paste(folderSectors,'\\',regionDat[i],'Sector\\',sep = ""),layer = terrMargin$layer)
      
      misplacedmargin <- intersect(terrMargin,wrongGrids)
      offmarginVec[which(margStats$No==regionDat[i])] <- gLength(misplacedmargin)
      
      dfmargin<-SpatialLinesDataFrame(misplacedmargin, data.frame(id=1:length(misplacedmargin)))
      #writeOGR(dfmargin, dsn=paste(folderSectors,'\\',regionDat[i],'Sector\\',regionDat[i],'terrestrialmargin_misplaced.shp',sep = "") ,layer="dfmargin",driver="ESRI Shapefile", overwrite=T)
      writeOGR(grid, dsn=paste(folderSectors,'\\',regionDat[i],'Sector\\',regionDat[i],'Grid.shp',sep = "") ,layer="NA",driver="ESRI Shapefile", overwrite=T)
      
    }
    if(length(offmargin$ID)==0){
      offmarginVec[which(margStats$No==regionDat[i])] <- 0}
  }
  return(offmarginVec) 
}

offmarginVec <- marginErr(regionCW,offmarginVec)
offmarginVec <- marginErr(regionNW,offmarginVec)
offmarginVec <- marginErr(regionNO,offmarginVec)
offmarginVec <- marginErr(regionNE,offmarginVec)
offmarginVec <- marginErr(regionCE,offmarginVec)
offmarginVec <- marginErr(regionSE,offmarginVec)
offmarginVec <- marginErr(regionSW,offmarginVec)
offmarginVec <- marginErr(regionPGICNO,offmarginVec)
offmarginVec <- marginErr(regionPGICNE,offmarginVec)
offmarginVec <- marginErr(regionPGICCE,offmarginVec)
offmarginVec <- marginErr(PGICCW,offmarginVec)
offmarginVec <- marginErr(regionPGICSW,offmarginVec)
offmarginVec <- marginErr(PGICNW,offmarginVec)
offmarginVec <- marginErr(regionPGICSE,offmarginVec)

margStats$terrestrial_corrected <- margStats$terrestrial - offmarginVec
margStats$terrestrial_off <- offmarginVec

write.csv(margStats,file='C:\\Work\\Research\\tGISM\\Sectors\\marginstatistics.csv',row.names=FALSE)

##################
# Plot fraction plots for each subbasin
##################

# Dummy Plot

dataML <- data.frame(
  category = c('terrestrial_corrected','terrestrial_off','marine','lake'),
  count=c(as.vector(c(0.25,0.05,1/3,1/3))) 
)
dataML$fraction <- dataML$count / sum(dataML$count)
dataML$ymax <- cumsum(dataML$fraction)
dataML$ymin <- c(0, head(dataML$ymax, n=-1))
dataML$labelPosition <- (dataML$ymax + dataML$ymin) / 2
dataML$label <- paste0(round(dataML$fraction*100, digits = 1),"%")
pGM <- ggplot(dataML, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(color ='red',lwd=0) +
  scale_fill_manual(values=c(  "#008000","#0000FF","#554400","#BC8F8F")) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        panel.grid = element_blank(),
        title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.margin = unit(0,"null"),
        panel.border = element_blank()) +
  geom_text(aes(label = 'Region'), x=2,y=500,size=30,color ='black')


png(file=paste(output_basepath,'\\Figures\\Donuts\\dummy_marginlengths.png',sep = ""), res = 300,width=3600,height=3600,bg = "transparent")
par(bty = 'n') 
print(pGM,axes=FALSE)
dev.off() 

MarginDonutsGreenland <- function(regionDat, countMargin){
  
  marginDat <- margStats[match(regionDat,margStats[,1]),]
  marginDat <- marginDat[,c(1,2,12,13,4,5,6,7,8,9,10,11,3)]
  
  dataML <- data.frame(
    category = c('terrestrial_corrected','terrestrial_off','marine','lake'),
    count=c(as.vector(colSums(marginDat[,3:6]))) 
  )
  print(sqrt(sum(dataML$count)/1000000)*10)
  dataML$fraction <- dataML$count / sum(dataML$count)
  dataML$ymax <- cumsum(dataML$fraction)
  dataML$ymin <- c(0, head(dataML$ymax, n=-1))
  dataML$labelPosition <- (dataML$ymax + dataML$ymin) / 2
  dataML$label <- paste0(round(dataML$fraction*100, digits = 1),"%")
  print(dataML$label)
  print(sum(dataML$count)/1000)
  pGM <- ggplot(dataML, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect(color =colorMargins[countMargin],lwd=0) +
    #geom_text( x=c(4.3,4.3,4.3), aes(y=labelPosition, label=label), size=30) +
    scale_fill_manual(values=c(  "#008000","#0000FF","#554400","#BC8F8F")) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    theme_void() +
    theme(legend.position = "none",
          panel.background = element_rect(fill='transparent', color=NA), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA),
          #panel.background = element_blank(),
          panel.grid = element_blank(),
          title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.margin = unit(0,"null"),
          panel.border = element_blank())
        
    #geom_text(aes(label = nameMargins[countMargin]), x=2,y=500,size=30,color ='black')
  
  
  
  png(file=paste(output_basepath,'\\Figures\\Donuts\\',nameMargins[countMargin],'_marginlengths.png',sep = ""), res = 300,width=3600,height=3600,bg = "transparent")
  par(bty = 'n') 
  print(pGM,axes=FALSE)
  dev.off() 
}
MarginDonutsGreenland(regionCW, 1)
MarginDonutsGreenland(regionNW, 2)
MarginDonutsGreenland(regionNO, 3)
MarginDonutsGreenland(regionNE, 4)
MarginDonutsGreenland(regionCE, 5)
MarginDonutsGreenland(regionSE, 6)
MarginDonutsGreenland(regionSW, 7)
MarginDonutsGreenland(PGICCW, 8)
MarginDonutsGreenland(PGICNW, 9)
MarginDonutsGreenland(PGICNOW, 10)
MarginDonutsGreenland(PGICNOC, 11)
MarginDonutsGreenland(PGICNOE, 12)
MarginDonutsGreenland(PGICNEN, 13)
MarginDonutsGreenland(PGICNEC, 14)
MarginDonutsGreenland(PGICNES, 15)
MarginDonutsGreenland(PGICCEN, 16)
MarginDonutsGreenland(PGICCEC, 17)
MarginDonutsGreenland(PGICCES, 18)
MarginDonutsGreenland(PGICSEN, 19)
MarginDonutsGreenland(PGICSEC, 20)
MarginDonutsGreenland(PGICSES, 21)
MarginDonutsGreenland(PGICSWS, 22)
MarginDonutsGreenland(PGICSWN, 23)

##################
# Combine margin fractions from all regions and store
##################

# GrIS
for(i in 1:7){
  terrestrialMar <- ogrInfo(paste(output_basepath,'\\Sectors\\',i,'RegionalTerrMargin.shp',sep = ""))
  terrestrialMar <-   readOGR(dsn=paste(output_basepath,'\\Sectors\\',i,'RegionalTerrMargin.shp',sep = ""),layer = terrestrialMar$layer)
  if(i==1){regionalTerrMar <- terrestrialMar}
  if(i>1){regionalTerrMar <- rbind(terrestrialMar,regionalTerrMar)}
  lakeMar <- ogrInfo(paste(output_basepath,'\\Sectors\\',i,'RegionalLakMargin.shp',sep = ""))
  lakeMar <-   readOGR(dsn=paste(output_basepath,'\\Sectors\\',i,'RegionalLakMargin.shp',sep = ""),layer = lakeMar$layer)
  if(i==1){regionalLakeMar <- lakeMar}
  if(i>1){regionalLakeMar <- rbind(lakeMar,regionalLakeMar)}
  marineMar <- ogrInfo(paste(output_basepath,'\\Sectors\\',i,'RegionalMarMargin.shp',sep = ""))
  marineMar <-   readOGR(dsn=paste(output_basepath,'\\Sectors\\',i,'RegionalMarMargin.shp',sep = ""),layer = marineMar$layer)
  if(i==1){regionalMarMar <- marineMar}
  if(i>1){regionalMarMar <- rbind(marineMar,regionalMarMar)}

  errMar <- ogrInfo(paste(output_basepath,'\\Sectors\\',i,'RegionalErrorMargin.shp',sep = ""))
  errMar <-   readOGR(dsn=paste(output_basepath,'\\Sectors\\',i,'RegionalErrorMargin.shp',sep = ""),layer = errMar$layer)
  if(i==1){regionalErrMar <- errMar}
  if(i>1){regionalErrMar <- rbind(errMar,regionalErrMar)}

}
rgdal::writeOGR(regionalTerrMar,paste(output_basepath,'\\Sectors\\RegionalTerrMargin_GrIS.shp',sep = ""),layer = 'regionalTerrMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalMarMar,paste(output_basepath,'\\Sectors\\RegionalMarMargin_GrIS.shp',sep = ""),layer = 'regionalMarMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalLakeMar,paste(output_basepath,'\\Sectors\\RegionalLakMargin_GrIS.shp',sep = ""),layer = 'regionalLakMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalErrMar,paste(output_basepath,'\\Sectors\\RegionalErrorMargin_GrIS.shp',sep = ""),layer = 'regionalErrMargin',driver = 'ESRI Shapefile',overwrite_layer = T)

# PGIC
for(i in 8:23){
  terrestrialMar <- ogrInfo(paste(output_basepath,'\\Sectors\\',i,'RegionalTerrMargin.shp',sep = ""))
  terrestrialMar <-   readOGR(dsn=paste(output_basepath,'\\Sectors\\',i,'RegionalTerrMargin.shp',sep = ""),layer = terrestrialMar$layer)
  if(i==8){regionalTerrMar <- terrestrialMar}
  if(i>8){regionalTerrMar <- rbind(terrestrialMar,regionalTerrMar)}
  lakeMar <- ogrInfo(paste(output_basepath,'\\Sectors\\',i,'RegionalLakMargin.shp',sep = ""))
  lakeMar <-   readOGR(dsn=paste(output_basepath,'\\Sectors\\',i,'RegionalLakMargin.shp',sep = ""),layer = lakeMar$layer)
  if(i==8){regionalLakeMar <- lakeMar}
  if(i>8){regionalLakeMar <- rbind(lakeMar,regionalLakeMar)}
  marineMar <- ogrInfo(paste(output_basepath,'\\Sectors\\',i,'RegionalMarMargin.shp',sep = ""))
  marineMar <-   readOGR(dsn=paste(output_basepath,'\\Sectors\\',i,'RegionalMarMargin.shp',sep = ""),layer = marineMar$layer)
  if(i==8){regionalMarMar <- marineMar}
  if(i>8){regionalMarMar <- rbind(marineMar,regionalMarMar)}
  
  errMar <- ogrInfo(paste(output_basepath,'\\Sectors\\',i,'RegionalErrorMargin.shp',sep = ""))
  errMar <-   readOGR(dsn=paste(output_basepath,'\\Sectors\\',i,'RegionalErrorMargin.shp',sep = ""),layer = errMar$layer)
  if(i==8){regionalErrMar <- errMar}
  if(i>8){regionalErrMar <- rbind(errMar,regionalErrMar)}
}
rgdal::writeOGR(regionalTerrMar,paste(output_basepath,'\\Sectors\\RegionalTerrMargin_PGIC.shp',sep = ""),layer = 'regionalTerrMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalMarMar,paste(output_basepath,'\\Sectors\\RegionalMarMargin_PGIC.shp',sep = ""),layer = 'regionalMarMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalLakeMar,paste(output_basepath,'\\Sectors\\RegionalLakMargin_PGIC.shp',sep = ""),layer = 'regionalLakMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalErrMar,paste(output_basepath,'\\Sectors\\RegionalErrorMargin_PGIC.shp',sep = ""),layer = 'regionalErrMargin',driver = 'ESRI Shapefile',overwrite_layer = T)