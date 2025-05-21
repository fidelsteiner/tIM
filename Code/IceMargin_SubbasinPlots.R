################################################################################
# Produce figures 
#
# IceMargin_SubbasinPlots.R
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
output_basepath2 <- 'E:\\tGISM'

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
folderSectors <- paste('E:\\tGISM\\Sectors',sep = "")
bedM <- raster(paste(data_Bedmachine,'\\BedmachineRaster_proj.tif',sep=''))


projec <-'+proj=utm +datum=WGS84'

projec_utm <-'+proj=utm +zone=19 +north +datum=WGS84'
projec_27N <-'+proj=utm +zone=27N +datum=WGS84'

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


# Concept Plot Slope distribution
library(truncnorm)
png(file=paste(output_basepath,'\\Figures\\WorkFlow\\dummy_slopedistributions.png',sep = ""), res = 300,width=1200,height=1200)
par(mar = c(4.5, 4.5, 0.5, 0.5))
plot(ecdf(rtruncnorm(n=1000, a=5, b=85, mean=39.4, sd=12.09)),lwd=3, col='#ffbbadff',xlim=c(0,90),ylab = 'cumulative frequency',xlab='slope [°]',main='',cex.lab=1.5,cex.axis=1.5)
plot(ecdf(rtruncnorm(n=1000, a=5, b=85, mean=29.4, sd=22.09)),lwd=3, col='#feb0aefe',add=T)
plot(ecdf(rtruncnorm(n=1000, a=5, b=85, mean=19.4, sd=15.09)),lwd=3, col='#efbbadff',add=T)
plot(ecdf(rtruncnorm(n=1000, a=22, b=78, mean=19.4, sd=15.09)),lwd=3, col='#efbbadef',add=T)
plot(ecdf(rtruncnorm(n=1000, a=5, b=85, mean=29.4, sd=22.09)),lwd=3, col='#efbbadef',add=T)
plot(ecdf(rtruncnorm(n=1000, a=5, b=95, mean=59.4, sd=22.09)),lwd=3, col='#efbbadef',add=T)
plot(ecdf(rtruncnorm(n=3000, a=10, b=85, mean=29.4, sd=2.09)),lwd=3, col='#efbbadef',add=T)
dev.off()

# Read general margin statistics
marginStats <- 'C:\\Work\\Research\\tGISM\\Sectors\\marginstatistics.csv'

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

# Remove erroneous margin segments

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
# Plot fraction plots for each subbasin

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
  geom_rect(color ='red',lwd=8) +
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
    geom_rect(color =colorMargins[countMargin],lwd=8) +
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



# Plot slope statistics

MarginSlopesGreenland <- function(regionDat, countMargin){
  
med <- vector()
hig <- vector()
low <- vector()
cnt <- vector()

png(file=paste(output_basepath,'\\Figures\\SlopeGraphs\\',nameMargins[countMargin],'_marginslopes.png',sep = ""), res = 300,width=3600,height=3600)
par(bty = 'n') 
plot(c(0,0), c(0,0), type = "n",xlim = c(0 , 1), ylim = c( 0, 1),xaxs='i', xlab ='',ylab='',xaxt = "n",yaxt = "n") 
axis(side = 1, lwd = 4,labels =F)
axis(side = 2, lwd = 4,labels =F)
axis(side = 3, lwd = 4,labels =F)
axis(side = 4, lwd = 4,labels =F)
abline(a=0,b=1,col='black',lwd=8, lty=2,xaxt = "n",yaxt = "n")
abline(a=0,b=33/45,col='black',lwd=8, lty=2,xaxt = "n",yaxt = "n")

if(countMargin<10000){
for(i in 1:length(regionDat)){
    slopeStats <- read.csv(paste(folderSectors,'\\',regionDat[i],'Sector\\',regionDat[i],'_slopePercentiles.csv',sep=''))

  med[i] <- slopeStats[dim(slopeStats)[1],4]
  hig[i] <- slopeStats[dim(slopeStats)[1],6]
  low[i] <- slopeStats[dim(slopeStats)[1],2]
  cnt[i] <- slopeStats[dim(slopeStats)[1],8]
  abline(a=0,b=med[i]/45,col='red',lwd=1,xaxt = "n",yaxt = "n")
  abline(a=0,b=hig[i]/45,col='black',lwd=1,xaxt = "n",yaxt = "n")
  abline(a=0,b=low[i]/45,col='blue',lwd=1,xaxt = "n",yaxt = "n")
}
abline(a=0,b=sum((cnt/sum(cnt))*med)/45,col='red',lwd=8,lty=2,xaxt = "n",yaxt = "n")
text(0.3,0.9,paste0(round(sum((cnt/sum(cnt))*med),1),"°"),cex=10,col='red')
dev.off() 
}

png(file=paste(output_basepath,'\\Figures\\SlopeGraphs\\dummy_marginslopes.png',sep = ""), res = 300,width=3600,height=3600)
par(bty = 'n') 
plot(c(0,0), c(0,0), type = "n",xlim = c(0 , 1), ylim = c( 0, 1),xaxs='i', xlab ='',ylab='',xaxt = "n",yaxt = "n") 
axis(side = 1, lwd = 4,labels =F)
axis(side = 2, lwd = 4,labels =F)
axis(side = 3, lwd = 4,labels =F)
axis(side = 4, lwd = 4,labels =F)
abline(a=0,b=1,col='black',lwd=8, lty=2,xaxt = "n",yaxt = "n")
abline(a=0,b=33/45,col='black',lwd=8, lty=2,xaxt = "n",yaxt = "n")


  abline(a=0,b=23/45,col='red',lwd=1,xaxt = "n",yaxt = "n")
  abline(a=0,b=65/45,col='black',lwd=1,xaxt = "n",yaxt = "n")
  abline(a=0,b=8/45,col='blue',lwd=1,xaxt = "n",yaxt = "n")

abline(a=0,b=25/45,col='red',lwd=8,lty=2,xaxt = "n",yaxt = "n")
text(0.3,0.9,paste0('median',"°"),cex=6,col='black')
dev.off() 
}

MarginSlopesGreenland(regionCW, 1)
MarginSlopesGreenland(regionNW, 2)
MarginSlopesGreenland(regionNO, 3)
MarginSlopesGreenland(regionNE, 4)
MarginSlopesGreenland(regionCE, 5)
MarginSlopesGreenland(regionSE, 6)
MarginSlopesGreenland(regionSW, 7)
MarginSlopesGreenland(PGICCW, 8)
MarginSlopesGreenland(PGICNW, 9)
MarginSlopesGreenland(PGICNOW, 10)
MarginSlopesGreenland(PGICNOC, 11)
MarginSlopesGreenland(PGICNOE, 12)
MarginSlopesGreenland(PGICNEN, 13)
MarginSlopesGreenland(PGICNEC, 14)
MarginSlopesGreenland(PGICNES, 15)
MarginSlopesGreenland(PGICCEN, 16)
MarginSlopesGreenland(PGICCEC, 17)
MarginSlopesGreenland(PGICCES, 18)
MarginSlopesGreenland(PGICSEN, 19)
MarginSlopesGreenland(PGICSEC, 20)
MarginSlopesGreenland(PGICSES, 21)
MarginSlopesGreenland(PGICSWS, 22)
MarginSlopesGreenland(PGICSWN, 23)
MarginSlopesGreenland(regionPGICNO, 24)
MarginSlopesGreenland(regionPGICNE, 25)
MarginSlopesGreenland(regionPGICCE, 26)
MarginSlopesGreenland(regionPGICSE, 27)
MarginSlopesGreenland(regionPGICSW, 28)


# determine distributions of known cliffs by creating CDFs for (a) cliffs, (b) steep ramps and (c) shallow ramps. Those CDFs will then be compared 
# to all samples by a Kolmogorov Smirnoff Test to establish if the cell is either of the three

# Sample cells
# 301: cliff c(51963, 52256, 51962,52547, 52840, 53692, 51966)
# steep ramp c(5255, 53131, 53423)
# 46:
# shallow ramp: c(3093, 3707, 3708)
# 176: cliff c(7700, 7872, 7563, 7184, 7528, 6839, 7183, 7011, 6667, 6495, 6323, 6151,5979, 5980,5808)
# steep ramp c(7868, 7869, 7870, 7871, 7866, 8341, 8342, 8510, 8511)
# shallow ramp c(9346, 9347, 9518, 9519, 9691, 9692, 9693, 9864, 9865)

SubBasinID_1 <- 176
SubBasinID_2 <- 301
SubBasinID_3 <- 46

SubBasinID_4 <- 4
SubBasinID_5 <- 85
SubBasinID_6 <- 306
# CDF for cliffs

sPerc_1 <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID_1,'Sector\\',SubBasinID_1,'_slopePercentiles.csv',sep = ""))
sPerc_2 <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID_2,'Sector\\',SubBasinID_2,'_slopePercentiles.csv',sep = ""))
sPerc_3 <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID_3,'Sector\\',SubBasinID_3,'_slopePercentiles.csv',sep = ""))
sPerc_4 <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID_4,'Sector\\',SubBasinID_4,'_slopePercentiles.csv',sep = ""))
sPerc_5 <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID_5,'Sector\\',SubBasinID_5,'_slopePercentiles.csv',sep = ""))
sPerc_6 <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID_6,'Sector\\',SubBasinID_6,'_slopePercentiles.csv',sep = ""))


sGrid_1 <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID_1,'Sector\\',SubBasinID_1,'Grid.shp',sep = ""))
sGrid_1 <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID_1,'Sector\\',SubBasinID_1,'Grid.shp',sep = ""),layer = sGrid_1$layer)
sGrid_2 <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID_2,'Sector\\',SubBasinID_2,'Grid.shp',sep = ""))
sGrid_2 <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID_2,'Sector\\',SubBasinID_2,'Grid.shp',sep = ""),layer = sGrid_2$layer)
sGrid_3 <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID_3,'Sector\\',SubBasinID_3,'Grid.shp',sep = ""))
sGrid_3 <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID_3,'Sector\\',SubBasinID_3,'Grid.shp',sep = ""),layer = sGrid_3$layer)
sGrid_4 <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID_4,'Sector\\',SubBasinID_4,'Grid.shp',sep = ""))
sGrid_4 <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID_4,'Sector\\',SubBasinID_4,'Grid.shp',sep = ""),layer = sGrid_4$layer)
sGrid_5 <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID_5,'Sector\\',SubBasinID_5,'Grid.shp',sep = ""))
sGrid_5 <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID_5,'Sector\\',SubBasinID_5,'Grid.shp',sep = ""),layer = sGrid_5$layer)
sGrid_6 <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID_6,'Sector\\',SubBasinID_6,'Grid.shp',sep = ""))
sGrid_6 <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID_6,'Sector\\',SubBasinID_6,'Grid.shp',sep = ""),layer = sGrid_6$layer)


SubBasinID_test <- 301
sPerc_test <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID_test,'Sector\\',SubBasinID_test,'_slopePercentiles.csv',sep = ""))
sGrid_test <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID_test,'Sector\\',SubBasinID_test,'Grid.shp',sep = ""))
sGrid_test <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID_test,'Sector\\',SubBasinID_test,'Grid.shp',sep = ""),layer = sGrid_test$layer)


SubBasinID_test2 <- 301
sPerc_test2 <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID_test2,'Sector\\',SubBasinID_test2,'_slopePercentiles.csv',sep = ""))
sGrid_test2 <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID_test2,'Sector\\',SubBasinID_test2,'Grid.shp',sep = ""))
sGrid_test2 <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID_test2,'Sector\\',SubBasinID_test2,'Grid.shp',sep = ""),layer = sGrid_test2$layer)



CDF_cliff <- vector()
segk_cliff <- match(c(7700, 7872, 7528, 6839, 7183, 7011, 6667, 5979),sPerc_1$grid.ID)
for(i in 1:length(segk_cliff)){
  gridslVal <- unlist(read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID_1,'Sector\\',SubBasinID_1,'_',sPerc_1[segk_cliff[i],1],'_slopeValues.txt',sep = ""),fill = TRUE))
  gridslVal[gridslVal<5]<-NA
  #quantile(unlist(gridslVal), probs = c(0.05, 0.25, 0.5, 0.75, 0.95),na.rm=T)
  CDF_cliff <- c(CDF_cliff,unlist(gridslVal))
}

CDF_cliff2 <- vector()
segk_cliff <- match(c(51963, 52256),sPerc_2$grid.ID)
for(i in 1:length(segk_cliff)){
  gridslVal <- unlist(read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID_2,'Sector\\',SubBasinID_2,'_',sPerc_2[segk_cliff[i],1],'_slopeValues.txt',sep = ""),fill = TRUE))
  gridslVal[gridslVal<5]<-NA
  #quantile(unlist(gridslVal), probs = c(0.05, 0.25, 0.5, 0.75, 0.95),na.rm=T)
  CDF_cliff2 <- c(CDF_cliff2,unlist(gridslVal))
}

CDF_steep<-vector()
segk_steep <- match(c(7868, 7869, 7870, 7871, 7866),sPerc_1$grid.ID)

for(i in 1:length(segk_steep)){
  gridslVal <- unlist(read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID_1,'Sector\\',SubBasinID_1,'_',sPerc_1[segk_steep[i],1],'_slopeValues.txt',sep = ""),fill = TRUE))
  gridslVal[gridslVal<5]<-NA
  CDF_steep<-c(CDF_steep,unlist(gridslVal))
}

CDF_steep2<-vector()
segk_steep <- match(c(52255, 53131, 53423),sPerc_2$grid.ID)

for(i in 1:length(segk_steep)){
  gridslVal <- unlist(read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID_2,'Sector\\',SubBasinID_2,'_',sPerc_2[segk_steep[i],1],'_slopeValues.txt',sep = ""),fill = TRUE))
  gridslVal[gridslVal<5]<-NA
  CDF_steep2<-c(CDF_steep2,unlist(gridslVal))
}

CDF_shallow <- vector()
segk_shallow <- match(c(9346, 9347, 9518, 9519, 9691, 9692, 9693, 9864, 9865),sPerc_1$grid.ID)

for(i in 1:length(segk_shallow)){
  gridslVal <- unlist(read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID_1,'Sector\\',SubBasinID_1,'_',sPerc_1[segk_shallow[i],1],'_slopeValues.txt',sep = ""),fill = TRUE))
  gridslVal[gridslVal<5]<-NA
  CDF_shallow <- c(CDF_shallow,unlist(gridslVal))
}

CDF_shallow2 <- vector()
segk_shallow <- match(c(3053, 3128,3129, 3585, 3659, 3660),sPerc_3$grid.ID)

for(i in 1:length(segk_shallow)){
  gridslVal <- unlist(read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID_3,'Sector\\',SubBasinID_3,'_',sPerc_3[segk_shallow[i],1],'_slopeValues.txt',sep = ""),fill = TRUE))
  gridslVal[gridslVal<5]<-NA
  CDF_shallow2 <- c(CDF_shallow2,unlist(gridslVal))
}

CDF_shallow3 <- vector()
segk_terrain <- match(c(98943,98944, 99188, 99189,99191),sPerc_6$grid.ID)

for(i in 1:length(segk_terrain)){
  gridslVal <- unlist(read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID_6,'Sector\\',SubBasinID_6,'_',sPerc_6[segk_terrain[i],1],'_slopeValues.txt',sep = ""),fill = TRUE))
  gridslVal[gridslVal<5]<-NA
  CDF_shallow3  <- c(CDF_shallow3,unlist(gridslVal))
}

CDF_terrain_mountain <- vector()
segk_terrain <- match(c(178, 179),sPerc_5$grid.ID)

for(i in 1:length(segk_terrain)){
  gridslVal <- unlist(read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID_5,'Sector\\',SubBasinID_5,'_',sPerc_5[segk_terrain[i],1],'_slopeValues.txt',sep = ""),fill = TRUE))
  gridslVal[gridslVal<5]<-NA
  CDF_terrain_mountain  <- c(CDF_terrain_mountain,unlist(gridslVal))
}

#Total n used cells
length(which(!is.na(CDF_cliff)))+length(which(!is.na(CDF_cliff2))) + length(which(!is.na(CDF_terrain_mountain))) +
  length(which(!is.na(CDF_steep)))+length(which(!is.na(CDF_steep2))) +
  length(which(!is.na(CDF_shallow)))+length(which(!is.na(CDF_shallow2))) +length(which(!is.na(CDF_shallow3)))

png(file=paste(path_figs,'\\_CDF_slopes.png',sep=''), res = 300,width=2000,height=1800)  

plot(ecdf(CDF_cliff),xlim=c(0,90),xlab='slope β [°]',ylab='F(β)',main=NULL,lty=1,col='#ed0c0c',lwd=3,cex=2)
plot(ecdf(CDF_cliff2),add=T,col='#ed0c0c',lty=1,lwd=3,cex=2)

plot(ecdf(CDF_steep),add=T,col='#fc8101',lty=1,lwd=3,cex=2)
plot(ecdf(CDF_steep2),add=T,col='#fc8101',lty=1,lwd=3,cex=2)
#plot(ecdf(CDF_steep3),add=T,col='#fc8101',lty=1,lwd=3,cex=2)
#plot(ecdf(CDF_steep[which(CDF_steep>30)]),add=T,col='black',lty=2,lwd=1,cex=2)
plot(ecdf(CDF_shallow),add=T,col='#319f28',lty=1,lwd=3,cex=2)
plot(ecdf(CDF_shallow2),add=T,col='#319f28',lty=1,lwd=3,cex=2)
plot(ecdf(CDF_shallow3),add=T,col='#319f28',lty=1,lwd=3,cex=2)

plot(ecdf(CDF_terrain_mountain),add=T,col='#ed4c0c',lty=1,lwd=3,cex=2)
#plot(ecdf(CDF_terrain_mountain2),add=T,col='green',lty=1,lwd=2,cex=2)
abline(h=0.85,col="grey", lwd=2, lty=2)
abline(h=0.5,col="grey", lwd=2, lty=2)

#plot(ecdf(CDF_terrain_olsen),add=T,col='blue',lty=1,lwd=2,cex=2)

#plot(ecdf(CDF_steep[which(CDF_steep>20)]),add=T,col='black',lty=2,lwd=1,cex=2)
#plot(ecdf(CDF_cliff[which(CDF_cliff>20)]),add=T,col='red',lty=2,lwd=1,cex=2)
#plot(ecdf(CDF_terrain_mountain[which(CDF_terrain_mountain>20)]),add=T,col='green',lty=2,lwd=1,cex=2)

legend('bottomright',legend=c("near-vertical margin", "steep ramp",'shallow ramp'),col=c('#ed0c0c','#fc8101','#319f28'),lty=c(1),lwd=2, bty='n',cex=1.5)
dev.off()

# Slope violin plots per grid cell

# Produce base slope files from BedMachine (only run after new Grid Cells produced)

BaseSlope <- function(regionDat,LOC){
  
  for(subk in 1:length(regionDat)){ 
    
    SubBasinID <- regionDat[subk]
    sPerc <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_slopePercentiles.csv',sep = ""))
    sPerc <- sPerc[-dim(sPerc)[1],]
    
    numGrid <- dim(sPerc)[1]
    
    sGrid <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'Grid.shp',sep = ""))
    sGrid <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'Grid.shp',sep = ""),layer = sGrid$layer)
    
    sGridClip <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'GridClipped.shp',sep = ""))
    sGridClip <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'GridClipped.shp',sep = ""),layer = sGridClip$layer)
    
    bedM_subbasin <- crop(bedM,extent(spTransform(sGrid,CRS("+proj=longlat +datum=WGS84 +no_defs"))))
    
    slopeBase <- rep(NA, numGrid)
    elevBase <- rep(NA, numGrid)
    slopeGrid2 <- terrain(bedM_subbasin, opt="slope", unit="degrees", neighbors=8)
    
    for(segk in 1:numGrid){
      if(!is.null(intersect(extent(slopeGrid2),extent(spTransform(sGridClip[segk,],CRS("+proj=longlat +datum=WGS84 +no_defs")))))){
        elevGrid <- crop(bedM_subbasin,spTransform(sGridClip[segk,],CRS("+proj=longlat +datum=WGS84 +no_defs")))
        slopeGrid <- crop(slopeGrid2,spTransform(sGridClip[segk,],CRS("+proj=longlat +datum=WGS84 +no_defs")))
      }else{next}
      
      slopeBase[segk] <- cellStats(slopeGrid, 'mean')
      elevBase[segk] <- cellStats(elevGrid, 'mean')
      
      if(segk %% 100==0) {
        cat(paste0("basin: ", floor(segk/numGrid*100), "%"))}
    }
    write.csv(cbind(sPerc[,1],slopeBase),paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'slopeBase.csv',sep = ""),row.names=F)
    write.csv(cbind(sPerc[,1],elevBase),paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'elevBase.csv',sep = ""),row.names=F)
    
    
    
    if(subk %% 1==0) {
      cat(paste0("region: ", floor(subk/length(regionDat)*100), "%"))}
  }
}
BaseSlope(regionNO,'NO')
BaseSlope(regionNW,'NW')
BaseSlope(regionCW,'CW')
BaseSlope(regionNE,'NE')
BaseSlope(regionCE,'CE')
BaseSlope(regionSE,'SE')
BaseSlope(regionSW,'SW')
BaseSlope(PGICCW,'pCW')
BaseSlope(PGICNW,'pNW')
BaseSlope(regionPGICNO,'pNO')
BaseSlope(regionPGICNE,'pNE')
BaseSlope(regionPGICCE,'pCE')
BaseSlope(regionPGICSE,'pSE')
BaseSlope(regionPGICSW,'pSW')

library(twosamples)
MarginGridSections <- function(regionDat,LOC){
  steepLength <- vector()
  shallowLength <- vector()
  cliffLength <- vector()
  truecliffLength <- vector()
  
  for(subk in 1:length(regionDat)){ 
    
    SubBasinID <- regionDat[subk]
    sPerc <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_slopePercentiles.csv',sep = ""))
    offMargin <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_offmargin.csv',sep = ""))
    slopeBase <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'slopeBase.csv',sep = ""))
    elevBase <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'elevBase.csv',sep = ""))
    
    #sPerc <- sPerc[-match(offMargin$ID,sPerc$grid.ID),]
    sPerc <- sPerc[-dim(sPerc)[1],]
    
    numGrid <- dim(sPerc)[1]
    
    sGrid <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'Grid.shp',sep = ""))
    sGrid <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'Grid.shp',sep = ""),layer = sGrid$layer)
    
    sGridClip <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'GridClipped.shp',sep = ""))
    sGridClip <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'GridClipped.shp',sep = ""),layer = sGridClip$layer)
  
    shallowRamp <- vector()
    steepRamp <- vector()
    cliffRamp <- vector()
    truecliffRamp<- vector()
    slopeclass <- matrix(NA, nrow = numGrid,ncol = 4)
    
    # Save output for this subbasin
    tabSB <- matrix(NA, ncol = 8, nrow = numGrid+1)
    tabSB[,1] <- c(1:numGrid,'TOTAL')
    
    typeID <- rep(NA, numGrid)
    elevID <- rep(NA, numGrid)
    
    for(segk in 1:numGrid){
      if(!is.na(match(sGrid$id[segk],offMargin$ID))){
        typeID[segk] <- 6
        next}
      if(!is.na(slopeBase$slopeBase[segk])&slopeBase$slopeBase[segk]>30){typeID[segk] <- 5
        next}

      if(!is.na(match(as.numeric(sPerc[segk,1]),offMargin$ID))){typeID[segk]<-6}
      if(is.na(match(as.numeric(sPerc[segk,1]),offMargin$ID))){
      gridslVal <- read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_',sPerc[segk,1],'_slopeValues.txt',sep = ""),fill = TRUE)
      gridslVal[gridslVal<5]<-NA
      
      #if(!is.na(slopeBase$slopeBase[segk])&slopeBase$slopeBase[segk]>10){
      #  gridslVal<-gridslVal-slopeBase$slopeBase[segk]
      #  gridslVal[gridslVal<5]<-NA
      #}
      if(mean(unlist(gridslVal),na.rm=T)!='NaN'&length(which(!is.na(unlist(gridslVal))))>1000){
        
        #plot(ecdf(unlist(gridslVal)),add=T,lty=1,lwd=1,cex=2,col='grey')
        #plot(ecdf(CDF_cliff),add=T,col='red')
        #plot(ecdf(CDF_steep),add=T,col='green')
        
        ksCliff <- dts_stat(sort(unlist(gridslVal)),sort(CDF_cliff))
        ksCliff2 <- dts_stat(sort(unlist(gridslVal)),sort(CDF_cliff2))
        ksMountainCliff <- dts_stat(sort(unlist(gridslVal)),sort(CDF_terrain_mountain))
        ksSteep <- dts_stat(sort(unlist(gridslVal)),sort(CDF_steep))
        ksSteep2 <- dts_stat(sort(unlist(gridslVal)),sort(CDF_steep2))
        ksSteep3 <- dts_stat(sort(unlist(gridslVal)),sort(CDF_steep3))
        ksShallow <- dts_stat(sort(unlist(gridslVal)),sort(CDF_shallow))
        ksShallow2 <- dts_stat(sort(unlist(gridslVal)),sort(CDF_shallow2))
        ksShallow3 <- dts_stat(sort(unlist(gridslVal)),sort(CDF_shallow3))
        
        marginMorph <- which.min(c(ksShallow,ksShallow2,ksShallow3,ksSteep,ksSteep2,ksCliff,ksCliff2,ksMountainCliff))
        
        if(marginMorph<4){typeID[segk]=3}
        if(marginMorph==4|marginMorph==5){typeID[segk]=2}
        if(marginMorph==6|marginMorph==7){typeID[segk]=1}
        if(marginMorph==8){typeID[segk]=4}
          #if(min(ksCliff, ksSteep,ksMountainCliff)>dts_stat(sort(CDF_cliff),sort(CDF_shallow))){typeID[segk] <- 7}
          #typeID[segk] <- unlist(which(c(ksCliff, ksSteep) == min(ksCliff, ksSteep)))
        #}
      }else{
        typeID[segk] <- 6
      }  
        # store grid cells according to their steepness class
        #if(length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<10))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])>0.5){
        #  shallowRamp <- c(shallowRamp,sPerc[segk,1])
        #shallowLength <- c(shallowLength,sPerc[segk,7]*length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<10))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]))
          
        #}
        #if(length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<45&unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>20))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])>0.1){
        #  steepRamp <- c(steepRamp,sPerc[segk,1])
        #steepLength <- c(steepLength,sPerc[segk,7]*length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<45&unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>20))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]))
        #}
        #if(length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>45))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])>0.1){
        #  cliffRamp <- c(cliffRamp,sPerc[segk,1])
        #cliffLength <- c(cliffLength,sPerc[segk,7]*length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>65))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]))
        #}
      #  
        #if(length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>65))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])>0.05){
        #  truecliffRamp <- c(truecliffRamp,sPerc[segk,1])
        #  truecliffLength <- c(truecliffLength,sPerc[segk,7]*length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>65))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]))
        #  
        #}
        
        
        #tabSB[segk,2:6] <- quantile(unlist(gridslVal), probs = c(0.05, 0.25, 0.5, 0.75, 0.95),na.rm=T)
        #tabSB[segk,7] <- sPerc[segk,7]
        
        #tabSB[segk,8] <- sPerc[segk,8]
        
        
        #slopeclass[segk,1] <- length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<10))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])
        #slopeclass[segk,2] <- length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>=10 & unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<20))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])
        #slopeclass[segk,3] <- length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>=20 & unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]<45))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])
        #slopeclass[segk,4] <- length(which(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))]>=45))/length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])
        
        #wdata_violin = data.frame(
        #  location = factor(rep(c(SubBasinID),c(length(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])))),
        #  slope = c(unlist(gridslVal)[which(!is.na(unlist(gridslVal)))])
        #)
        
        #mu <- wdata_violin%>% 
        #  # arrange(order)  %>% 
        #  group_by(location) %>%
        #  summarise(grp.mean = median(slope,na.rm=T))
        
       # tryCatch({
          #png(file=paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_',sPerc[segk,1],'_violin.png',sep=''), res = 300,width=2000,height=900)  
          #par(mar=c(4,2,5,1),cex.lab=1.2,cex.axis=1.5)
          #print(ggplot(wdata_violin, aes(x = slope, y = location)) +
          #        geom_density_ridges(aes(fill = location)) +
          #        scale_fill_manual(values = brewer.pal(n = 5, name = "Dark2")) +
          #        geom_vline(aes(xintercept = grp.mean, color = location),
          #                   data = mu, linetype = "dashed") +
          #        scale_x_continuous(name=expression('Slope [°]'), limits=c(0, 95)) + 
          #        #scale_y_continuous(name='', limits=c(1)) +
          #        theme_bw() +
          #        theme(axis.text.x = element_text(size=14, angle=270),
          #              axis.text.y = element_text(size=14, angle=0)) +
          #        theme(axis.title.y =element_blank(),
          #              axis.text.y=element_blank()) +
          #        theme(legend.position="none"))
          #dev.off()
        #}, error=function(e){})
      }
      if(segk %% 100==0) {
        cat(paste0("basin: ", floor(segk/numGrid*100), "%"))}
    }
    
    sGrid$type <- typeID[1:length(sGrid$type)]
    sGrid$elev <- elevBase$elevBase
    
    rgdal::writeOGR(sGrid,paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'Grid.shp',sep = ""),layer = 'sGrid',driver = 'ESRI Shapefile',overwrite_layer = T)
    
    png(file=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_slopeclasses.png',sep=''), res = 300,width=2000,height=900) 
    barplot(t(slopeclass),col=brewer.pal(n = 4, name = "Reds"),cex.names=1,
            names.arg=sPerc[1:numGrid,1],ylim=c(0,1))
    legend("top", fill = brewer.pal(n = 4, name = "Reds"), legend = c('<10','10-20','20-45','>45'), 
           horiz = TRUE, inset = c(0,-0.5), xpd = TRUE)
    dev.off()
    

    if(subk %% 1==0) {
      cat(paste0("region: ", floor(subk/length(regionDat)*100), "%"))}
  }
 
 # write.csv(steepLength,paste(output_basepath,'\\Sectors\\steepGrids\\',LOC,'_steeplengths.csv',sep = ""))
 # write.csv(cliffLength,paste(output_basepath,'\\Sectors\\steepGrids\\',LOC,'_clifflengths.csv',sep = ""))
 # write.csv(shallowLength,paste(output_basepath,'\\Sectors\\steepGrids\\',LOC,'_shallowlengths.csv',sep = ""))
  #write.csv(truecliffLength,paste(output_basepath,'\\Sectors\\steepGrids\\',LOC,'_trueclifflengths.csv',sep = ""))
  
}



MarginGridSections(regionCW,'CW')

MarginGridSections(regionNE,'NE')

MarginGridSections(regionCE,'CE')

MarginGridSections(regionSE,'SE')

MarginGridSections(regionSW,'SW')

MarginGridSections(PGICCW,'pCW')

MarginGridSections(PGICNW,'pNW')

MarginGridSections(regionPGICNO,'pNO')

MarginGridSections(regionPGICCE,'pCE')

MarginGridSections(regionPGICSE,'pSE')

MarginGridSections(regionPGICSW,'pSW')

MarginGridSections(regionNO,'NO')

MarginGridSections(regionNW,'NW')

MarginGridSections(regionPGICNE,'pNE')

dataML <- data.frame(
  category = c('cliff','steepcliff','steep','shallow','NA'),
  count=c(as.vector(c(0.1,0.05,0.5,0.3,0.1))) 
)
dataML$fraction <- dataML$count / sum(dataML$count)
dataML$ymax <- cumsum(dataML$fraction)
dataML$ymin <- c(0, head(dataML$ymax, n=-1))
dataML$labelPosition <- (dataML$ymax + dataML$ymin) / 2
dataML$label <- paste0(round(dataML$fraction*100, digits = 1),"%")
pGM <- ggplot(dataML, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(color ='black',lwd=4) +
  #geom_text( x=c(4.3,4.3,4.3), aes(y=labelPosition, label=label), size=30) +
  scale_fill_manual(values=c(  "#ed0c0cff","#BC8F8F","#319f28ff","#fc8101ff","#ed4c0cff")) +
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
        panel.border = element_blank()) +
  
  geom_text(aes(label = 'Region'), x=2,y=500,size=30,color ='black')


png(file=paste(output_basepath,'\\Figures\\Donuts\\Dummy_marginSlopes.png',sep = ""), res = 300,width=3600,height=3600,bg = "transparent")
par(bty = 'n') 
print(pGM,axes=FALSE)
dev.off() 

RegionalSlopes <- function(regionDat,countMargin){
  steepLength <- vector()
  shallowLength <- vector()
  cliffLength <- vector()
  steepcliffLength <- vector()
  
  totalLength <- vector()
  
  for(subk in 1:length(regionDat)){ 
    
    SubBasinID <- regionDat[subk]
    offMargin <- read.csv(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_offmargin.csv',sep = ""))
    DEMList <- read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'ArcticDEMTiles.txt',sep = ""))
    
    sGrid <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'Grid.shp',sep = ""))
    sGrid <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'Grid.shp',sep = ""),layer = sGrid$layer)
   
    sGridClip <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'GridClipped.shp',sep = ""))
    sGridClip <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'GridClipped.shp',sep = ""),layer = sGridClip$layer)
    
    #Extract elevations from Cliff Cells
    #for(kD in 1:length(DEMList)){
    #  DEMH <- raster(DEMList[2,])
    #}
    
    terrMargin <- ogrInfo(paste(folderSectors,'\\',regionDat[subk],'Sector\\',regionDat[subk],'terrestrialmargin.shp',sep = ""))
    terrMargin <- readOGR(dsn=paste(folderSectors,'\\',regionDat[subk],'Sector\\',sep = ""),layer = terrMargin$layer)
    
    if(file.exists(paste(folderSectors,'\\',regionDat[subk],'Sector\\',regionDat[subk],'marinemargin.shp',sep = ""))){
    marMargin <- ogrInfo(paste(folderSectors,'\\',regionDat[subk],'Sector\\',regionDat[subk],'marinemargin.shp',sep = ""))
    marMargin <- readOGR(dsn=paste(folderSectors,'\\',regionDat[subk],'Sector\\',sep = ""),layer = marMargin$layer)
    if(regionDat[subk]==237){marMargin<-marMargin[,1]}
    if(exists('regionalMarMargin')){regionalMarMargin <- rbind(regionalMarMargin,marMargin)}else{regionalMarMargin <- marMargin}}
    
    if(file.exists(paste(folderSectors,'\\',regionDat[subk],'Sector\\',regionDat[subk],'lakemargin.shp',sep = ""))){
    lakMargin <- ogrInfo(paste(folderSectors,'\\',regionDat[subk],'Sector\\',regionDat[subk],'lakemargin.shp',sep = ""))
    lakMargin <- readOGR(dsn=paste(folderSectors,'\\',regionDat[subk],'Sector\\',sep = ""),layer = lakMargin$layer)
    if(exists('regionalLakMargin')){regionalLakMargin <- rbind(regionalLakMargin,lakMargin)}else{regionalLakMargin <- lakMargin}}
    
    
    if(subk>1){regionalTerrMargin <- rbind(regionalTerrMargin,terrMargin)}
    if(subk==1){regionalTerrMargin <- terrMargin}

    # Save Margins specific to type
    if(length(which(sGrid$type==1))>0){
      
      if(!is.na(mean(over(sGrid[which(sGrid$type==1),],terrMargin)[,1],na.rm=T))){
        cliffMargin <- intersect(sGrid[which(sGrid$type==1),],terrMargin)
      if(exists('regionalcliffMargin')){regionalcliffMargin <- rbind(regionalcliffMargin,cliffMargin)}
      if(!exists('regionalcliffMargin')){regionalcliffMargin <- cliffMargin}
    }}
    
    if(length(which(sGrid$type==2))>0){
      if(!is.na(over(sGrid[which(sGrid$type==2),],terrMargin))[1,1]){
      steepMargin <- intersect(sGrid[which(sGrid$type==2),],terrMargin)
      
      if(exists('regionalsteepMargin')){regionalsteepMargin <- rbind(regionalsteepMargin,steepMargin)}
      if(!exists('regionalsteepMargin')){regionalsteepMargin <- steepMargin}
    }}
    
    if(length(which(sGrid$type==3))>0){
      if(!is.na(over(sGrid[which(sGrid$type==3),],terrMargin))[1,1]){
      shallowMargin <- intersect(sGrid[which(sGrid$type==3),],terrMargin)
      
      if(exists('regionalshallowMargin')){regionalshallowMargin <- rbind(regionalshallowMargin,shallowMargin)}
      if(!exists('regionalshallowMargin')){regionalshallowMargin <- shallowMargin}
      }}
    
    if(length(which(sGrid$type==4))>0){
      if(!is.na(over(sGrid[which(sGrid$type==4),],terrMargin))[1,1]){
        steepCliffMargin <- intersect(sGrid[which(sGrid$type==4),],terrMargin)
        
        if(exists('regionalsteepCliffMargin')){regionalsteepCliffMargin <- rbind(regionalsteepCliffMargin,steepCliffMargin)}
        if(!exists('regionalsteepCliffMargin')){regionalsteepCliffMargin <- steepCliffMargin}
      }}
    
    if(length(which(sGrid$type==5))>0){
      if(!is.na(over(sGrid[which(sGrid$type==5),],terrMargin))[1,1]){
        ignoredMargin <- intersect(sGrid[which(sGrid$type==5),],terrMargin)
        
        if(subk>1&exists('regionalignoredMargin')){regionalignoredMargin <- rbind(regionalignoredMargin,ignoredMargin)}
        if(subk>1&!exists('regionalignoredMargin')){regionalignoredMargin <- ignoredMargin}
        if(subk==1){regionalignoredMargin <- ignoredMargin}
      }}
    
    if(length(offMargin$ID)>0){
      errorMargin <- intersect(sGrid[match(offMargin$ID,sGrid$id),],terrMargin)
      if(subk>1){regionalerrorMargin <- rbind(regionalerrorMargin,errorMargin)}
      if(subk==1){regionalerrorMargin <- errorMargin}
    }

    
    if(dim(   sGrid[which(sGrid$type==1),])[1]>0&!is.na(over(sGrid[which(sGrid$type==1),],terrMargin))[1,1]){
      cliffRamp <- intersect(terrMargin,sGrid[which(sGrid$type==1),])
      cliffLength[subk] <- gLength(cliffRamp)
    }else{cliffLength[subk] <- 0}
    
    if(dim(   sGrid[which(sGrid$type==4),])[1]>0&!is.na(over(sGrid[which(sGrid$type==4),],terrMargin))[1,1]){
      steepcliffRamp <- intersect(terrMargin,sGrid[which(sGrid$type==4),])
      steepcliffLength[subk] <- gLength(steepcliffRamp)
    }else{steepcliffLength[subk] <- 0}

    tryCatch({shallowRamp <- intersect(terrMargin,sGrid[which(sGrid$type==3),])
    shallowLength[subk] <- gLength(shallowRamp)
    },
             error = function(e){    shallowLength[subk] <- 0})
    
    if(dim(   sGrid[which(sGrid$type==2),])[1]>0&!is.na(over(sGrid[which(sGrid$type==2),],terrMargin))[1,1]){
      steepRamp <- intersect(terrMargin,sGrid[which(sGrid$type==2),])
      steepLength[subk] <- gLength(steepRamp)
    }else{steepLength[subk] <- 0}
    
    totalLength[subk] <- gLength(terrMargin)

    if(subk %% 1==0) {
      cat(paste0("region: ", floor(subk/length(regionDat)*100), "%"))}
  }
  
  cliffLength[which(is.na(cliffLength))] <- 0
  steepcliffLength[which(is.na(cliffLength))] <- 0
  steepLength[which(is.na(steepLength))] <- 0
  shallowLength[which(is.na(shallowLength))] <- 0
  
  dataMorphology <- cbind(regionDat,cliffLength, steepcliffLength,steepLength, shallowLength,totalLength)
  
  write.csv(dataMorphology,file=paste('C:\\Work\\Research\\tGISM\\Sectors\\',countMargin,'morphologystatistics.csv'))
  rgdal::writeOGR(regionalTerrMargin,paste('C:\\Work\\Research\\tGISM\\Sectors\\',countMargin,'RegionalTerrMargin.shp',sep = ""),layer = 'regionalTerrMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
  rgdal::writeOGR(regionalMarMargin,paste('C:\\Work\\Research\\tGISM\\Sectors\\',countMargin,'RegionalMarMargin.shp',sep = ""),layer = 'regionalMarMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
  rgdal::writeOGR(regionalLakMargin,paste('C:\\Work\\Research\\tGISM\\Sectors\\',countMargin,'RegionalLakMargin.shp',sep = ""),layer = 'regionalLakMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
  regionalerrorMargin <- regionalerrorMargin[,-5]
  rgdal::writeOGR(regionalerrorMargin,paste('C:\\Work\\Research\\tGISM\\Sectors\\',countMargin,'RegionalErrorMargin.shp',sep = ""),layer = 'regionalerrorMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
  rm(regionalerrorMargin)
  regionalignoredMargin <- regionalignoredMargin[,-5]
  rgdal::writeOGR(regionalignoredMargin,paste('C:\\Work\\Research\\tGISM\\Sectors\\',countMargin,'RegionalIgnoredMargin.shp',sep = ""),layer = 'regionalignoredMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
  rm(regionalignoredMargin)
  regionalcliffMargin <- regionalcliffMargin[,-5]
  rgdal::writeOGR(regionalcliffMargin,paste('C:\\Work\\Research\\tGISM\\Sectors\\',countMargin,'RegionalCliffMargin.shp',sep = ""),layer = 'regionalcliffMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
  rm(regionalcliffMargin)
  
  if(exists('regionalsteepCliffMargin')){regionalsteepCliffMargin <- regionalsteepCliffMargin[,-5]
  rgdal::writeOGR(regionalsteepCliffMargin,paste('C:\\Work\\Research\\tGISM\\Sectors\\',countMargin,'RegionalsteepCliffMargin.shp',sep = ""),layer = 'regionalsteepCliffMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
  rm(regionalsteepCliffMargin)}
  regionalsteepMargin <- regionalsteepMargin[,-5]
  rgdal::writeOGR(regionalsteepMargin,paste('C:\\Work\\Research\\tGISM\\Sectors\\',countMargin,'RegionalSteepMargin.shp',sep = ""),layer = 'regionalsteepMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
  rm(regionalsteepMargin)
  regionalshallowMargin <- regionalshallowMargin[,-5]
  rgdal::writeOGR(regionalshallowMargin,paste('C:\\Work\\Research\\tGISM\\Sectors\\',countMargin,'RegionalShallowMargin.shp',sep = ""),layer = 'regionalshallowMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
  rm(regionalshallowMargin)
  
  dataML <- data.frame(
    category = c('cliff','cliff2','steep','shallow','NA'),
    count=c(as.vector(colSums(cbind(cliffLength,steepcliffLength, steepLength, shallowLength,totalLength - (cliffLength+steepcliffLength+ steepLength+ shallowLength))))) 
  )
  dataML$fraction <- dataML$count / sum(dataML$count)
  dataML$ymax <- cumsum(dataML$fraction)
  dataML$ymin <- c(0, head(dataML$ymax, n=-1))
  dataML$labelPosition <- (dataML$ymax + dataML$ymin) / 2
  dataML$label <- paste0(round(dataML$fraction*100, digits = 1),"%")
  pGM <- ggplot(dataML, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect(color =colorMargins[countMargin],lwd=4) +
    #geom_text( x=c(4.3,4.3,4.3), aes(y=labelPosition, label=label), size=30) +
    scale_fill_manual(values=c(  "#ed0c0cff","#ed4c0cff","#BC8F8F","#319f28ff","#fc8101ff")) +
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
  
  print(sqrt(sum(dataML$count)/1000000)*10)
  print(dataML$label)
  print(sum(dataML$count)/1000)
  print(dataML$count)
  
  png(file=paste(output_basepath,'\\Figures\\Donuts\\',nameMargins[countMargin],'_marginSlopes.png',sep = ""), res = 300,width=3600,height=3600,bg = "transparent")
  par(bty = 'n') 
  print(pGM,axes=FALSE)
  dev.off() 
  
}

RegionalSlopes(regionCW, 1)
RegionalSlopes(regionNW, 2)
RegionalSlopes(regionNO, 3)
RegionalSlopes(regionNE, 4)
RegionalSlopes(regionCE, 5)
RegionalSlopes(regionSE, 6)
RegionalSlopes(regionSW, 7)
RegionalSlopes(PGICCW, 8)
RegionalSlopes(PGICNW, 9)
RegionalSlopes(PGICNOW, 10)
RegionalSlopes(PGICNOC, 11)
RegionalSlopes(PGICNOE, 12)
RegionalSlopes(PGICNEN, 13)
RegionalSlopes(PGICNEC, 14)
RegionalSlopes(PGICNES, 15)
RegionalSlopes(PGICCEN, 16)
RegionalSlopes(PGICCEC, 17)
RegionalSlopes(PGICCES, 18)
RegionalSlopes(PGICSEN, 19)
RegionalSlopes(PGICSEC, 20)
RegionalSlopes(PGICSES, 21)
RegionalSlopes(PGICSWS, 22)
RegionalSlopes(PGICSWN, 23)


# Combine margin fractions from all regions
for(i in 1:7){
  terrestrialMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalTerrMargin.shp',sep = ""))
  terrestrialMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalTerrMargin.shp',sep = ""),layer = terrestrialMar$layer)
  if(i==1){regionalTerrMar <- terrestrialMar}
  if(i>1){regionalTerrMar <- rbind(terrestrialMar,regionalTerrMar)}
  lakeMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalLakMargin.shp',sep = ""))
  lakeMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalLakMargin.shp',sep = ""),layer = lakeMar$layer)
  if(i==1){regionalLakeMar <- lakeMar}
  if(i>1){regionalLakeMar <- rbind(lakeMar,regionalLakeMar)}
  marineMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalMarMargin.shp',sep = ""))
  marineMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalMarMargin.shp',sep = ""),layer = marineMar$layer)
  if(i==1){regionalMarMar <- marineMar}
  if(i>1){regionalMarMar <- rbind(marineMar,regionalMarMar)}

  errMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalErrorMargin.shp',sep = ""))
  errMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalErrorMargin.shp',sep = ""),layer = errMar$layer)
  if(i==1){regionalErrMar <- errMar}
  if(i>1){regionalErrMar <- rbind(errMar,regionalErrMar)}
  
  cliffMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalCliffMargin.shp',sep = ""))
  cliffMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalCliffMargin.shp',sep = ""),layer = cliffMar$layer)
  if(i==1){regionalCliffMar <- cliffMar}
  if(i>1){regionalCliffMar <- rbind(cliffMar,regionalCliffMar)}
  
  steepcliffMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalsteepCliffMargin.shp',sep = ""))
  steepcliffMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalsteepCliffMargin.shp',sep = ""),layer = steepcliffMar$layer)
  if(i==1){regionalsteepCliffMar <- steepcliffMar}
  if(i>1){regionalsteepCliffMar <- rbind(steepcliffMar,regionalsteepCliffMar)}
  
  steepMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalSteepMargin.shp',sep = ""))
  steepMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalSteepMargin.shp',sep = ""),layer = steepMar$layer)
  if(i==1){regionalSteepMar <- steepMar}
  if(i>1){regionalSteepMar <- rbind(steepMar,regionalSteepMar)}
 
  shallowMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalShallowMargin.shp',sep = ""))
  shallowMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalShallowMargin.shp',sep = ""),layer = shallowMar$layer)
  if(i==1){regionalShallowMar <- shallowMar}
  if(i>1){regionalShallowMar <- rbind(shallowMar,regionalShallowMar)}
  
}
rgdal::writeOGR(regionalTerrMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalTerrMargin_GrIS.shp',sep = ""),layer = 'regionalTerrMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalMarMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalMarMargin_GrIS.shp',sep = ""),layer = 'regionalMarMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalLakeMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalLakMargin_GrIS.shp',sep = ""),layer = 'regionalLakMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalErrMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalErrorMargin_GrIS.shp',sep = ""),layer = 'regionalErrMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalCliffMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalCliffMargin_GrIS.shp',sep = ""),layer = 'regionalCliffMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalsteepCliffMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalsteepCliffMargin_GrIS.shp',sep = ""),layer = 'regionalsteepCliffMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalSteepMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalSteepMargin_GrIS.shp',sep = ""),layer = 'regionalSteepMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalShallowMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalShallowMargin_GrIS.shp',sep = ""),layer = 'regionalShallowMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
for(i in 8:23){
  terrestrialMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalTerrMargin.shp',sep = ""))
  terrestrialMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalTerrMargin.shp',sep = ""),layer = terrestrialMar$layer)
  if(i==8){regionalTerrMar <- terrestrialMar}
  if(i>8){regionalTerrMar <- rbind(terrestrialMar,regionalTerrMar)}
  lakeMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalLakMargin.shp',sep = ""))
  lakeMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalLakMargin.shp',sep = ""),layer = lakeMar$layer)
  if(i==8){regionalLakeMar <- lakeMar}
  if(i>8){regionalLakeMar <- rbind(lakeMar,regionalLakeMar)}
  marineMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalMarMargin.shp',sep = ""))
  marineMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalMarMargin.shp',sep = ""),layer = marineMar$layer)
  if(i==8){regionalMarMar <- marineMar}
  if(i>8){regionalMarMar <- rbind(marineMar,regionalMarMar)}
  
  errMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalErrorMargin.shp',sep = ""))
  errMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalErrorMargin.shp',sep = ""),layer = errMar$layer)
  if(i==8){regionalErrMar <- errMar}
  if(i>8){regionalErrMar <- rbind(errMar,regionalErrMar)}
  
  cliffMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalCliffMargin.shp',sep = ""))
  cliffMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalCliffMargin.shp',sep = ""),layer = cliffMar$layer)
  if(i==8){regionalCliffMar <- cliffMar}
  if(i>8){regionalCliffMar <- rbind(cliffMar,regionalCliffMar)}
  
  steepcliffMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalsteepCliffMargin.shp',sep = ""))
  steepcliffMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalsteepCliffMargin.shp',sep = ""),layer = steepcliffMar$layer)
  if(i==8){regionalsteepCliffMar <- steepcliffMar}
  if(i>8){regionalsteepCliffMar <- rbind(steepcliffMar,regionalsteepCliffMar)}
  
  steepMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalSteepMargin.shp',sep = ""))
  steepMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalSteepMargin.shp',sep = ""),layer = steepMar$layer)
  if(i==8){regionalSteepMar <- steepMar}
  if(i>8){regionalSteepMar <- rbind(steepMar,regionalSteepMar)}
  
  shallowMar <- ogrInfo(paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalShallowMargin.shp',sep = ""))
  shallowMar <-   readOGR(dsn=paste('C:\\Work\\Research\\tGISM\\Sectors\\',i,'RegionalShallowMargin.shp',sep = ""),layer = shallowMar$layer)
  if(i==8){regionalShallowMar <- shallowMar}
  if(i>8){regionalShallowMar <- rbind(shallowMar,regionalShallowMar)}
  
}
rgdal::writeOGR(regionalTerrMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalTerrMargin_PGIC.shp',sep = ""),layer = 'regionalTerrMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalMarMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalMarMargin_PGIC.shp',sep = ""),layer = 'regionalMarMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalLakeMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalLakMargin_PGIC.shp',sep = ""),layer = 'regionalLakMar',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalErrMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalErrorMargin_PGIC.shp',sep = ""),layer = 'regionalErrMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalCliffMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalCliffMargin_PGIC.shp',sep = ""),layer = 'regionalCliffMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalsteepCliffMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalsteepCliffMargin_PGIC.shp',sep = ""),layer = 'regionalsteepCliffMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalSteepMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalSteepMargin_PGIC.shp',sep = ""),layer = 'regionalSteepMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(regionalShallowMar,paste('C:\\Work\\Research\\tGISM\\Sectors\\RegionalShallowMargin_PGIC.shp',sep = ""),layer = 'regionalShallowMargin',driver = 'ESRI Shapefile',overwrite_layer = T)
