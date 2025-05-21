################################################################################
# Validate steep ice margin with high resolution DEMs from field sites
#
# ArctciDEM_Margin_Extractor.R
#
# ReadMe: 
#
# https://github.com/fidelsteiner/tGISM
#
# Input:
# Validation on two sites (Red Rock, Mittivakat) by comparing the Arctic DEM against available Pleiades imagery
# 
#
# Created:          2021/10/15
# Latest Revision:  2023/01/23
#
# Jakob F Steiner | jakob@x-hydrolab.org | x-hydrolab.org 
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
library(hydroGOF)


##################
# Multiple subsites for analysis of the margin
##################
GlacierName <- 'RedRock'

fnDEM_1 <- 'PGO_2017-08-15_RedRockCliff_DEM_2m_StereoProj.tif'
fnDEM_2 <- 'c0_reg_dem_ProjectRaster1.tif'

fnDEM_MV_1 <- 'E:\\GeospatialData\\DEMData\\ValidationDEMs\\Mittivakat\\transfer_20608_files_611c06f9\\Mittivakat_Merged.tif'
fnDEM_MV_2 <- 'E:\\GeospatialData\\DEMData\\ValidationDEMs\\Mittivakat\\transfer_20608_files_611c06f9\\2018-08-09_Mittivakkat\\PGO\\2018-08-09_Mittivakkat\\PGO_2018-08-09_Mittivakkat_DEM_2m.tif'
fnDEM_MV_3 <- 'E:\\GeospatialData\\DEMData\\ValidationDEMs\\Mittivakat\\transfer_20608_files_611c06f9\\2018-08-15_Mittivakkat\\PGO\\2018-08-15_Mittivakkat\\PGO_2018-08-15_Mittivakkat_DEM_2m.tif'

fnDEM_MV_Arctic1 <- 'E:\\ArcticDEM\\DEMs\\14_44_1_1_2m_v4.1_dem.tif'
fnDEM_MV_Arctic2 <- 'E:\\ArcticDEM\\DEMs\\14_44_1_2_2m_v4.1_dem.tif'
fnDEM_MV_Arctic3 <- 'E:\\ArcticDEM\\DEMs\\14_44_2_1_2m_v4.1_dem.tif'
fnDEM_MV_Arctic4 <- 'E:\\ArcticDEM\\DEMs\\14_44_2_2_2m_v4.1_dem.tif'


fnDEM_RR_Arctic1 <- 'E:\\ArcticDEM\\DEMs\\27_35_2_2_2m_v4.1_dem.tif'

fnMarginOutline_RR <- 'RedRock_ManualBuffer_200m'
fnMarginLine_RR <- 'RedRock_ManualLine'

fnMarginOutline_MV <- 'Mittivakat_ManualBuffer_200m'
fnMarginLine_MV <- 'Mittivakat_ManualLine'

##################
# Paths
##################

path_RR <-'C:\\Work\\Research\\tGISM\\MarginValidation\\ArcticDEMValidation\\RedRock'
path_MV <-'C:\\Work\\Research\\tGISM\\MarginValidation\\ArcticDEMValidation\\Mittivakat'

path_figs<-'C:\\Work\\Research\\tGISM\\Figures'
pathDEM_PleiadesRR<-'E:\\GeospatialData\\DEMData\\ValidationDEMs\\RedRock_Pleiades\\Pleiades\\DEM\\2017-08-15_RedRockCliff\\PGO\\2017-08-15_RedRockCliff'
pathDEM_PleiadesMV_1<-'E:\\GeospatialData\\DEMData\\ValidationDEMs\\Mittivakat\\transfer_20608_files_611c06f9\\2018-08-08_Mittivakkat\\PGO\\2018-08-08_Mittivakkat'
pathDEM_PleiadesMV_2<-'E:\\GeospatialData\\DEMData\\ValidationDEMs\\Mittivakat\\transfer_20608_files_611c06f9\\2018-08-09_Mittivakkat\\PGO\\2018-08-09_Mittivakkat'
pathDEM_PleiadesMV_3<-'E:\\GeospatialData\\DEMData\\ValidationDEMs\\Mittivakat\\transfer_20608_files_611c06f9\\2018-08-15_Mittivakkat\\PGO\\2018-08-15_Mittivakkat'

pathDEM_Arctic<-'G:\\ArcticDEM'

projec<-'+proj=utm +datum=WGS84'

projec_utm<-'+proj=utm +zone=19N +datum=WGS84'
projec_27N<-'+proj=utm +zone=27N +datum=WGS84'

##################
# Red Rock Site
##################

# read Margin line
ogrInfo(path_RR,fnMarginLine_RR)
margin_pEXTL<-readOGR(dsn=path_RR,layer=fnMarginLine_RR)
#projection(margin_pEXT)<-projec
marginValLength <- rgeos::gLength(margin_pEXTL)

# read Margin Buffer
ogrInfo(path_RR,fnMarginOutline_RR)
margin_pEXT<-readOGR(dsn=path_RR,layer=fnMarginOutline_RR)
#projection(margin_pEXT)<-projec
margin_pEXT<-SpatialPolygons(margin_pEXT@polygons,proj4string=margin_pEXT@proj4string)

#######################
# Read DEM Data
#######################
DEM_temp <- raster(pathDEM_PleiadesRR&'/'&fnDEM_1)
#projection(DEM_temp) <- projec_utm

DEM_Arctic <- raster(fnDEM_RR_Arctic1)
#projection(DEM_Arctic) <- projec_utm

margin_pEXT <- spTransform(margin_pEXT,CRS('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))

DEM_masked <- mask(crop(DEM_temp,margin_pEXT),margin_pEXT)
DEMArctic_masked <- mask(crop(DEM_Arctic,margin_pEXT),margin_pEXT)

RR_slope <- terrain(DEM_masked, opt="slope", unit="degrees", neighbors=8)

RR_slopeArctic <- terrain(DEMArctic_masked, opt="slope", unit="degrees", neighbors=8)

set.seed(1234)
wdata_IceCap = data.frame(
  DEM = factor(rep(c("ArcticDEM", "Pleiades"), each=length(which(!is.na(RR_slope[]))))),
  slope = c(RR_slopeArctic[which(!is.na(RR_slope[]))], RR_slope[which(!is.na(RR_slope[]))])
)
head(wdata_IceCap, 4)

mu <- wdata_IceCap %>% 
  group_by(DEM) %>%
  summarise(grp.mean = median(slope,na.rm=T))

png(file=path_figs&'\\RedRock_ICeCap_DEMmatch.png', res = 300,width=3600,height=900)
par(mar=c(4,2,5,1),cex.lab=2,cex.axis=2)
ggplot(wdata_IceCap, aes(x = slope, y = DEM)) +
  geom_density_ridges(aes(fill = DEM)) +
  scale_fill_manual(values = c("#E7B800","#00AFBB")) +
  geom_vline(aes(xintercept = grp.mean, color = DEM),
             data = mu, linetype = "dashed",size = 1.2) +
  scale_x_continuous(name=expression('Slope [°]'), limits=c(0, 85)) +
  theme(axis.text.x = element_text(size=24, angle=270),
        axis.text.y = element_text(size=24, angle=0)) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank()) +
  theme(legend.text=element_text(size=rel(1)),legend.position="none")
theme_bw()
dev.off()


#ECDFs
ecdfArctic <- ecdf(RR_slopeArctic[])
ecdfPleiades <- ecdf(RR_slope[])

png(file=path_figs&'\\ECDF_RedRock.png', res = 160,width=1200,height=900)
par(mar=c(5,5,1,1),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE))
plot(sort(RR_slopeArctic[]),ecdfArctic(sort(RR_slopeArctic[])),type='l',xlim=c(0, 90),lwd=2, main='', panel.first=grid(equilogs=T),ylab = 'probability of exceedance [-]',xlab = 'slope [°]',do.points=F,verticals = T,col.01line = NULL,ylim=c(0,1))
points(sort(RR_slope[]),ecdfPleiades(sort(RR_slope[])),type='l',col='red',lwd=2,col.01line = NULL)
legend('bottomright',bty='n',col=c('black','red'),lty = 1,lwd = 2,legend=c('ArcticDEM', 'Pléiades'),cex=1.5)
dev.off()

# error statistics between slope rasters
RR_slope_res <- resample(RR_slope,RR_slopeArctic)

medPleiades <- median(RR_slope[],na.rm=T)
medArctic <- median(RR_slopeArctic[],na.rm=T)

medDev <- median(RR_slope_res[] - RR_slopeArctic[],na.rm=T)
biasDev <- (sum(RR_slope_res[] - RR_slopeArctic[],na.rm=t))/length(which(!is.na(RR_slope_res[] - RR_slopeArctic[])))
NSEDev <- NSE(RR_slope_res[],RR_slopeArctic[],na.rm=T,fun=NULL)
corDEV <- cor(RR_slope_res[which(!is.na(RR_slope_res[] - RR_slopeArctic[]))],RR_slopeArctic[which(!is.na(RR_slope_res[] - RR_slopeArctic[]))],method = 'spearman')

#######################
# Read Mittivakkat Site
#######################

# read Margin line
ogrInfo(path_MV,fnMarginLine_MV)
margin_pEXTL<-readOGR(dsn=path_MV,layer=fnMarginLine_MV)
#projection(margin_pEXT)<-projec
marginValLength <- rgeos::gLength(margin_pEXTL)

# read Margin Buffer
ogrInfo(path_MV,fnMarginOutline_MV)
margin_pEXT<-readOGR(dsn=path_MV,layer=fnMarginOutline_MV)
#projection(margin_pEXT)<-projec
margin_pEXT<-SpatialPolygons(margin_pEXT@polygons,proj4string=margin_pEXT@proj4string)

#######################
# Read DEM Data
#######################
DEM_temp_1 <- raster(fnDEM_MV_1)
#DEM_temp_2 <- raster(fnDEM_MV_2)
#DEM_temp_3 <- raster(fnDEM_MV_3)
#DEM_temp <- merge(DEM_temp_1,DEM_temp_2,DEM_temp_3)
#projection(DEM_temp) <- projec_utm

DEM_Arctic1 <- raster(fnDEM_MV_Arctic1)
DEM_Arctic3 <- raster(fnDEM_MV_Arctic3)
#projection(DEM_Arctic) <- projec_utm
margin_pEXT <- spTransform(margin_pEXT,CRS('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))
DEM_Arctic1 <- crop(DEM_Arctic1,margin_pEXT)
DEM_Arctic3 <- crop(DEM_Arctic3,margin_pEXT)
DEM_Arctic <- merge(DEM_Arctic1,DEM_Arctic3) 


DEM_masked <- mask(crop(DEM_temp_1,margin_pEXT),margin_pEXT)
MV_slope <- terrain(DEM_masked, opt="slope", unit="degrees", neighbors=8)

DEMArctic_masked <- mask(crop(DEM_Arctic,margin_pEXT),margin_pEXT)
MV_slopeArctic <- terrain(DEMArctic_masked, opt="slope", unit="degrees", neighbors=8)

set.seed(1234)
wdata_IceCap = data.frame(
  DEM = factor(rep(c("ArcticDEM", "Pleiades"), each=length(which(!is.na(MV_slope[]))))),
  slope = c(MV_slopeArctic[which(!is.na(MV_slope[]))], MV_slope[which(!is.na(MV_slope[]))])
)
head(wdata_IceCap, 4)

mu <- wdata_IceCap %>% 
  group_by(DEM) %>%
  summarise(grp.mean = median(slope,na.rm=T))

png(file=path_figs&'\\Mittivakat_ICeCap_DEMmatch.png', res = 300,width=3600,height=900)
par(mar=c(4,2,5,1),cex.lab=1.2,cex.axis=1.5)
ggplot(wdata_IceCap, aes(x = slope, y = DEM)) +
  geom_density_ridges(aes(fill = DEM)) +
  scale_fill_manual(values = c("#E7B800","#00AFBB")) +
  geom_vline(aes(xintercept = grp.mean, color = DEM),
             data = mu, linetype = "dashed",size = 1.2) +
  scale_x_continuous(name=expression('Slope [°]'), limits=c(0, 80)) +
  theme(axis.text.x = element_text(size=24, angle=270),
        axis.text.y = element_text(size=24, angle=0)) +
  theme(axis.title.y =element_blank(),
        axis.text.y=element_blank()) +
  theme(legend.text=element_text(size=rel(1)),legend.position="none")
theme_bw()
dev.off()


#ECDFs
ecdfArctic <- ecdf(MV_slopeArctic[])
ecdfPleiades <- ecdf(MV_slope[])

png(file=path_figs&'\\ECDF_Mittivakat.png', res = 160,width=1200,height=900)
par(mar=c(5,5,1,1),cex.lab=1.5,cex.axis=1.5)
layout(matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE))
plot(sort(MV_slopeArctic[]),ecdfArctic(sort(MV_slopeArctic[])),type='l',xlim=c(0, 90),lwd=2, main='', panel.first=grid(equilogs=T),ylab = 'probability of exceedance [-]',xlab = 'slope [°]',do.points=F,verticals = T,col.01line = NULL,ylim=c(0,1))
points(sort(MV_slope[]),ecdfPleiades(sort(MV_slope[])),type='l',col='red',lwd=2,col.01line = NULL)
legend('bottomright',bty='n',col=c('black','red'),lty = 1,lwd = 2,legend=c('ArcticDEM', 'Pléiades'),cex=1.5)
dev.off()

# error statistics between slope rasters
MV_slope_res <- resample(MV_slope,MV_slopeArctic)

medPleiades <- median(MV_slope[],na.rm=T)
medArctic <- median(MV_slopeArctic[],na.rm=T)

medDev <- median(MV_slope_res[] - MV_slopeArctic[],na.rm=T)
biasDev <- (sum(MV_slope_res[] - MV_slopeArctic[],na.rm=t))/length(which(!is.na(MV_slope_res[] - MV_slopeArctic[])))
NSEDev <- NSE(MV_slope_res[],MV_slopeArctic[],na.rm=T,fun=NULL)
corDEV <- cor(MV_slope_res[which(!is.na(MV_slope_res[] - MV_slopeArctic[]))],MV_slopeArctic[which(!is.na(MV_slope_res[] - MV_slopeArctic[]))],method = 'spearman')
