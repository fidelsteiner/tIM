################################################################################
# Spatial Analysis of the Ice Margin over Greenland / Global Code
#
#
# ReadMe:
# https://github.com/fidelsteiner/tGISM
#
# Created:          2018/02/05
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
#install.packages('pacman')
library(pacman)
p_load(spatialEco,ncdf4,gstat,rgdal,rgeos,maptools,raster,foreach,lubridate,compare,colorRamps,data.table,circular,parallel,snowfall,truncnorm,rlecuyer,forecast,rasterVis,R.utils,zoo,rlist,ggplot2)
library(sf)
library(sp)

##################
# paths
##################

data_basepath <- 'E:\\GeospatialData'                     # path for heavy input data
output_basepath <- 'C:\\Work\\Research\\tGISM'            # path for code/output   

path_outlines <- paste(data_basepath,'\\PROMICE_database\\Glacier Outline',sep = "")      # PROMICE Data (margin)
fnMarginOutline <- 'PROMICdiss'
path_outlines2 <- paste(data_basepath,'\\Outlines\\',sep = "")  # NOT IN USE REMOVE
path_cci <- paste(data_basepath,'\\cci_data\\outlines',sep = "")                          # CCI Outline data Ice Sheet
fnMarginOutline_cci <- 'cci_go_greenland_icesheet_2000'
path_outlines_RGI <- paste(data_basepath,'\\Periphery',sep = "")                          # RGI outlines          
fnMarginOutline_RGI <- '05_rgi60_GreenlandPeriphery'
pathDEM_topo <- paste(data_basepath,'\\Bathymetry',sep = "")                              # GEBCO Bathymetry data
path_figs <- paste(data_basepath,'\\Figures',sep = "")                                    # Output path for figures
path_mankoff <- paste(data_basepath,'\\PROMICE_database\\MankoffData\\sector_D.csv',sep = "")
path_mouginot <- 'E:\\Research\\OtherResearch\\Greenland\\Mouginot\\'      # not used REMOVE
path_noel <- paste(data_basepath,'\\NoelData',sep = "")                                   # SMB Data

path_sectorMargins <- paste(data_basepath,'\\Sectors\Sectors\SectorMargins',sep = "")     # split margins across all sectors

# Bedmachine data
bedmachineRAW <- 'E:\\Research\\OtherResearch\\Greenland\\Bedmachine\\5000000451838\\160281892\\BedMachineGreenland-2017-09-20.nc'

# basin data
GLbasins <- 'E:\\Research\\OtherResearch\\Greenland\\Mouginot\\drainagebasins\\GRE_Basins_IMBIE2_v1.3.shp'
# sector data
GLsectors <- 'E:\\Research\\OtherResearch\\Greenland\\Mouginot\\drainagebasins\\MouginotRignot\\doi_10.7280_D1WT11__v1\\Greenland_Basins_PS_v1.4.2.shp'

#read PROMICE margin  including all ice caps (1980s)
ogrInfo(path_outlines,fnMarginOutline)
margin_pEXT<-readOGR(dsn=path_outlines,layer=fnMarginOutline)
margin_pEXT<-SpatialPolygons(margin_pEXT@polygons,proj4string=margin_pEXT@proj4string)
projection(margin_pEXT)<-CRS("+init=epsg:4326")
margin_pEXT_transf<-spTransform(margin_pEXT,crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))

margin_actual <- margin_pEXT_transf[1]
margin_glaciers <- margin_pEXT_transf[2]
margin_peary <- margin_pEXT_transf[3]
margin_icecaps <- margin_pEXT_transf[4]

#read CCI Margin (Rastner2012)
margin_CCI_sf = sf::read_sf(path_cci&"\\cci_go_greenland_icesheet_2000.shp")  # assumes in project root
margin_CCI_sf = sf::st_transform(margin_CCI_sf, crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))

margin_CCI <- sf::as_Spatial(st_geometry(margin_CCI_sf), IDs = as.character(1:nrow(margin_CCI_sf)))
margin_CCI_nonunataks <- remove.holes(st_as_sf(margin_CCI))     # remove nunataks from the dataset

#read Rastner margin (Rastner2012)
ogrInfo(path_outlines_RGI,fnMarginOutline_RGI)
margin_RGI<-readOGR(dsn=path_outlines_RGI,layer=fnMarginOutline_RGI)
margin_RGI<-SpatialPolygons(margin_RGI@polygons,proj4string=margin_RGI@proj4string)
projection(margin_RGI)<-CRS("+init=epsg:4326")
margin_RGI_transf<-spTransform(margin_RGI,crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))
margin_RGI_transf <- margin_RGI_transf %>% gBuffer(byid=TRUE, width=0)  # clean intersections/loops

# read basin outlines
ogrInfo(GLbasins)
basins_pEXT<-readOGR(dsn=GLbasins)
basins_pEXT<-SpatialPolygons(basins_pEXT@polygons,proj4string=basins_pEXT@proj4string)
projection(basins_pEXT)<-CRS("+init=epsg:4326")
basins_pEXT_transf<-spTransform(basins_pEXT,crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))

##################
#Remove marginal glaciers from the PROMICE dataset (multiple hours of computing on single core)
if(file.exists(path_outlines_RGI&'\\PROMICE_Rastner_clean.shp')){
  ogrInfo(path_outlines_RGI,'PROMICE_Rastner_clean')
  cleanmargin <- readOGR(dsn=path_outlines_RGI,layer='PROMICE_Rastner_clean')
  cleanmargin <- SpatialPolygons(cleanmargin@polygons,proj4string=cleanmargin@proj4string)
} else{
croppedmargin <- rgeos::gDifference(margin_actual,margin_RGI_transf)
writeOGR(as(croppedmargin, "SpatialPolygonsDataFrame"), dsn=path_outlines_RGI&'\\PROMICE_Rastner.shp',layer=1, driver="ESRI Shapefile")

croppedmargin_decoupled <- disaggregate(croppedmargin)
cleanmargin <- croppedmargin_decoupled[which.max(area(croppedmargin_decoupled))] # Remove small 'islands'
writeOGR(as(cleanmargin, "SpatialPolygonsDataFrame"), dsn=path_outlines_RGI&'\\PROMICE_Rastner_clean.shp',layer=1, driver="ESRI Shapefile")
}
cleanmargin_nonunataks <- remove.holes(cleanmargin)   # remove nunataks from the dataset
margin_line <- as(cleanmargin,'SpatialLines') # Convert to Spatial Line for processing
margin_line_nonuna <- as(cleanmargin_nonunataks,'SpatialLines') # Convert to Spatial Line for processing

# read sector outlines (Mouginot2019)
ogrInfo(GLsectors)
sectors_pEXT<-readOGR(dsn=GLsectors)
sectors_pEXT<-SpatialPolygons(sectors_pEXT@polygons,proj4string=sectors_pEXT@proj4string)
projection(sectors_pEXT)<-CRS("+init=epsg:4326")
#sectors_pEXT_transf<-spTransform(sectors_pEXT,crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))
projection(sectors_pEXT)<-crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')

# Read Bedmachine (Citation!)
bedmachine_bed <- raster(bedmachineRAW, varname ='bed')
bedmachine_bed[bedmachine_bed<=-1000]<-NA
projection(bedmachine_bed) <- '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

coarse_bed <- aggregate(bedmachine_bed,fact=1)

################### RACMO data from Brice Noel (Noel2019/2016)
################### 
fullRACMO1km <- raster(path_noel&'\\SMB\\smb_rec.1958-2019.BN_RACMO2.3p2_FGRN055_GrIS.MM.nc')
fullRACMOgrid <- raster(path_noel&'\\Icemask_Topo_Iceclasses_lon_lat_average_1km_GrIS.nc',varname='Promicemask')


################### 
# read and process velocities (Mouginot2019)
################### 
ncin_x <- raster(path_mouginot&'Velocities\\2010_2017\\vel_2016-07-01_2017-06-31.nc',varname='VX')
ncin_y <- raster(path_mouginot&'Velocities\\2010_2017\\vel_2016-07-01_2017-06-31.nc',varname='VY')
velmap <- sqrt(ncin_x^2 + ncin_y^2)
projection(velmap)<-crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
projection(ncin_x)<-crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
projection(ncin_y)<-crs('+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')
velmap[velmap[]<0] <- NA
velmap_bed <- raster::resample(velmap,coarse_bed)
xvel <- raster::resample(ncin_x,coarse_bed)
yvel <- raster::resample(ncin_y,coarse_bed)

dir_v <- atan2(yvel,xvel) * 180 / pi
velmap_bed[coarse_bed[]<0]<-NA
xvel[coarse_bed[]<0]<-NA
yvel[coarse_bed[]<0]<-NA
dir_v[coarse_bed[]<0]<-NA

writeRaster(xvel, path_outlines2&'\\dumper\\xvel.tif','GTiff',overwrite=T)
writeRaster(yvel, path_outlines2&'\\dumper\\yvel.tif','GTiff',overwrite=T)
writeRaster(dir_v, path_outlines2&'\\dumper\\dir_v.tif','GTiff',overwrite=T)
writeRaster(velmap_bed, path_outlines2&'\\dumper\\mag.tif','GTiff',overwrite=T)

# extract velocities along margin
velmap_coarse <- aggregate(velmap,fact= 20)
bedmachine_coarse <- aggregate(bedmachine_bed,fact=20)
velmap_bed <- raster::resample(bedmachine_coarse,velmap_coarse)
margin_line <- as(margin_actual,'SpatialLines')
velMarg <- rasterize(margin_line,velmap_bed)
velMarg_holed <- boundaries(velMarg, type='inner', classes=FALSE, directions=8,asNA=T)

velMarg_holed[which(velMarg_holed[]>0&!is.na(velMarg_holed[]))] <- 2
velMarg_holed[which(bedmachine_coarse[]<=0&!is.na(velMarg_holed[]))] <- 3



# Read Mankoff Fluxes
s.sf <- st_read(GLsectors)
s.sf$NAME
mankoffFlux <- read.csv(path_mankoff,sep=',',header=T)

fluxVal <- vector()
for(i in 1:length(s.sf$NAME)){
  
  if(c(levels(s.sf$NAME)[i]) %in% colnames(mankoffFlux)==TRUE){
    fluxVal[i] <- mean(mankoffFlux[,c(levels(s.sf$NAME)[i])])
  }
  else{
    fluxVal[i] <- NA
  }
}


# Read Mankoff Fluxes
p<-sectors_pEXT

( pid <- sapply(slot(p, "polygons"), function(x) slot(x, "ID")) )
( p.df <- data.frame( ID=1:length(p), row.names = pid) )  
p <- SpatialPolygonsDataFrame(p, p.df)
p@data$flux <- fluxVal
colVals <- round(p@data$flux*length(pal(100)))
colVals[is.na(colVals)]<-1
colRamp <- c('#ffffff',col=pal(99))

png(file=path_figs&'\\iceDischarge.png', res = 160,width=1800,height=1800)
par(mar=c(2,2,2,2),cex.lab=1.2,cex.axis=1.2)
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(p,col=colRamp[colVals])
legend_image <- as.raster(matrix(colRamp, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'discharge [Gt yr-1]')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0,38,l=5))
rasterImage(rev(legend_image[0:100]), 0, 0, 1,1)
dev.off()

##############
# Overlay subbasins on margin
##############

# Test one subbasin
ltermFrac <- matrix(NA,nrow=length(sectors_pEXT),ncol=2)

#for(i in 1:length(sectors_pEXT))
testfunc <- function(i){
#if(is.na(over(sectors_pEXT[i],margin_RGI_transf))){
sector_line <- as(sectors_pEXT[i],'SpatialLines')

if(!is.null(sector_line)){
cleanmargin_line <- as(margin_CCI,'SpatialLines')
sector_buff <- buffer(sectors_pEXT[i], width = 1000, dissolve = T)
cleanmargin_sector <- crop(cleanmargin_line,extent(sector_buff)+c(1000*c(1,1,-1,-1)))

cleanmargin_line_nonuna <- as(margin_CCI_nonunataks,'SpatialLines')
cleanmargin_sector_nonuna <- crop(cleanmargin_line_nonuna,extent(sector_buff)+c(1000*c(1,1,-1,-1)))

#extraFrame <- c(-1,1,-1,1)
#if(extent(cleanmargin_sector)[2]<0){extraFrame[2] <- 1}
#if(extent(cleanmargin_sector)[3]<0){extraFrame[3] <- 1}

#opt<-crop(sector_line,extent(cleanmargin_sector)+c(1000*extraFrame))
#buff<- buffer(cleanmargin_sector, width = 1000, dissolve = T)
#opt2 <- raster::intersect(opt,buff)

#buff1 <- buffer(sector_line, width = 2000, dissolve = T)
#buff1 <- remove.holes(buff1)
#marginCCIcrop <- crop(margin_CCI_nonunataks,extent(sector_buff)+c(1000*extraFrame))
                      
#opp<-gDifference(buff1,marginCCIcrop)

if(!is.null(cleanmargin_sector)){
margin_sector <- raster::intersect(cleanmargin_sector,sector_buff)
row.names(margin_sector) <- "1"

df<-SpatialLinesDataFrame(margin_sector, data.frame(id=1:length(margin_sector)))
writeOGR(df, dsn=path_sectorMargins&'\\marginSector'&i&'.shp' ,layer="id",driver="ESRI Shapefile")

#marg_ras <- rasterize(margin_sector,coarse_bed)
#marg_ras[which(coarse_bed[]>0&!is.na(marg_ras[]))] <- 2
#marg_ras[which(coarse_bed[]<=0&!is.na(marg_ras[]))] <- 3
#marg_ras[which(marg_ras[]==1)] <- NA
#writeRaster(marg_ras, path_sectorMargins&'\\marginSector'&i&'.tif','GTiff',overwrite=T)

#ltermFrac[i,1] <- length(which(marg_ras[]==2)) / length(which(marg_ras[]>1))

#tabData <- as.table(ltermFrac[i,1])
#rownames(tabData) <- c('fraction')

#write.table(tabData, file = path_sectorMargins&'\\temp\\marginSector'&i&'.txt')

}
if(!is.null(cleanmargin_sector_nonuna)){
  margin_sector_nonuna <- raster::intersect(cleanmargin_sector_nonuna,sector_buff)
  row.names(margin_sector_nonuna) <- "1"
  
  df<-SpatialLinesDataFrame(margin_sector_nonuna, data.frame(id=1:length(margin_sector_nonuna)))
  writeOGR(df, dsn=path_sectorMargins&'\\marginSector_nonuna'&i&'.shp' ,layer="id",driver="ESRI Shapefile")

  #marg_ras_nonuna <- rasterize(margin_sector_nonuna,coarse_bed)
  #marg_ras_nonuna[which(coarse_bed[]>0&!is.na(marg_ras_nonuna[]))] <- 2
  #marg_ras_nonuna[which(coarse_bed[]<=0&!is.na(marg_ras_nonuna[]))] <- 3
  #marg_ras_nonuna[which(marg_ras_nonuna[]==1)] <- NA
  #writeRaster(marg_ras_nonuna, path_sectorMargins&'\\marginSector_nonuna'&i&'.tif','GTiff',overwrite=T)
  
  #ltermFrac[i,2] <- length(which(marg_ras_nonuna[]==2)) / length(which(marg_ras_nonuna[]>1))
  
  #tabData <- as.table(ltermFrac[i,2])
  #rownames(tabData) <- c('fraction')

  #write.table(tabData, file = path_sectorMargins&'\\temp\\marginSector_nonuna'&i&'.txt')
}
}
#}
print(i)
}

# number of workers (leave one to keep system more responsive)
nworkers <- detectCores()-4

# Run Model and save results in blocks to 

  # run entire function parallelized
sfInit(parallel=T,cpu=nworkers, type="SOCK")                                                # initialize cluster
loadlib <- lapply(c('raster','truncnorm','rgdal'),                                                  # load required packages in cluster
                    function(x) sfLibrary(x, character.only=T))
sfExportAll()                                                                               # transfer all variables to clusters
sfClusterSetupRNG()                                                                         # set up random number generator
  
#ModelRunSeries <- seq(1,length(sectors_pEXT),1)[-c(6,7,53,54,106,109,159,160,211,212)];
ModelRunSeries <- c(253,seq(234,237,1),205,206,184,seq(177,179,1),seq(153,158,1),99)
tout.multi <- system.time(
    out <- sfLapply(ModelRunSeries, testfunc)
)
sfStop() 



p<-sectors_pEXT

( pid <- sapply(slot(p, "polygons"), function(x) slot(x, "ID")) )
( p.df <- data.frame( ID=1:length(p), row.names = pid) )  
p <- SpatialPolygonsDataFrame(p, p.df)
ltermFrac2<- seq(1,260,1)*0 + 1
ltermFrac2[1:length(ltermFrac[,1])] <- ltermFrac[,1]
p@data$frac <- ltermFrac2
pal <- colorRampPalette(c("white", "red"))

colRamp <- c('#ffffff',col=pal(99))
colVals <- round(p@data$frac*length(pal(100)))


png(file=path_figs&'\\relativeTerrestrialMargin.png', res = 160,width=1800,height=1800)
par(mar=c(2,2,2,2),cex.lab=1.2,cex.axis=1.2)
layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
plot(p,col=colRamp[colVals])
legend_image <- as.raster(matrix(colRamp, ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'terrestrial fraction')
text(x=1.5, y = seq(0,1,l=5), labels = seq(0.3,1,l=5))
rasterImage(rev(legend_image[30:100]), 0, 0, 1,1)
dev.off()

p@data$flux <- fluxVal
colVals <- round(p@data$flux*length(pal(100)))
colVals[is.na(colVals)]<-1

pal <- terrain.colors(1000)

png(file=path_figs&'\\IceMargin.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
#layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))
plot(bedmachine_bed,xlab = 'Westing [m]', ylab ='Northing [m]',legend.args=list(text=expression('Elevation [m]')),col=pal)
plot(margin_line,add=TRUE)
grid(NULL,NULL)
dev.off()

# Rasterize Margin to Bedmachine
rasterized_margin_icesheet <- rasterize(margin_actual,coarse_bed)
# Remove the inside cells (ice sheet)
rasterized_margin_holed <- rasterized_margin_icesheet
rasterized_margin_holed <- boundaries(rasterized_margin_icesheet, type='inner', classes=FALSE, directions=8,asNA=T)
writeRaster(rasterized_margin, path_outlines2&'\\RasterizedMargin_icesheet_promicee150m.tif','GTiff')


rasterized_margin_holed[which(rasterized_margin_holed[]==0)] <- NA
rasterized_margin_holed[which(coarse_bed[]>0&!is.na(rasterized_margin_holed[]))] <- 2
rasterized_margin_holed[which(coarse_bed[]<=0&!is.na(rasterized_margin_holed[]))] <- 3
rasterized_margin_holed[which(rasterized_margin_holed[]==1)] <- NA
marginasPoints <- rasterToPoints(rasterized_margin_holed)

spdf <- SpatialPointsDataFrame(marginasPoints[,1:2],as.data.frame(marginasPoints))
projection(spdf) <- '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'

writeOGR(obj=spdf, dsn=path_outlines2&"\\margin_icesheet_150m.shp", layer='margin', driver="ESRI Shapefile") # this is in geographical projection

percLand <- 100 * sum(rasterized_margin_holed[] == 2,na.rm=T) / (sum(rasterized_margin_holed[] == 2,na.rm=T) + sum(rasterized_margin_holed[] == 3,na.rm=T))
percLand_icecaps <- 100 * sum(rasterized_margin_holed_promice[] == 2,na.rm=T) / (sum(rasterized_margin_holed_promice[] == 2,na.rm=T) + sum(rasterized_margin_holed_promice[] == 3,na.rm=T))

rbPal <- colorRampPalette(c('red','blue'))
col <- rbPal(2)[as.numeric(cut(marginasPoints[,3],breaks = 2))]
png(file=path_figs&'\\IceMargin_landocean.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
#layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))
plot(coarse_bed,xlab = 'Westing [m]', ylab ='Northing [m]',legend.args=list(text=expression('Elevation [m]')),col=pal)
points(marginasPoints,col=col)
#plot(rasterized_margin_holed,add=TRUE,legend=F,col=rainbow(2,alpha=0.8))
#plot(margin_pEXT_transf,add=T)
grid(NULL,NULL)
dev.off()

col <- rbPal(2)[as.numeric(cut(marginasPoints_promice[,3],breaks = 2))]
png(file=path_figs&'\\IceMargin_landocean_plusicecaps.png', res = 160,width=1800,height=1800)
par(mar=c(5,7,2,4),cex.lab=1.2,cex.axis=1.2)
#layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE))
plot(coarse_bed,xlab = 'Westing [m]', ylab ='Northing [m]',legend.args=list(text=expression('Elevation [m]')),col=pal)
points(marginasPoints_promice,add=T,col=col)
#plot(rasterized_margin_holed,add=TRUE,legend=F,col=rainbow(2,alpha=0.8))
plot(margin_pEXT_transf,add=T)
grid(NULL,NULL)
dev.off()

#levelplot(op2)

#bin <- matrix(0,11,2)
#startK<-0
#for(iExtent in 1:length(seq(extent(op2)[1],extent(op2)[2],150000))){
#  kount <- 150000
#e <- extent(extent(op2)[1] + startK,extent(op2)[1] + startK + kount,-3392600,-632600)
#startK <- startK + kount + 1
#allCells <- extract(op2, e)
#bin[iExtent,] <- cbind(sum(allCells == 2,na.rm=T),sum(allCells == 3,na.rm=T))
#}

#barplot(bin[,1]/rowSums(bin)*100)