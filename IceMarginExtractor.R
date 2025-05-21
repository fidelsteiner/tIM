################################################################################
# Extract Ice Margin Morphology for separate segments for individual analysis
#
# IceMarginExtractor.R
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
##################
# Multiple subsites for analysis of the margin
##################

data_basepath <- 'E:\\GeospatialData'                     # path for heavy input data
data_ArcticDEM <- 'E:\\ArcticDEM'                         # path for ArcticDEM data
output_basepath <- 'C:\\Work\\Research\\tGISM'            # path for code/output   
output_basepath2 <- 'E:\\tGISM'            # path for code/output   

##################
# Paths
##################

path_outlines <- 'E:\\Research\\OtherResearch\\Greenland\\Code\\SectorMargins\\ManualRedRock' # Not used, remove
path_figs <- paste(output_basepath,'\\Figures',sep = "")
pathDEM_Pleiades <- paste(data_basepath,'DEMData\\ValidationDEMs\\RedRock_Pleiades\\Pleiades\\DEM\\2017-08-15_RedRockCliff\\PGO\\2017-08-15_RedRockCliff',sep = "")
pathDEM_Arctic <- data_ArcticDEM
data_Bedmachine <- 'E:\\GeospatialData\\Bedmachine\\5000000451838\\160281892' # path Bedmachine
data_lakes <- 'E:\\GeospatialData\\Lakes\\Greenland' # path Lakes

projec <-'+proj=utm +datum=WGS84'

projec_utm <-'+proj=utm +zone=19 +north +datum=WGS84'
projec_27N <-'+proj=utm +zone=27N +datum=WGS84'

##################
# Loop through DEMs and margin segments for individual processing
##################

# read Bedmachine for later processing
bedM <- raster(paste(data_Bedmachine,'\\BedmachineRaster_proj.tif',sep=''))

# read lakes for later processing
poplak <- ogrInfo(paste(data_lakes,'\\20170101-ESACCI-L3S_GLACIERS-IML-MERGED-fv1_corrected.shp',sep = ""))
margin_lak <- readOGR(dsn=paste(paste(data_lakes,'\\20170101-ESACCI-L3S_GLACIERS-IML-MERGED-fv1_corrected.shp',sep = "")),layer = poplak$layer)
margin_lak <- spTransform(margin_lak, CRS('+proj=longlat +datum=WGS84 +no_defs'))

# read peripheral glaciers for later processing
poppg <- ogrInfo('E:\\GeospatialData\\Periphery\\05_rgi60_GreenlandPeriphery_Dissolved.shp')
margin_pgic <- readOGR(dsn='E:\\GeospatialData\\Periphery\\05_rgi60_GreenlandPeriphery_Dissolved.shp',layer = poppg$layer)
margin_pgic <- spTransform(margin_pgic, CRS('+proj=longlat +datum=WGS84 +no_defs'))

# read complete ice sheet outline for peripheral intersection (2nd outline are all merged subbasins to avoid overlapping margins between peripheral outlines and ice sheet margin being counted twice)
popish <- ogrInfo('E:\\GeospatialData\\cci_data\\outlines\\cci_go_greenland_WGS_fixedintersections.shp')
margin_popish <- readOGR(dsn='E:\\GeospatialData\\cci_data\\outlines\\cci_go_greenland_WGS_fixedintersections.shp',layer = popish$layer)
#margin_popish <- spTransform(margin_popish, CRS('+proj=longlat +datum=WGS84 +no_defs'))
popish2 <- ogrInfo('E:\\GeospatialData\\cci_data\\outlines\\cci_go_greenland_WGS_fixedintersections_mergedsubbasins.shp')
margin_popish2 <- readOGR(dsn='E:\\GeospatialData\\cci_data\\outlines\\cci_go_greenland_WGS_fixedintersections_mergedsubbasins.shp',layer = popish2$layer)
margin_popish2 <- spTransform(margin_popish, CRS('+proj=longlat +datum=WGS84 +no_defs'))
icemargin_total <- spTransform(margin_popish,CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
icemargin_total <- icemargin_total[,3]
names(icemargin_total)[names(icemargin_total) == "area"] <- "Area"
icemargin_total2 <- spTransform(margin_pgic,CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
icemargin_total2 <- icemargin_total2[,9]
margTotal <- rbind(icemargin_total,icemargin_total2)

margSegments <- list.files(paste(output_basepath,'\\Sectors\\SectorBuffers_nuna',sep = ""), pattern='.shp$', all.files=FALSE, full.names=FALSE) # List of segments of the margin with buffer around
margSectors <- list.files(paste(output_basepath,'\\Sectors\\SectorMargins_nuna',sep = ""), pattern='.shp$', all.files=FALSE, full.names=FALSE) # List of segments of the margin without buffer around

DEMSegments <- list.files(paste(data_ArcticDEM,'\\DEMs',sep = ""), pattern='.tif$', all.files=FALSE, full.names=TRUE) # list of all DEMs around the margin of Greenland
DEMSegments_names <- sub('\\.tif$','',list.files(paste(data_ArcticDEM,'\\DEMs',sep = ""), pattern='\\.tif$', all.files=FALSE, full.names=F)) # list of all DEMs around the margin of Greenland

# Prepare slope and roughness maps from base DEMs

for(DEMsl in 1:length(DEMSegments)){

  if(!file.exists(paste('E:\\ArcticDEM\\DEMs\\slopes\\',DEMSegments_names[DEMsl],'_slope.tif',sep = ""))){
  tryCatch({
    DEM_ <- raster(DEMSegments[DEMsl])
    marg <- crop(margtotal,extent(DEM_)+ c(-500,500,-500,500))
    DEM_ <- crop(DEM_,extent(marg) + c(-500,500,-500,500))
    slope_ <- terrain(DEM_, opt="slope", unit="degrees", neighbors=8)
    writeRaster(slope_, filename=paste('E:\\ArcticDEM\\DEMs\\slopes\\',DEMSegments_names[DEMsl],'_slope.tif',sep = ""), format="GTiff", overwrite=TRUE)#}
    
  }, error=function(e){
    tryCatch({
      DEM_ <- raster(DEMSegments[DEMsl])
      marg <- crop(margtotal,extent(DEM_)+ c(-1000,1000,-1000,1000))
      DEM_ <- crop(DEM_,extent(marg) + c(-1000,1000,-1000,1000))
      slope_ <- terrain(DEM_, opt="slope", unit="degrees", neighbors=8)
      writeRaster(slope_, filename=paste('E:\\ArcticDEM\\DEMs\\slopes\\',DEMSegments_names[DEMsl],'_slope.tif',sep = ""), format="GTiff", overwrite=TRUE)#}
      
    }, error=function(e){
      tryCatch({
        DEM_ <- raster(DEMSegments[DEMsl])
        marg <- crop(margtotal,extent(DEM_)+ c(-5000,5000,-5000,5000))
        DEM_ <- crop(DEM_,extent(marg) + c(-1000,1000,-1000,1000))
        slope_ <- terrain(DEM_, opt="slope", unit="degrees", neighbors=8)
        writeRaster(slope_, filename=paste('E:\\ArcticDEM\\DEMs\\slopes\\',DEMSegments_names[DEMsl],'_slope.tif',sep = ""), format="GTiff", overwrite=TRUE)#}
        
      }, error=function(e){
        DEM_ <- raster(DEMSegments[DEMsl])
        slope_ <- terrain(DEM_, opt="slope", unit="degrees", neighbors=8)
        writeRaster(slope_, filename=paste('E:\\ArcticDEM\\DEMs\\slopes\\',DEMSegments_names[DEMsl],'_slope.tif',sep = ""), format="GTiff", overwrite=TRUE)#}
        
      })
    })
  })}

  removeTmpFiles(h=1)

  if(!file.exists(paste('E:\\ArcticDEM\\DEMs\\TRI\\',DEMSegments_names[DEMsl],'_TRI.tif',sep = ""))){
    DEM_ <- raster(DEMSegments[DEMsl])
    slope_ <- raster(paste('E:\\ArcticDEM\\DEMs\\slopes\\',DEMSegments_names[DEMsl],'_slope.tif',sep = ""))
    DEM_ <- crop(DEM_,extent(slope_))
    TRI_ <- terrain(DEM_, opt="TRI") # terrain ruggedness index (Wilson et al.2007)
    writeRaster(TRI_, filename=paste('E:\\ArcticDEM\\DEMs\\TRI\\',DEMSegments_names[DEMsl],'_TRI.tif',sep = ""), format="GTiff", overwrite=TRUE)}
  
  if(!file.exists(paste('E:\\ArcticDEM\\DEMs\\TRI\\',DEMSegments_names[DEMsl],'_TRI.txt',sep = ""))){
    tryCatch({
      #DEM_ <- raster(DEMSegments[DEMsl])
      TRI_ <- raster(paste('E:\\ArcticDEM\\DEMs\\TRI\\',DEMSegments_names[DEMsl],'_TRI.tif',sep = "")) # terrain ruggedness index (Wilson et al.2007)

      marg <- crop(margTotal,extent(TRI_)+ c(-5000,5000,-5000,5000))
      TRI_ <- crop(TRI_,extent(marg)+ c(-10000,10000,-10000,10000))
      #DEM_ <- crop(DEM_,extent(slope_))
      cat(paste0("Masking: ", floor(DEMsl/length(DEMSegments)*100), "%"))
      TRI_ <- mask(TRI_,marg,inverse=T)

      cat(paste0("Masked: ", floor(DEMsl/length(DEMSegments)*100), "%"))
      writeRaster(TRI_, filename=paste('E:\\ArcticDEM\\DEMs\\TRI\\',DEMSegments_names[DEMsl],'_TRI.tif',sep = ""), format="GTiff", overwrite=TRUE)
      tabProgressTRI <- 1
      write.table(tabProgressTRI, file = paste('E:\\ArcticDEM\\DEMs\\TRI\\',DEMSegments_names[DEMsl],'_TRI.txt',sep = ""), sep = "\t",row.names = T)
    }, error=function(e){})}
  
      if(DEMsl %% 1==0) {
    cat(paste0("DEM progress: ", floor(DEMsl/length(DEMSegments)*100), "%"))}
}

# Core Processing

# 301, 302, 303, 304, 305, 306, 307, 308, 315, 313, 312, 311, 310, 309 shear margin Excel

for(i in 1:length(margSegments)){
SubBasinID <- as.numeric(regmatches(margSegments[i], gregexpr("[[:digit:]]+", margSegments[i]))) # original ID of Greenland subbasin
  
if(!file.exists(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'progress_1k.txt',sep = ""))){
  # read in the 100 m buffer along the margin produced from the outline
popOgr <- ogrInfo(paste(output_basepath,'\\Sectors\\SectorBuffers_nuna\\',margSegments[i],sep = ""))
margin_buf <- readOGR(dsn=paste(output_basepath,'\\Sectors\\SectorBuffers_nuna',sep = ""),layer = popOgr$layer)
margin_con <- SpatialPolygons(margin_buf@polygons,proj4string=margin_buf@proj4string)

# read in the actual margin segments
popOgr <- ogrInfo(paste(output_basepath,'\\Sectors\\SectorMargins_nuna\\',margSectors[i],sep = ""))
margin_sec_buf <- readOGR(dsn=paste(output_basepath,'\\Sectors\\SectorMargins_nuna\\',margSectors[i],sep = ""),layer = popOgr$layer)
if(SubBasinID>299){
  margin_sec_con <- margin_sec_buf
  margin_sec_con <- spTransform(margin_sec_con,CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  margin_con <- spTransform(margin_con,CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  margin_sec_con <- st_cast( st_as_sf(margin_sec_con, feature = seq_len(length(terrestrial_margin_sec_con))), "MULTILINESTRING")
  margin_sec_con <- as_Spatial(margin_sec_con)
  
}else{margin_sec_con <- SpatialLines(margin_sec_buf@lines,proj4string=margin_sec_buf@proj4string)}

# produce polylines without nunataks for statistics
marg_nonuna <- st_as_sf(margin_sec_con, feature = seq_len(length(margin_sec_con)))

if(SubBasinID>299){
  library(nngeo)
  marg_nonuna <- st_cast(marg_nonuna, "MULTIPOLYGON")
  marg_nonuna_holes <- marg_nonuna
  marg_nonuna <- st_remove_holes(marg_nonuna)
  marg_nonuna <- st_cast(marg_nonuna, "MULTILINESTRING")
  marg_nonuna <- as_Spatial(marg_nonuna)
  
  marg_nonuna_holes <- st_cast(marg_nonuna_holes, "MULTILINESTRING")
  marg_nonuna_holes <- as_Spatial(marg_nonuna_holes)
  removed_holes <- gDifference(marg_nonuna_holes,marg_nonuna)
  removed_holes <- spTransform(removed_holes, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
}else{
  ## cast to LINESTRING 
  marg_nonuna <- st_cast(marg_nonuna, "LINESTRING")
  marg_nonuna_area <- st_polygonize(marg_nonuna)
  marg_nunas <- which(as.vector(st_area(marg_nonuna_area))>0)
if(!is.null(length(marg_nunas))&length(marg_nunas)>0){marg_nonuna <- as_Spatial(marg_nonuna[-marg_nunas,])}
if(length(marg_nunas)==0){marg_nonuna <- as_Spatial(marg_nonuna)} 
if(is.null(length(marg_nunas))){marg_nonuna <- as_Spatial(marg_nonuna)}}

tab_marginStats <- read.csv(paste(output_basepath,'\\Sectors\\marginstatistics.csv',sep=''))
tab_marginStats[which(tab_marginStats[,1]==SubBasinID),1] <- SubBasinID

# if subbasin part of ice sheet remove intersection with peripheral glaciers
if(SubBasinID<261){
peripheral_subbasin <- crop(margin_pgic,extent(spTransform(margin_con,CRS("+proj=longlat +datum=WGS84 +no_defs")))+c(-2,2,-2,2))
peripheral_subbasin <- rgeos::gBuffer(spTransform(peripheral_subbasin,CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")),width=100)
peripheral_overlap <- gIntersection(peripheral_subbasin, margin_sec_con)
peripheral_overlap_nonuna <- gIntersection(peripheral_subbasin,marg_nonuna)
}

# if subbasin part of peripheral glaciers and ice caps remove intersection with ice sheet
if(SubBasinID>299){
  peripheral_subbasin <- crop(margin_popish2,extent(spTransform(margin_con,CRS("+proj=longlat +datum=WGS84 +no_defs")))+c(-2,2,-2,2))
  peripheral_subbasin <- rgeos::gBuffer(spTransform(peripheral_subbasin,CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")),width=100)
  peripheral_overlap <- gIntersection(peripheral_subbasin, margin_sec_con)
  peripheral_overlap_nonuna <- peripheral_overlap
}

if(is.null(peripheral_overlap)){
  margin_pgiccontact <- 0}
if(!is.null(peripheral_overlap)){
  margin_pgiccontact <- rgeos::gLength(peripheral_overlap)
  peripheral_overlap <- spTransform(peripheral_overlap, CRSobj='+proj=longlat +datum=WGS84 +no_defs')
}

if(is.null(peripheral_overlap_nonuna)){
  margin_pgiccontact_nonuna <- 0}
if(!is.null(peripheral_overlap_nonuna)){
  margin_pgiccontact_nonuna <- rgeos::gLength(peripheral_overlap_nonuna)
  peripheral_overlap_nonuna <- spTransform(peripheral_overlap_nonuna, CRSobj='+proj=longlat +datum=WGS84 +no_defs')
}

# analyse the margin properties with respect to the bed topography
bedM_subbasin <- crop(bedM, extent(spTransform(margin_con,CRS("+proj=longlat +datum=WGS84 +no_defs"))))
bedM_subbasin[bedM_subbasin>10] <- NA
if(SubBasinID==260){
  xy<-xyFromCell(bedM_subbasin,1:length(bedM_subbasin))
  bedM_subbasin[which(xy[,1]<=-37)]<-NA
  bedM_subbasin[which(xy[,2]<=66)]<-NA
}
if(SubBasinID>299){
  bedM_subbasin[bedM_subbasin<=-200] <- NA
  bedM_subbasin <- aggregate(bedM_subbasin,fact=2)
}
if(SubBasinID==301){
  xy<-xyFromCell(bedM_subbasin,1:length(bedM_subbasin))
  bedM_subbasin[which(xy[,1]>=-65&xy[,2]>=76.5)]<-NA
  bedM_subbasin[which(xy[,1]>=-66&xy[,1]<=-60&xy[,2]>=75.7)]<-NA

}

# convert to polygons
if(is.nan(cellStats(bedM_subbasin,'mean',na.rm=T))){
  bedM_subbasin[1] <- 0
  projection(bedM_subbasin) <- '+proj=longlat +datum=WGS84 +no_defs'
}
bedM_polygon <- rasterToPolygons(bedM_subbasin, dissolve=TRUE)
projection(bedM_polygon) <- '+proj=longlat +datum=WGS84 +no_defs'
margin_sec_con <- spTransform(margin_sec_con, CRSobj=projection(bedM_polygon))
marg_nonuna <- spTransform(marg_nonuna, CRSobj=projection(bedM_polygon))

marine_overlap <- gIntersection(margin_sec_con,bedM_polygon)
marine_overlap_nonuna <- gIntersection(marg_nonuna,bedM_polygon)

if(is.null(marine_overlap)){
  margin_marinecontact <- 0}
if(!is.null(marine_overlap)){
  margin_marinecontact <- rgeos::gLength(spTransform(marine_overlap, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")))}

if(is.null(marine_overlap_nonuna)){
  margin_marinecontact_nonuna <- 0}
if(!is.null(marine_overlap_nonuna)){
  margin_marinecontact_nonuna <- rgeos::gLength(spTransform(marine_overlap_nonuna, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")))}

# analyse the margin properties with respect to lakes
lake_overlap <- gIntersection(margin_sec_con, margin_lak)
lake_overlap_nonuna <- gIntersection(marg_nonuna, margin_lak)

if(is.null(lake_overlap)){
  margin_lakcontact <- 0}
if(!is.null(lake_overlap)){
  margin_lakcontact <- rgeos::gLength(spTransform(lake_overlap, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")))
}

if(is.null(lake_overlap_nonuna)){
  margin_lakcontact_nonuna <- 0}
if(!is.null(lake_overlap_nonuna)){
  margin_lakcontact_nonuna <- rgeos::gLength(spTransform(lake_overlap_nonuna, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")))
}

margin_con <- spTransform(margin_con, CRSobj=projection(bedM_polygon))
terrestrial_margin_con <- gDifference(margin_con,bedM_polygon)
terrestrial_margin_con <- gDifference(terrestrial_margin_con,margin_lak)
if(!is.null(peripheral_overlap)){
  terrestrial_margin_con <- gDifference(terrestrial_margin_con,peripheral_overlap)}

terrestrial_margin_sec_con <- gDifference(margin_sec_con,bedM_polygon)
terrestrial_margin_sec_con <- gDifference(terrestrial_margin_sec_con,margin_lak)
if(!is.null(peripheral_overlap)){
  terrestrial_margin_sec_con <- gDifference(terrestrial_margin_sec_con,peripheral_overlap)}

terrestrial_margin_nonuna <- gDifference(marg_nonuna,bedM_polygon)
terrestrial_margin_nonuna <- gDifference(terrestrial_margin_nonuna,margin_lak)
if(!is.null(peripheral_overlap_nonuna)){
  terrestrial_margin_nonuna <- gDifference(terrestrial_margin_nonuna,peripheral_overlap_nonuna)}

# reproject margin to stereo projection for processing with DEMs
terrestrial_margin_sec_con <- spTransform(terrestrial_margin_sec_con, CRSobj="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
terrestrial_margin_con <- spTransform(terrestrial_margin_con, CRSobj="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
terrestrial_margin_nonuna <- spTransform(terrestrial_margin_nonuna, CRSobj="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# stop loop if terrestrial margin is length=0
if(is.null(terrestrial_margin_sec_con)){
  tab_marginStats[which(tab_marginStats[,1]==SubBasinID),2] <- 0 + margin_lakcontact + margin_marinecontact
  tab_marginStats[which(tab_marginStats[,1]==SubBasinID),3] <- 0
  tab_marginStats[which(tab_marginStats[,1]==SubBasinID),7] <- 0 + margin_lakcontact_nonuna + margin_marinecontact_nonuna
  tab_marginStats[which(tab_marginStats[,1]==SubBasinID),8] <- 0
  dir.create(file.path(paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector',sep = "")))
  if(!is.null(marine_overlap)){
    dfM<-SpatialLinesDataFrame(marine_overlap, data.frame(id=1:length(marine_overlap)))
    writeOGR(dfM, dsn=paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'marinemargin.shp',sep = "") ,layer="marine_overlap",driver="ESRI Shapefile")}
  if(!is.null(lake_overlap)){
    dfL<-SpatialLinesDataFrame(lake_overlap, data.frame(id=1:length(lake_overlap)))
    writeOGR(dfL, dsn=paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'lakemargin.shp',sep = "") ,layer="lake_overlap",driver="ESRI Shapefile")}
  if(!is.null(peripheral_overlap)){
    dfP<-SpatialLinesDataFrame(peripheral_overlap, data.frame(id=1:length(peripheral_overlap)))
    writeOGR(dfP, dsn=paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'peripheralconnector.shp',sep = "") ,layer="peripheral_overlap",driver="ESRI Shapefile")}
  
  write.csv(tab_marginStats,paste(output_basepath,'\\Sectors\\marginstatistics.csv',sep=''), row.names = FALSE)
  next
}

# read in margins of manually corrected subbasins
if(SubBasinID %in% c(4,6,7,8,9,11,15,19,20,21,22,26,27,32,34,36,37,38,39,41,42,43,45,46,50,52,53,55,56,57,58,60,61,67,68,70,73,74,76,77,84,86,88,95,96,100,103,104,105,106,108,109,111,113,114,115,117,120,121,126,127,128,129,131,132,134,135,139,142,141,144,147,162,164,165,169,172,180,188,189,190,191,196,199,200,
                     201,202,203,207,208,209,212,213,214,215,237,238,239,243,245,249,250,251,255,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315)){
  terrMargin <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'terrestrialmargin.shp',sep = ""))
  terrMargin <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',sep = ""),layer = terrMargin$layer)
  terrestrial_margin_sec_con <- as(terrMargin,'SpatialLines')
  terrestrial_margin_sec_con <- spTransform(terrestrial_margin_sec_con,CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"))
  terrestrial_margin_sec_con <- gDifference(terrestrial_margin_sec_con,spTransform(margin_lak, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")))
  if(!is.null(peripheral_overlap)){
    if(SubBasinID %in% c(306,307, 308, 309, 310, 311, 141)){}else{
    terrestrial_margin_sec_con <- gDifference(terrestrial_margin_sec_con,spTransform(peripheral_overlap, CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")))}}
  dir.create(file.path(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector',sep = "")))
  dfT<-SpatialLinesDataFrame(terrestrial_margin_sec_con, data.frame(id=1:length(terrestrial_margin_sec_con)))
  
  writeOGR(dfT, dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'terrestrialmargin.shp',sep = "") ,layer="terrestrial_margin_sec_con",driver="ESRI Shapefile")
  

  if(SubBasinID>299){
    terrestrial_margin_nonuna <- gDifference(terrestrial_margin_sec_con,removed_holes)
  }else{
    ## cast to LINESTRING 
    marg_nonuna <- st_as_sf(terrestrial_margin_sec_con, feature = seq_len(length(terrestrial_margin_sec_con)))
    marg_nonuna <- st_cast(marg_nonuna, "LINESTRING")
    marg_nonuna_area <- st_polygonize(marg_nonuna)
    marg_nunas <- which(as.vector(st_area(marg_nonuna_area))>0)
    if(!is.null(length(marg_nunas))&length(marg_nunas)>0){terrestrial_margin_nonuna <- as_Spatial(marg_nonuna[-marg_nunas,])}
    if(length(marg_nunas)==0){terrestrial_margin_nonuna <- as_Spatial(marg_nonuna)} 
    if(is.null(length(marg_nunas))){terrestrial_margin_nonuna <- as_Spatial(marg_nonuna)}}
  if(SubBasinID!=201&!is.null(marine_overlap)){
  marMargin <- ogrInfo(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'marinemargin.shp',sep = ""))
  marMargin <- readOGR(dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',sep = ""),layer = marMargin$layer)
  marMargin <- spTransform(as(marMargin,'SpatialLines'), CRSobj="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
  margin_marinecontact <- rgeos::gLength(marMargin)
  marine_overlap <- marMargin
  
  if(SubBasinID>299){
    marMargin_nonuna <- gDifference(marMargin,removed_holes)
  }else{
    ## cast to LINESTRING 
    marg_nonuna <- st_as_sf(marMargin, feature = seq_len(length(marMargin)))
    marg_nonuna <- st_cast(marg_nonuna, "LINESTRING")
    marg_nonuna_area <- st_polygonize(marg_nonuna)
    marg_nunas <- which(as.vector(st_area(marg_nonuna_area))>0)
    if(!is.null(length(marg_nunas))&length(marg_nunas)>0){marMargin_nonuna <- as_Spatial(marg_nonuna[-marg_nunas,])}
    if(length(marg_nunas)==0){marMargin_nonuna <- as_Spatial(marg_nonuna)} 
    if(is.null(length(marg_nunas))){marMargin_nonuna <- as_Spatial(marg_nonuna)}}}
  margin_marinecontact_nonuna <- rgeos::gLength(marMargin_nonuna)
}else{
  
  dir.create(file.path(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector',sep = "")))
  dfT<-SpatialLinesDataFrame(terrestrial_margin_sec_con, data.frame(id=1:length(terrestrial_margin_sec_con)))
  
  writeOGR(dfT, dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'terrestrialmargin.shp',sep = "") ,layer="terrestrial_margin_sec_con",driver="ESRI Shapefile")
  if(!is.null(marine_overlap)){
    dfM<-SpatialLinesDataFrame(marine_overlap, data.frame(id=1:length(marine_overlap)))
    
    writeOGR(dfM, dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'marinemargin.shp',sep = "") ,layer="marine_overlap",driver="ESRI Shapefile")}
}

if(!is.null(lake_overlap)){
  dfL<-SpatialLinesDataFrame(lake_overlap, data.frame(id=1:length(lake_overlap)))
  writeOGR(dfL, dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'lakemargin.shp',sep = "") ,layer="lake_overlap",driver="ESRI Shapefile")}
if(!is.null(peripheral_overlap)){
  dfP<-SpatialLinesDataFrame(peripheral_overlap, data.frame(id=1:length(peripheral_overlap)))
  writeOGR(dfP, dsn=paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'peripheralconnector.shp',sep = "") ,layer="peripheral_overlap",driver="ESRI Shapefile")}


tab_marginStats[which(tab_marginStats[,1]==SubBasinID),2] <- rgeos::gLength(terrestrial_margin_sec_con) + margin_lakcontact + margin_marinecontact
tab_marginStats[which(tab_marginStats[,1]==SubBasinID),5] <- margin_lakcontact
tab_marginStats[which(tab_marginStats[,1]==SubBasinID),4] <- margin_marinecontact
tab_marginStats[which(tab_marginStats[,1]==SubBasinID),3] <- rgeos::gLength(terrestrial_margin_sec_con)
tab_marginStats[which(tab_marginStats[,1]==SubBasinID),6] <- margin_pgiccontact

tab_marginStats[which(tab_marginStats[,1]==SubBasinID),7] <- rgeos::gLength(terrestrial_margin_nonuna) + margin_lakcontact_nonuna + margin_marinecontact_nonuna
tab_marginStats[which(tab_marginStats[,1]==SubBasinID),10] <- margin_lakcontact_nonuna
tab_marginStats[which(tab_marginStats[,1]==SubBasinID),9] <- margin_marinecontact_nonuna
tab_marginStats[which(tab_marginStats[,1]==SubBasinID),8] <- rgeos::gLength(terrestrial_margin_nonuna)
tab_marginStats[which(tab_marginStats[,1]==SubBasinID),11] <- margin_pgiccontact_nonuna

write.csv(tab_marginStats,paste(output_basepath,'\\Sectors\\marginstatistics.csv',sep=''), row.names = FALSE)

remove(bedM_subbasin)
remove(bedM_polygon)

#######################
# Split up margin into cells for ease of processing
#######################

if(SubBasinID<300){

  # remove small nunataks for slope calculations
marg_nonuna <- st_as_sf(terrestrial_margin_sec_con, feature = seq_len(length(terrestrial_margin_sec_con)))
marg_nonuna <- st_cast(marg_nonuna, "LINESTRING")
marg_nonuna_area <- st_polygonize(marg_nonuna)
marg_nunas <- which(as.vector(st_area(marg_nonuna_area))>0&as.vector(st_area(marg_nonuna_area))<250000)   

if(length(marg_nunas)>0){terrestrial_margin_sec_con <- as_Spatial(marg_nonuna[-marg_nunas,])}
}
                          

if(SubBasinID>299){
  marg_nonuna <- st_as_sf(terrestrial_margin_sec_con, feature = seq_len(length(terrestrial_margin_sec_con)))
  marg_nonuna <- st_cast(marg_nonuna, "MULTILINESTRING")
  marg_nonuna_area <- st_polygonize(marg_nonuna)
  marg_nonuna_area <- st_collection_extract(marg_nonuna_area)
  marg_nunas <- which(as.vector(st_area(marg_nonuna_area))>0&as.vector(st_area(marg_nonuna_area))<250000)   
  nunasToberemoved <- st_cast(marg_nonuna_area[marg_nunas,],'MULTILINESTRING')
  
  #if(length(marg_nunas)>0){terrestrial_margin_con <- gDifference(as_Spatial(marg_nonuna),as_Spatial(nunasToberemoved))}
  
}


bufferedTerr <- gBuffer(terrestrial_margin_sec_con,width=101)
#remove all buffer area that is not close to remaining terrestrial margin
if(SubBasinID==306|308|303){
  terrestrial_margin_con <- gIntersection(spTransform(terrestrial_margin_con,CRSobj="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"),bufferedTerr)
  
}else{
  terrestrial_margin_con <- gIntersection(spTransform(terrestrial_margin_sec_con,CRSobj="+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"),bufferedTerr)
}
grid <- raster(extent(terrestrial_margin_con), resolution = c(1000,1000), crs = proj4string(terrestrial_margin_con))


grid <- raster::extend(grid, c(1,1))
gridPolygon <- rasterToPolygons(grid)

gridPolygon$id <- 1:nrow(gridPolygon)

if(SubBasinID>299){
  #intersectGridClipped <- raster::intersect(gridPolygon, gBuffer(terrestrial_margin_con, byid=TRUE, width=0.1))
  intersectGridClipped <- raster::intersect(gridPolygon, terrestrial_margin_con)
  # merge sections within the same grid cell that were saved separately but with same ID
  n_occur <- data.frame(table(intersectGridClipped$id))
  n_duplicate <- n_occur[n_occur$Freq > 1,]
  if(dim(n_duplicate)[1]>0){
  manip <- intersectGridClipped
  for(dup in 1:dim(n_duplicate)[1]){
    dupls <- which(manip$id==n_duplicate$Var1[dup])
    dupls2 <- which(intersectGridClipped$id==n_duplicate$Var1[dup])
    spatPol <- aggregate(intersectGridClipped[dupls2,])
    spatPol <- SpatialPolygonsDataFrame(spatPol, data.frame(id=1:length(spatPol)))
    spatPol$layer <- NA
    spatPol$id <- as.numeric(as.character(n_duplicate$Var1[dup]))
    manip <- manip[-dupls,]
    manip <- rbind(manip,spatPol)
    dup <- dup+2
  }
  intersectGridClipped <- manip}
  
}else{intersectGridClipped <- raster::intersect(gridPolygon, terrestrial_margin_con)
# merge sections within the same grid cell that were saved separately but with same ID
n_occur <- data.frame(table(intersectGridClipped$id))
n_duplicate <- n_occur[n_occur$Freq > 1,]
if(dim(n_duplicate)[1]>0){
manip <- intersectGridClipped
for(dup in 1:dim(n_duplicate)[1]){
  dupls <- which(manip$id==n_duplicate$Var1[dup])
  dupls2 <- which(intersectGridClipped$id==n_duplicate$Var1[dup])
  spatPol <- aggregate(intersectGridClipped[dupls2,])
  spatPol <- SpatialPolygonsDataFrame(spatPol, data.frame(id=1:length(spatPol)))
  spatPol$layer <- NA
  spatPol$id <- as.numeric(as.character(n_duplicate$Var1[dup]))
  manip <- manip[-dupls,]
  manip <- rbind(manip,spatPol)
  dup <- dup+2
}
intersectGridClipped <- manip}
}

intersectGrid <- gridPolygon[gridPolygon$id %in% intersectGridClipped$id, ]
#dfG<-SpatialPolygonsDataFrame(intersectGrid, data.frame(id=1:length(intersectGrid)))
#dfGC<-SpatialPolygonsDataFrame(intersectGridClipped, data.frame(id=1:length(intersectGridClipped)))
rgdal::writeOGR(intersectGridClipped,paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'GridClipped.shp',sep = ""),layer = 'intersectGridClipped',driver = 'ESRI Shapefile',overwrite_layer = T)
rgdal::writeOGR(intersectGrid,paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'Grid.shp',sep = ""),layer = 'intersectGrid',driver = 'ESRI Shapefile',overwrite_layer = T)

# loop through all DEM outlines to see which ones overlap
if(file.exists(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'ArcticDEMTiles.txt',sep = ""))){
  DEMfils <- read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'ArcticDEMTiles.txt',sep = ""))
  DEMveclist <- match(DEMfils$x,DEMSegments)
}else{
DEMveclist <- vector()
for(k in 1:length(DEMSegments)){
  demdummy <- raster(paste(data_ArcticDEM,'\\DEMs\\',sub('\\.tif$','',basename(DEMSegments))[k],'.tif',sep=''))
    if(gIntersects(terrestrial_margin_con,as(extent(demdummy), "SpatialPolygons"))){
    DEMveclist <- c(DEMveclist, k)
  }
}}

# Extract data from DEMs/slope maps
projecDEMSegments <- list.files(paste(data_ArcticDEM,'\\DEMs',sep=''),pattern='.tif$', all.files=FALSE, full.names=TRUE)
projecslopeSegments <- list.files(paste(data_ArcticDEM,'\\DEMs\\slopes',sep=''),pattern='.tif$', all.files=FALSE, full.names=TRUE)
projecTRISegments <- list.files(paste(data_ArcticDEM,'\\DEMs\\TRI',sep=''),pattern='.tif$', all.files=FALSE, full.names=TRUE)

slopeVals_segment <- vector(mode='list', length=1)
TRIVals_segment <- vector(mode='list', length=1)
lakeVals_segment <- vector(mode='list', length=1)
for(k in 1:dim(intersectGridClipped)[1]){
  slopeVals_segment[[k]] <- NA
  TRIVals_segment[[k]] <- NA
  lakeVals_segment[[k]] <- NA
}

# Export the base ArcticDEM tiles used
write.table(projecDEMSegments[DEMveclist], file = paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'ArcticDEMTiles.txt',sep = ""), sep = "\t",row.names = T)

#}}

for(DEMk in 1 : length(DEMveclist)){
  #DEM_margin <- crop(raster(projecDEMSegments[DEMveclist[DEMk]]),extent(intersectGridClipped))
  #DEM_ <- mask(crop(DEM_margin,crop(terrestrial_margin_con,extent(DEM_margin))),crop(terrestrial_margin_con,extent(DEM_margin)))
  #slope_ <- terrain(DEM_margin, opt="slope", unit="degrees", neighbors=8)
  slope_ <- raster(projecslopeSegments[DEMveclist[DEMk]])

  
  # check which grid cells overlap
  overlapGC<-NA
  for(k in 1:dim(intersectGridClipped)[1]){
    if(!is.null(intersect(ext(intersectGridClipped[k,]),ext(slope_)))){
      overlapGC <- c(overlapGC,k)
    }
  }
  
  #for(k in 1:dim(intersectGridClipped)[1]){
  for(kp in 2:length(overlapGC)){
    k <- overlapGC[kp]
    if(!file.exists(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_',intersectGridClipped$id[k],'_slopeValues.txt',sep = ""))){
    
    if(!is.na(sum(over(spTransform(margin_lak,CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")),intersectGrid[k,])))){
      lakeVals_segment[[k]] <- sum(area(crop(spTransform(margin_lak,CRS("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")),intersectGrid[k,])))
    }else{
      lakeVals_segment[[k]] <- 0}
      tempsegmentMargin <- crop(slope_,extent(intersectGridClipped[k,]))
      segmentMargin <- mask(tempsegmentMargin,intersectGridClipped[k,])
      
      if(file.exists(paste('E:\\ArcticDEM\\DEMs\\TRI\\',DEMSegments_names[DEMveclist[DEMk]],'_TRI.txt',sep = ""))){
        #TRI_ <- raster(paste('E:\\ArcticDEM\\DEMs\\TRI\\',DEMSegments_names[DEMveclist[DEMk]],'_TRI.tif',sep = ""))
        #tempTRIMargin <- crop(TRI_,intersectGrid[k,])
      }
      if(!file.exists(paste('E:\\ArcticDEM\\DEMs\\TRI\\',DEMSegments_names[DEMveclist[DEMk]],'_TRI.txt',sep = ""))){
        #TRI_ <- raster(paste('E:\\ArcticDEM\\DEMs\\TRI\\',DEMSegments_names[DEMveclist[DEMk]],'_TRI.tif',sep = ""))
        #tempTRIMargin <- crop(TRI_,intersectGrid[k,])
        #TRIMarg <- crop(terrestrial_margin_sec_con,extent(intersectGrid[k,])+c(-2000,2000,-2000,2000))
        #tempTRIMargin <- mask(tempTRIMargin,TRIMarg,inverse=T)
      }
      #segmentMargin[segmentMargin[]<5]<-NA
      slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],segmentMargin[!is.na(segmentMargin[])])
      #TRIVals_segment[[k]] <- c(TRIVals_segment[[k]],tempTRIMargin[!is.na(tempTRIMargin[])])
    #} else {
    #  slopeVals_segment[[k]] <- c(slopeVals_segment[[k]],NA)
    #}

    if(length(slopeVals_segment)<dim(intersectGridClipped)[1]){
    slopeVals_segment <- append(slopeVals_segment, NA, after = k)}
  }
    if(k %% 1==0) {
      cat(paste0("subgrid: ", floor(k/dim(intersectGridClipped)[1]*100), "%"))}
  }
  if(DEMk %% 1==0) {
    cat(paste0("subbasin progress: ", floor(DEMk/length(DEMveclist)*100), "%"))}
  removeTmpFiles(h=1)
}

# Save output for this subbasin
tabSB <- matrix(NA, ncol = 8, nrow = length(slopeVals_segment)+1)
tabSB[,1] <- c(intersectGrid$id,'TOTAL')
#tabSB_TRI <- matrix(NA, ncol = 3, nrow = length(TRIVals_segment)+1)
#tabSB_TRI[,1] <- c(intersectGrid$id,'TOTAL')
#tabSB_lake <- matrix(NA, ncol = 2, nrow = length(lakeVals_segment)+1)
#tabSB_lake[,1] <- c(intersectGrid$id,'TOTAL')


for(segk in 1:length(slopeVals_segment)){
  if(!file.exists(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_',intersectGridClipped$id[segk],'_slopeValues.txt',sep = ""))){
    
  write(slopeVals_segment[[segk]], paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_',intersectGridClipped$id[segk],'_slopeValues.txt',sep = ""))}else{
  slopeVals_segment[segk]<-read.table(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_',intersectGridClipped$id[segk],'_slopeValues.txt',sep = ""),fill = TRUE)}

  tabSB[segk,2:6] <- quantile(unlist(slopeVals_segment[segk]), probs = c(0.05, 0.25, 0.5, 0.75, 0.95),na.rm=T)
tryCatch({
  tabSB[segk,7] <- gLength(crop(terrestrial_margin_sec_con,intersectGrid[segk,]))
}, error=function(e){
  tabSB[segk,7] <- 0
  })
tabSB[segk,8] <- length(which(!is.na(slopeVals_segment[[segk]])))
#tabSB_TRI[segk,2] <- mean(unlist(TRIVals_segment[[segk]]),na.rm=T)
#tabSB_lake[segk,2] <- unlist(lakeVals_segment[[segk]])
#tabSB_TRI[segk,3] <- length(which(!is.na(TRIVals_segment[[segk]])))
#}else{
#slopeVals_segment[[segk]] <- table.read(paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_',intersectGridClipped$id[segk],'_slopeValues.txt',sep = ""))
#tabSB[segk,8] <- length(which(!is.na(slopeVals_segment[[segk]])))
#}
}
finStat <- length(slopeVals_segment)+1
tabSB[finStat,2:6] <- quantile(unlist(slopeVals_segment), probs = c(0.05, 0.25, 0.5, 0.75, 0.95),na.rm=T)
tabSB[finStat,8] <- length(which(!is.na(unlist(slopeVals_segment))))
#tabSB_TRI[finStat,2] <- mean(unlist(TRIVals_segment),na.rm=T)
#tabSB_TRI[finStat,3] <- length(which(!is.na(unlist(TRIVals_segment))))
#tabSB_lake[finStat,2] <- sum(unlist(lakeVals_segment),na.rm=T)

colnames(tabSB) <- c('grid ID',0.05, 0.25, 0.5, 0.75, 0.95,'length','n')
write.csv(tabSB,paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_slopePercentiles.csv',sep = ""), row.names = FALSE)

#colnames(tabSB_TRI) <- c('grid ID','TRI','n')
#write.csv(tabSB_TRI,paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_TRI.csv',sep = ""), row.names = FALSE)

#colnames(tabSB_lake) <- c('grid ID','area')
#write.csv(tabSB_lake,paste(output_basepath,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'_lakeareas.csv',sep = ""), row.names = FALSE)

# just to track progress
tabProgress <- 1
write.table(tabProgress, file = paste(output_basepath2,'\\Sectors\\',SubBasinID,'Sector\\',SubBasinID,'progress_1k.txt',sep = ""), sep = "\t",row.names = T)


if(i %% 1==0) {
  cat(paste0("greenland wide progress: ", floor(i/length(margSegments)*100), "%"))}
}
}