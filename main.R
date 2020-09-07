### REMOTE SENSING INPUTS DOWNLOADER USING RGEE

# preliminaries... assuming you already installed rgee
rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rgdal,sf,raster,cptcity,rgee,foreach,doParallel)
ee_Initialize(email='arnanaraza2006@gmail.com',drive = TRUE)

# global variables
main.dir <- 'D:/RGEE/'
setwd(main.dir)

## Funciton to download satellite data within polygons from csv/shp 
DL <- function(year, aoi, satellite, outdir, rsl){
  
  # Year filter
  year2 <- as.character(year+1)
  start <- ee$Date(paste0(year,"-01-01"))
  end <- ee$Date(paste0(year2,"-01-01"))
  
  # Use block polys Define a region of interest with sf
  aoi_sf <- st_as_sf(aoi)
  st_crs(aoi_sf) <-  4326
  ee_roi <- aoi_sf %>%
    st_geometry() %>%
    sf_as_ee()
  
  # Get ALOS-PALSAR data and take the mean of image collection
  if (satellite == 'ALOS'){
    sat.col <-ee$ImageCollection("JAXA/ALOS/PALSAR/YEARLY/SAR") #2011-2015 missing 
    filter <- sat.col$filterBounds(ee_roi)$filterDate(start,end)
    b <- c("HH", "HV")
    img <- filter$select(b)$mean()
  }
  
  # Get Landsat data including overlapping years between two Landsat satellites
  else if (satellite == 'LANDSAT'){
    
    # Filter out poor quality pixels
    getQABits <- function(image, qa) {
      qa <- sum(2^(which(rev(unlist(strsplit(as.character(qa), "")) == 1))-1))
      image$bitwiseAnd(qa)$lt(1)}
    bnames <- c('swir2', 'swir1', 'nir', 'red', 'pixel_qa')
    
    if (year >= 2013){
      b <- c("B7", "B5", "B5", "B4","pixel_qa")
      b1 <- c("B7", "B5", "B4", "B3",'pixel_qa')
      sat.col0 <-ee$ImageCollection("LANDSAT/LC08/C01/T1_SR")$select(b,bnames) 
      sat.col1 <- ee$ImageCollection("LANDSAT/LE07/C01/T1_SR")$select(b1,bnames)  
      sat.col <- sat.col0$merge(sat.col1)}
    
    else if (year < 2013){
      b <- c("B7", "B5", "B4", "B3",'pixel_qa')
      sat.col1 <- ee$ImageCollection("LANDSAT/LE07/C01/T1_SR")$select(b,bnames)  
      sat.col0 <-ee$ImageCollection("LANDSAT/LT05/C01/T1_SR")$select(b,bnames)
      sat.col <- sat.col0$merge(sat.col1)
    } else{print('wrong year')}
    
    # clean the bands and calculate VIs
    clean_bnds <- function(img) {
      bnds <- img$select(bnames[-5])
      qa <- img$select("pixel_qa")
      quality_mask <- getQABits(qa, "00000100001")
      bnds %>%
        ee$Image$updateMask(quality_mask) %>%
        ee$Image$copyProperties(img, list("system:time_start"))}
    
    clean_ndvi <- function(img) {
      ndvi_values <- img$normalizedDifference(bnames[3:4])
      qa <- img$select("pixel_qa")
      quality_mask <- getQABits(qa, "00000100001")
      ndvi_values %>%
        ee$Image$updateMask(quality_mask) %>%
        ee$Image$copyProperties(img, list("system:time_start"))}
    
    clean_ndwi <- function(img) {
      ndwi_values <- img$normalizedDifference(bnames[c(3,2)])
      qa <- img$select("pixel_qa")
      quality_mask <- getQABits(qa, "00000100001")
      ndwi_values %>%
        ee$Image$updateMask(quality_mask) %>%
        ee$Image$copyProperties(img, list("system:time_start"))}
    
    clean_ndwi1 <- function(img) {
      ndwi_values1 <- img$normalizedDifference(bnames[c(3,1)])
      qa <- img$select("pixel_qa")
      quality_mask <- getQABits(qa, "00000100001")
      ndwi_values1 %>%
        ee$Image$updateMask(quality_mask) %>%
        ee$Image$copyProperties(img, list("system:time_start"))}
    
    # Create a yearly composite
    ndvi <- sat.col$filterBounds(ee_roi)$filterDate(start,end)$map(clean_ndvi)$max()
    ndwi1 <- sat.col$filterBounds(ee_roi)$filterDate(start,end)$map(clean_ndwi1)$max()
    img <- ee$Image(ndvi)$addBands(ndwi1)
  }
  
  
  # Polygon/Block-level download
  dir.create(file.path(main.dir, outdir))
  outdir <- paste0(main.dir,outdir)
  
  r_list <- list()
  
  for (i in 1:length(aoi)){
    dsn <-paste0(outdir,  satellite, '_', year,'.tif')
    print(dsn)
    a <- aoi[i,]@bbox
    g <- c(a[1],a[2],a[3],a[4])
    geometry <- ee$Geometry$Rectangle(
      coords = g,
      proj = "EPSG:4326",
      geodesic = FALSE)
    
    ee_raster <- ee_as_raster(
      image = img,
      region = geometry,
      dsn = dsn,
      scale = rsl,
      via = "getInfo",
      maxPixels=1e+100)
    r_list <- c(r_list, ee_raster)
  }
  return(r_list)
}

#download RS inputs in raster blocks 
landsat_yrs <- c(2000:2018)
plt <- read.csv(paste0(main.dir,'data/sample_sites.csv'))
coordinates(plt) <- ~long+lat

source(paste0(main.dir,'scripts/MakeBlockPolygon.R')) #point to square polygon function
aoi <- MakeBlockPolygon(plt, 0.1, 1)
aoi1 <- lapply(landsat_yrs, function(x) DL(x, aoi[1,], 'LANDSAT', 
                                        paste0('results/PH_',plt$site[1],'_100m/'),100))  #resampled to 100m (faster demo)

#loop the aoi polygons for ndvi time series 
for (i in 1:nrow(aoi)){
  lapply(landsat_yrs, function(x) DL(x, aoi[i,], 'LANDSAT', 
                                          paste0('results/PH_',plt$site[i],'_100m/'),100))
}

