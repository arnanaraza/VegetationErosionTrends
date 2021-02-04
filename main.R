### GEE DATA DOWNLOADER USING RGEE

# preliminaries... assuming you already installed rgee
rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rgdal,sf,raster,sf,plyr,dplyr,cptcity,bfast,devtools,remotes,zoo,
               rgee,foreach,doParallel,rlist,gdalUtils,ranger)
ee_Initialize(email='arnanaraza2006@gmail.com') ####!!!!!!!!

# global variables
main.dir <- 'C:/VegetationErosionTrends/'
setwd(main.dir)

## Function to download satellite data in areas of interest (aoi)
DL <- function(year=2007, aoi=aoi, satellite='LANDSAT', outdir='results/PH_dam1_100m', rsl=30){
  
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
    bnames <- c('swir2', 'swir1', 'nir', 'red', 'pixel_qa') #bands of interest
    
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
    img <- sat.col$filterBounds(ee_roi)$filterDate(start,end)$map(clean_ndvi)$max()
    ndwi1 <- sat.col$filterBounds(ee_roi)$filterDate(start,end)$map(clean_ndwi1)$max()
  #  img <- ee$Image(ndvi)$addBands(ndwi1) #did not include ndwi for faster demo
  }
  
  
  
  # Block-level download
  dir.create(file.path(main.dir, outdir))
  outdir <- paste0(main.dir,outdir)
  
  r.list <- list()
  
  lf <- list.files(outdir,satellite)
  omit <- c('aux', paste0(rsl, 'm')) #omit aux and vrt file
  t <- unique(grep(paste(omit,collapse="|"), lf, value=T))
  l <- setdiff(lf,t)
  l <- length(l)#for interrupted downloads
  if (nrow(aoi) == l){ 
    idx <- l}else{idx <- l + 1}
  
  for (i in idx:length(aoi)){
    dsn <-paste0(outdir,  satellite, '_', i)
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
    r.list <- c(r.list, ee_raster)
  }
  return(r.list)
}

##################################################

#download RS inputs in blocks 
landsat_yrs <- c(2000:2014)
plt <- read.csv(paste0(main.dir,'data/td.03.ha.csv')) #assuming point data (and not polygons) are available
coordinates(plt) <- ~long+lat

source(paste0(main.dir,'scripts/MakeBlockPolygon.R')) #point to square polygon function
aoi <- MakeBlockPolygon(plt, 0.01, 2) #10x10km blocks
aoi1 <- lapply(landsat_yrs, function(x) DL(x, aoi, 'LANDSAT', 
                                        paste0('results/PH_NFI_',x,'_100m/'),100))  #resampled to 100m (faster demo)


TS_Values <- function(plt=plt,country='AFR',yrs=landsat_yrs,satellite='LANDSAT',rsl=100,aoi=aoi){
  
  # mosaic and make vrt
  TS <- function(country, yr, satellite){
    folder_dir <- paste0(country,'_',yr,'_',rsl,'m')
    setwd(paste0(main.dir, 'results/',folder_dir))
    r_files <- list.files(getwd(), satellite)
    r_files <- r_files[!grepl('vrt|aux', r_files)]
    fname <- paste0(folder_dir,'_',satellite,".vrt")
    gdalbuildvrt(gdalfile = r_files, # uses all tiffs in the current folder
                 output.vrt = fname, overwrite=T)
    brick(fname)
  }
  
  # raster brick of all block data at regional scale
  mosaic_l <- lapply(landsat_yrs, function(x) TS(country,x, satellite))
  if(satellite == 'LANDSAT'){
    bnames <- c('ndvi')
  }
  mosaic_l <- lapply(mosaic_l, setNames, nm = bnames)   #stack of aoi yearly mosaic!
  mosaic_l <- lapply(mosaic_l, function(x) crop(x, aoi))
  mosaic_l <- lapply(mosaic_l, function(x) mask(x, aoi))
  
  #plot-level valuetable assembly
  val <- lapply(1:length(mosaic_l),  #should extract all years for all points!
                function(x) extract(mosaic_l[[x]], plt))                             
  yrs_rep <- do.call(rbind, replicate(nrow(plt), coredata(as.data.frame(landsat_yrs)), 
                                      simplify = FALSE))
  ts_rep <- do.call(rbind, replicate(length(mosaic_l), 
                                     coredata(plt), simplify = FALSE))
  
  #add temporal covs of plot time series data
  YearlyStack <- function(r_list,b){
    bnd <- lapply(r_list, function(x) stack(x[[b]]))
    bnd1 <- do.call(stack, bnd)
    
    nc <- detectCores()
    cl <- makeCluster(nc-1)
    registerDoParallel(cl, nc)
    
    stat <- foreach(i=1:nrow(aoi), .combine='merge', .errorhandling = 'remove',
                    .packages='raster', .export='aoi') %dopar% {
                      bnd_s <- crop(bnd1, aoi[i,]) ###!!!
                      time <- 1:nlayers(bnd_s)
                      lmod <- function(x) { lm(x ~ time)$coefficients[2] }
                      stack(calc(bnd_s, fun = mean, na.rm = T),calc(bnd_s, fun = min, na.rm = T),
                            calc(bnd_s, fun = max, na.rm = T),calc(bnd_s, fun = sd, na.rm = T),
                            calc(bnd_s, fun = lmod)) ### too much NA pixels!
                    }
    stopCluster(cl)
    names(stat)
    if (b==5){names(stat) <- c('nbr_mean', 'nbr_min', 'nbr_max', 'nbr_sd','nbr_slp')}
    if(b==6){names(stat) <- c('ndvi_mean', 'ndvi_min', 'ndvi_max', 'ndvi_sd', 'ndvi_slp')}
    if(b==7){names(stat) <- c('nbmi_mean', 'ndmi_min', 'ndmi_max', 'ndmi_sd', 'ndmi_slp')}
    
    return(stat)
  }
  start.time <- Sys.time()
  stats_l <- lapply(1, function(x) YearlyStack(mosaic_l,x)) #stats of ndvi 
  end.time <- Sys.time()
  print(end.time - start.time)

  ndvi_temp <- do.call(rbind, replicate(length(landsat_yrs), coredata(extract(stats_l[[1]], plt)), 
                                        simplify = FALSE))
  #add texture ndvi
  r=mosaic_l[[4]]$ndvi
  tex <- glcm(r, window = c(3,3), shift=c(1,1),
             statistics = c("mean", "variance", "homogeneity", "contrast", "entropy"))
  
  #make one data frame
  vt <- data.frame(ts_rep, yrs_rep,ldply(val,data.frame), ndvi_temp,
                   as.data.frame(extract(tex, plt)))
  covs <- lapply(mosaic_l, function(x) stack(x,stats_l[[1]]))
  bnames <- c("ndvi","ndvi_mean", "ndvi_min", "ndvi_max","ndvi_sd" , "ndvi_slp")
  names(vt)[c(17:22)] <- bnames
  head(vt)
  covs <- lapply(covs, setNames, nm = bnames)  
  list(covs, vt)
}

vtCovs <- TS_Values(plt,'PH_NFI',landsat_yrs,'LANDSAT',100, aoi)

covs <- vtCovs[[1]]
vt <- vtCovs[[2]]
y <-as.numeric(vt [,menu(names(vt), title="which column is the latitude")]) #11 10 56
x <-as.numeric(vt [,menu(names(vt), title="which column is the longitude")]) 
agb <-as.numeric(vt [,menu(names(vt),  title="which column is your AGB")]) 

vt <- data.frame(agb,vt[,c('ndvi','glcm_mean' ,'glcm_variance',
                           'glcm_homogeneity', 'glcm_contrast' ,'glcm_entropy' )])
vt <- na.omit(vt)
head(vt)
vt1 <- vt[1:nrow(plt),]
rf.mod <- ranger(vt1$agb ~ ., data=vt1,
                 importance='permutation', save.memory = T)
rf.mod
mean(rf.mod$predictions)
importance(rf)
