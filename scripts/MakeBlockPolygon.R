# make Polygon containing x, y and aligning to AGB map pixel  
# or larger cell over which AGB is aggregated
SRS <- CRS('+init=epsg:4326')

Block <- function(x, y, size){
  xll <- size * x %/% size
  yll <- size * y %/% size
  pol0 <- Polygon(cbind(c(xll, xll+size, xll+size, xll, xll),
                        c(yll, yll, yll+size, yll+size, yll)))
  pol1 <- Polygons(list(pol0), "pol")
  return(SpatialPolygons(list(pol1), proj4string=SRS))
}

MakeBlockPolygon <- function(plots, aggr = NULL,  minPlots = 1){
           
  # aggregate if aggr != NULL
  if(!is.null(aggr)){
    # aggregate to aggr degree cells
    plots$Xnew <- aggr * (0.5 + plots@coords[,1] %/% aggr)                   
    plots$Ynew <- aggr * (0.5 + plots@coords[,2] %/% aggr)                   
    
    plotsTMP <- aggregate(plots$Xnew, list(plots$Xnew, plots$Ynew), 
                          mean, na.rm=T)
    names(plotsTMP) <- c("POINT_X","POINT_Y", 'x')
    
    # only keep plots satisfying minPlots criterion
    if(minPlots > 1){
      blockCOUNT <- aggregate(plots@coords[,1], list(plots$Xnew, plots$Ynew), 
                              function(x) length(na.omit(x)))
      ndx <- which(blockCOUNT$x >= minPlots)
      plotsTMP <- plotsTMP[ndx,]
    }
    plots <- plotsTMP
    rsl <- aggr
    print(nrow(plots))
  } 
  
  # sample forest fraction and AGB data per cell/plot
  nc <- detectCores()
  cl <- makeCluster(nc-1)
  registerDoParallel(cl, nc)
  
  FFAGB <- foreach(i=1:nrow(plots), .combine='rbind', .errorhandling = 'remove',
                   .packages='raster', .export=c('Block', 'SRS')) %dopar% {
                                                   pol <- Block(plots$POINT_X[i],plots$POINT_Y[i], rsl)
                                                   id <- data.frame(ID=1:length(pol), row.names='pol')
                                                   SpatialPolygonsDataFrame(pol,id)
                                                 }
  stopCluster(cl)
  return(FFAGB)
}


#aus_plt1 <- invDasymetry(aus_plt, 0.25, 5)
#head(aus_plt1)
#aus_plt2 <- setdiff(as.data.frame(aus_plt1[c('x','y')]), as.data.frame(aoi[c('x', 'y')]))
#id <- data.frame()
#aus_plt3 <- SpatialPolygonsDataFrame(aus_plt2, id)