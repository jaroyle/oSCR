extract.rast <- function(ss, rast, mult=1, cov.name="val.1",func=median){
  for(i in 1:length(ss)){
    tmpS <- ss[[i]][,c("X","Y")]*mult
    id <- factor(1:nrow(tmpS))
    tmpR <- cbind(tmpS,id)
    aa <- rasterFromXYZ(tmpR, crs=projection(rast))
    bb <- rasterToPolygons(aa)
    r1 <- extract(rast,bb,fun=func)
    ss[[i]][,cov.name] <- r1
  }
  return(ss)
}
