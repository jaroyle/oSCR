print.scrFrame <- function(scrFrame){
  
  caphist.dimensions <- sapply(scrFrame$caphist,dim)
  
  if(nrow(caphist.dimensions)==2)
    caphist.dimensions <- rbind(caphist.dimensions,1)
  
  ave.caps <- numeric(ncol(caphist.dimensions)) 
  ave.spat.caps <- numeric(ncol(caphist.dimensions))
  mmdm.session <- numeric(ncol(caphist.dimensions))
  for(i in 1:length(scrFrame$caphist)){
    max.dist <- NULL
    ave.caps[i] <- mean(apply(apply(scrFrame$caphist[[i]],c(1,3),sum),1,sum))
    ave.spat.caps[i] <- mean(apply(apply(scrFrame$caphist[[i]],c(1,2),sum)>0,1,sum))
    for (j in 1:nrow(scrFrame$caphist[[i]])) {
      if(dim(scrFrame$caphist[[i]])[3]>1){
        where <- apply(scrFrame$caphist[[i]][j, , ], 1, sum) > 0
      }else{
        where <- scrFrame$caphist[[i]][j, , ] > 0
      }
      if (sum(where) > 1)
        max.dist <- c(max.dist, max(0, dist(scrFrame$traps[[i]][where, c("X", "Y")]), na.rm = T))
    }
    mmdm.session[i] <- mean(max.dist[max.dist > 0], na.rm = T)
  }
  caphist.dimensions <- rbind(caphist.dimensions,
                              ave.caps,
                              ave.spat.caps,
                              mmdm.session)
  
  dn <- list(c("n individuals", "n traps", "n occasions", "avg caps", "avg spatial caps", "mmdm"),
             paste("S",1:ncol(caphist.dimensions),sep=""))
    
  dimnames(caphist.dimensions) <- dn
  
  if(length(scrFrame$caphist)>1){
    cd <- as.data.frame(caphist.dimensions)
    cat("",fill=TRUE)
    print(round(cd[1:3,]))
    cat("",fill=TRUE)
    print(round(cd[4:6,],2))
  
    cat("",fill=TRUE)
    cat("Pooled MMDM: ",round(scrFrame$mmdm,2),fill=TRUE)
  }else{
    cd <- as.data.frame(caphist.dimensions)
    cat("",fill=TRUE)
    print(round(t(t(cd[1:3,,drop=FALSE]))))
    cat("",fill=TRUE)
    print(round(cd[4:6,,drop=FALSE],2))
  }
}
