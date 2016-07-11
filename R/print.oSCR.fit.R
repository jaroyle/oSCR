print.oSCR.fit <- function(x, burn=NULL, ...){
  if("parameters" %in% names(x$outStats)){
    if("model" %in% names(x$call)){
      mod <- x$call[["model"]]
    }else{
      mod <- list(D~1,p0~1,sig~1,asu~1)
    }
    tmpFit <- cbind(x$outStats[,3],
                    x$outStats[,4],
                    x$outStats[,3]/x$outStats[,4],
                    2*(1 - pnorm(abs(x$outStats[,3]/x$outStats[,4]))))
    rownames(tmpFit) <- x$outStats[,1]
    colnames(tmpFit) <-c("Estimate","SE","z","P(>|z|)")
    cat(" Model: ", paste(mod)[-1],fill=TRUE)
    cat(" Run time: ", x$proctime," minutes",fill=TRUE)
    cat(" AIC: ", x$AIC,fill=TRUE)
    cat(" ",fill=TRUE)
    cat("Summary table:","\n")
    print(round(tmpFit,3))
    cat("*Density is per pixel density")
    
  }else{
    if("model" %in% names(x$call)){
      mod <- x$call[["model"]]
    }else{
      mod <- list(D~1,p0~1,sig~1,asu~1)
    }
    tmpFit <- x$outStats
    cat("Model: ", paste(mod)[-1],fill=TRUE)
    cat("Run time: ", x$proctime," minutes",fill=TRUE)
    cat("AIC: ", x$AIC,fill=TRUE)
    cat(" ",fill=TRUE)
    cat("Summary table:","\n")
    print(round(tmpFit,3))
    cat("*Density is per pixel density")
  }
}
