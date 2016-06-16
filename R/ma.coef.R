ma.coef <- function(ms){

  ests <- ms$coef.tab[,-1]
  ests[is.na(ests)] <- 0
  wts <- ms$aic.tab$weight
  vi <- t(ests!=0) %*% wts 
  
  ma.beta.shrink <- colSums(ests * wts)
  tmp <- matrix(ma.beta.shrink,dim(ests)[1],dim(ests)[2],byrow=T)
  if(ms$se){
    vars <- ms$se.tab[,-1]
    vars[is.na(vars)] <- 0
    ma.se.shrink <- colSums(wts*sqrt(vars^2 + (ests-tmp)^2))  
  }else{
    ma.se.shrink <- NA
  }

  ma.beta <- numeric(ncol(ests))
  ma.se <- numeric(ncol(ests))

  for(i in 1:ncol(ests)){
    rmv <- which(ests[,i] !=0)
    rescale.wt <- wts[rmv]/sum(wts[rmv])  
    ma.beta[i] <- sum(ests[rmv,i] * rescale.wt)
    if(ms$se){
      ma.se[i] <- sum(rescale.wt * sqrt(vars[rmv,i]^2 + (ests[rmv,i]-ma.beta[i])^2))
    }else{
      ma.se[i] <- NA
    }
  }
  ma.coef <- data.frame(colnames(ests),
                        cbind(ma.beta,
                              ma.se,
                              ma.beta.shrink,
                              ma.se.shrink,
                              vi))
  rownames(ma.coef) <- NULL
  colnames(ma.coef) <- c("Parameter", "Estimate", paste("Std. Error",sep=""), 
                         "Estimate*", paste("Std. Error*",sep=""), "RVI")
  class(ma.coef) <- c("ma.coef")
  return(ma.coef)
}

