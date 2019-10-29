collapse.k <- function(scrFrame){
  
  scrFrame$caphist <- sapply(scrFrame$caphist,
                             function(x) apply(x,c(1,2),sum))
  for(i in 1:length(scrFrame$caphist)){
    ch <- scrFrame$caphist[[i]]
    tmp <- array(ch, c(nrow(ch),ncol(ch),1))
    scrFrame$caphist[[i]] <- tmp 
  }

  if(!is.null(scrFrame$trapCovs)){
    message("\nImportant: Time varying trap covaraites are no longer useful.\n")
    varocc <- all(scrFrame$occasions == scrFrame$occasions[1])
    vartrap <- numeric(length(scrFrame$trapCovs))
    for(i in 1:length(scrFrame$trapCovs)){
      scrFrame$trapCovs[[i]] <- list(scrFrame$trapCovs[[i]][[1]])
      eff <- apply(scrFrame$trapOperation[[i]],1,sum)
      scrFrame$trapCovs[[i]][[1]]$effort <- eff
      vartrap[i] <- all(eff == eff[1])
    }
    if(!varocc | any(vartrap == FALSE)){
      message("Important: Traps have uequal effort, use 'offset(log(effort))' is detection model.\n")
    }
  }else{
    scrFrame$trapCovs <- list()
    for(i in 1:length(scrFrame$caphist)){
      eff <- apply(scrFrame$trapOperation[[i]],1,sum)
      effdf <- data.frame(effort = eff)
      scrFrame$trapCovs[[i]] <- list(effdf)
    }
  }
  
  for(i in 1:length(scrFrame$trapOperation)){
    scrFrame$trapOperation[[i]] <- scrFrame$trapOperation[[i]][,1,drop=FALSE]
    scrFrame$trapOperation[[i]][] <- 1
  }
  
  scrFrame$occasions[] <- 1
  return(scrFrame)
}
