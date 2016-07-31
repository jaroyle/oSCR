make.scrFrame <- function(caphist, traps, indCovs=NULL, 
                          trapCovs=NULL, trapOperation=NULL, type="scr"){
  
  #must have caphist and traps
  if(any(is.null(caphist),is.null(traps)))
    stop("caphist and trap must be provided")
  
  #caphist
  if(!is.list(caphist))
    stop("caphist must be a list")
  n.sessions <- length(caphist)  
  caphist.dimensions <- sapply(caphist,dim)
  
  if(nrow(caphist.dimensions)==2)
    caphist.dimensions <- rbind(caphist.dimensions,1)
  
  for(i in 1:n.sessions){
    caphist[[i]] <- array(caphist[[i]], dim=caphist.dimensions[,i])
    all.zero <- apply(apply(caphist[[i]],c(1,3),sum),1,sum)
    if(any(all.zero==0))
      stop("At least one individual has an all-zero encounter history")
  }
  
  #indCovs
  if(!is.null(indCovs)){
    if(!is.list(indCovs))
      stop("indCovs must be a list")
    if(any(!sapply(indCovs,is.data.frame)))
      stop("indCovs must be a list of dataframes")
    if(length(indCovs) != length(caphist))
      stop("number of sessions in indCovs does not match capphist")
    
    check.dim <- sapply(indCovs,nrow) 
    if(any(check.dim!=caphist.dimensions[1,]))
      stop("number of individuals in indCovs does not match caphist")
    if(!("rmv" %in% indCovs[[1]])){
      for(i in 1:length(indCovs)){
        indCovs[[i]]$removed <- dim(caphist[[i]])[3]
      }
    }
  }else{
    indCovs <- list()
    for(i in 1:length(caphist)){
      indCovs[[i]] <- data.frame(removed = rep(dim(caphist[[i]])[3],dim(caphist[[i]])[1]))
    }
  }
  
  #traps
  if(!is.list(traps))
    stop("traps must be a list")
  #if(any(!sapply(traps,is.data.frame)))
  #  stop("traps must be a list of dataframes")
  if(length(traps)!=length(caphist))
    stop("number of sessions in traps does not match caphist")
  
  check.dim <- sapply(traps,nrow) 
  if(!all(check.dim==caphist.dimensions[2,]))
    stop("number of traps does not match caphist")
  
  #trapCovs
  if(!is.null(trapCovs)){
    if(!is.list(trapCovs))
      stop("trapCovs must be a list")
    if(any(!sapply(trapCovs,is.list)))
      stop("trapCovs must be a list of lists")
    if(any(!unlist(sapply(trapCovs,function(x)sapply(x,is.data.frame)))))
      stop("trapCovs must be a list of dataframes")
    if(length(trapCovs) != length(caphist))
      stop("number of sessions in trapCovs does not match capphist")
    #check.dim <- sapply(trapCovs,function(x)sapply(x,nrow))
    check.dim <- lapply(trapCovs,function(x)sapply(x,nrow))
    for(i in 1:length(check.dim)){
      if(!all(check.dim[[i]]==caphist.dimensions[2,i]))
        stop("number of traps does not match caphist")
    }
  }
  
  #trapOperation
  if(!is.null(trapOperation)){
    if(!is.list(trapOperation))
      stop("trapOperation must be a list")
    #if(any(!sapply(trapOperation,is.data.frame)))
    #  stop("trapOperation must be a list of dataframes")
    if(length(trapOperation) != length(caphist))
      stop("number of sessions in trapOperation does not match capphist")
    check.dim <- sapply(trapOperation,nrow) 
    if(!all(check.dim==caphist.dimensions[2,]))
      stop("number of traps does not match caphist")
  }
  
  #mean maximum distance moved
  max.dist <- NULL
  for (i in 1:length(caphist)) {
    for (j in 1:nrow(caphist[[i]])){
      if(dim(caphist[[i]])[3]>1){
        where <- apply(caphist[[i]][j, , ], 1, sum) > 0
      }else{
        where <- caphist[[i]][j, , ] > 0
      }
      if (sum(where) > 1)
        max.dist <- c(max.dist, max(0, dist(traps[[i]][where, c("X", "Y")]), na.rm = T))
    }
  }
  mmdm <- mean(max.dist[max.dist > 0], na.rm = T)
  
  
  scrFrame <- list("caphist" = caphist,
                   "traps" = traps,
                   "indCovs" = indCovs,
                   "trapCovs" = trapCovs,
                   "trapOperation" = trapOperation,
                   "occasions" = caphist.dimensions[3,],
                   "type" = type,
                   "mmdm" = mmdm)
  
  class(scrFrame) <- "scrFrame"  
  return(scrFrame)
}
