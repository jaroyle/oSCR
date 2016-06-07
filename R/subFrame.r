 subFrame <- function(scrFrame,subs){
   
   caphist <- scrFrame$caphist[subs]
   traps <- scrFrame$traps
   if(!is.null(scrFrame$indCovs)){
     indCovs <- scrFrame$indCovs[subs]
   }else{
     indCovs <- scrFrame$indCovs[subs]
   }
   if(!is.null(scrFrame$trapCovs)){
     trapCovs <- scrFrame$trapCovs[subs]
   }else{
     trapCovs <- scrFrame$trapCovs[subs]
   }
   if(!is.null(scrFrame$trapOperation)){
     trapOperation <- scrFrame$trapOperation[subs]
   }else{
     trapOperation <- scrFrame$trapOperation[subs]
   }
   
   caphist.dimensions <- sapply(caphist,dim)
   
   type <- scrFrame$type
   
   #mean maximum distance moved - mmdm
   max.dist <- NULL
   
   for (i in 1:length(caphist)) {
     for (j in 1:nrow(caphist[[i]])) {
       where <- apply(caphist[[i]][j, , ], 1, sum) > 0
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

################################################################################
