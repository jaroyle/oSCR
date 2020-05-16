 subFrame <- function(scrFrame,subs){
   
   caphist <- scrFrame$caphist[subs]
   traps <- scrFrame$traps[subs]

   if(!is.null(scrFrame$indCovs)){
     indCovs <- scrFrame$indCovs[subs]
   }else{
     indCovs <- NULL
   }
   
   if(!is.null(scrFrame$trapCovs)){
     trapCovs <- scrFrame$trapCovs[subs]
   }else{
     trapCovs <- NULL
   }
   
   if(!is.null(scrFrame$sigCovs)){
     pick.rows <- scrFrame$sigCovs$session %in% subs
     sigCovs <- scrFrame$sigCovs[pick.rows,,drop=FALSE]
   }else{
     sigCovs <- data.frame(session = factor(1:length(subs)))
   }
   
   if(!is.null(scrFrame$trapOperation)){
     trapOperation <- scrFrame$trapOperation[subs]
   }else{
     trapOperation <- scrFrame$trapOperation[subs]
   }

   if(!is.null(scrFrame$telemetry)){
     telemetry <- scrFrame$telemetry
   }else{
     telemetry <- NULL
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
   mdm <- max(max.dist,na.rm=T)
   
   scrFrame <- list("caphist" = caphist,
                    "traps" = traps,
                    "indCovs" = indCovs,
                    "trapCovs" = trapCovs,
                    "sigCovs" = sigCovs,
                    "trapOperation" = trapOperation,
                    "occasions" = caphist.dimensions[3,],
                    "type" = type,
                    "mmdm" = mmdm,
                    "mdm" = mdm,
                    "telemetry" = telemetry)
   
   class(scrFrame) <- "scrFrame"  
   return(scrFrame)
 }

################################################################################
