plot.ssDF <- function(ssDF, scrFrame=NULL, collapse=FALSE){
  
  if(!collapse){
  for(i in 1:length(ssDF)){
    ss.xy <- ssDF[[i]][,c("X","Y")]
    plot(ss.xy, pch=16, cex=0.5, col="grey", asp = TRUE)
    if(!is.null(scrFrame)){
      sf.xy <- scrFrame$traps[[i]][,c("X","Y")]  
      points(sf.xy, pch=15, col=4)
    }
  }}    
  if(collapse){
  ss.xy <- NULL  
  sf.xy <- NULL
  for(i in 1:length(ssDF)){
      ss.xy <- rbind(ss.xy, ssDF[[i]][,c("X","Y")])
      if(!is.null(scrFrame)){
        sf.xy <- rbind(sf.xy, scrFrame$traps[[i]][,c("X","Y")])  
      }
  }
  plot(ss.xy, pch=16, cex=0.5, col="grey", asp = TRUE)
  if(!is.null(scrFrame))
     points(sf.xy, pch=15, col=4)
  }   
}  
