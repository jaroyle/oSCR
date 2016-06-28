plot.ssDF <- function(ssDF, scrFrame=NULL, collapse=FALSE){
  op <- par(no.readonly=TRUE)
  par(oma=c(0,0,0,0),mar=c(0,0,0,0))
  if(!collapse){
  for(i in 1:length(ssDF)){
    ss.xy <- ssDF[[i]][,c("X","Y")]
    plot(ss.xy, pch=16, cex=0.5, col="grey", asp = TRUE,axes=FALSE)
    if(!is.null(scrFrame)){
      sf.xy <- scrFrame$traps[[i]][,c("X","Y")]  
      points(sf.xy, pch=15, cex=0.7, col=4)
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
  plot(ss.xy, pch=16, cex=0.5, col="grey", asp = TRUE, axes=FALSE)
  if(!is.null(scrFrame))
     points(sf.xy, pch=15, cex=0.7, col=4)
  }   
  #par(op)
}  
