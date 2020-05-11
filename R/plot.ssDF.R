plot.ssDF <- function(ssDF, scrFrame=NULL, collapse=FALSE, spider=FALSE){
  if(!collapse){
  for(i in 1:length(ssDF)){
    ss.xy <- ssDF[[i]][,c("X","Y")]
    plot(ss.xy, pch=16, cex=0.5, col="grey", asp = TRUE,axes=FALSE,xlab="",ylab="")
    if(!is.null(scrFrame)){
      sf.xy <- scrFrame$traps[[i]][,c("X","Y")]  
      points(sf.xy, pch=15, cex=0.7, col=4)
      if(spider & !is.null(scrFrame)){
        spiderplot(scrFrame=scrFrame,session=i,add = TRUE)
      }
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
  plot(ss.xy, pch=16, cex=0.5, col="grey", asp = TRUE, axes=FALSE,xlab="",ylab="")
  if(!is.null(scrFrame))
     points(sf.xy, pch=15, cex=0.7, col=4)
  }   
}  
