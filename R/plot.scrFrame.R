plot.scrFrame<-function(scrFrame, ask=TRUE){
  
  all.ind.xy <- list()
  mean.loc <- list()
  mu.x <- list()
  mu.y <- list()
  
  for(s in 1:length(scrFrame$caphist)){
    s.ind.xy <- NULL
    s.ind <- NULL
    tmp.ch <- scrFrame$caphist[[s]]
    tmp.tr <- scrFrame$traps[[s]]
    
    for(i in 1:nrow(tmp.ch)){
      pick <- apply(apply(tmp.ch,c(1,2),sum)>1,2,sum)  
      s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),pick),c("X","Y")])
      s.ind <- c(s.ind,rep(i,sum(pick)))
      
      }
    all.ind.xy[[s]] <- data.frame(ind = s.ind, 
                                  x   = s.ind.xy[,1], 
                                  y   = s.ind.xy[,2])
    mu.x[[s]] <- tapply(all.ind.xy[[s]]$x,all.ind.xy[[s]]$ind,mean)
    mu.y[[s]] <- tapply(all.ind.xy[[s]]$y,all.ind.xy[[s]]$ind,mean)
  }
  for(s in 1:length(scrFrame$caphist)){
    plot(scrFrame$traps[[s]][,c("X","Y")])
    for(i in 1:nrow(tmp.ch)){
      segments(mu.x[[s]][i],mu.x[[s]][i], 
               all.ind.xy[[s]]$x, all.ind.xy[[s]]$y)
    }
  }  
}
