plot.scrFrame<-function(scrFrame, ask=TRUE){
  op <- par(no.readonly=TRUE)
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
      pick <- apply(tmp.ch[i,,]>0,1,sum)  
      s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),pick),c("X","Y")])
      s.ind <- c(s.ind,rep(i,sum(pick)))
      
    }
    all.ind.xy[[s]] <- data.frame(ind = s.ind, 
                                  x   = s.ind.xy[,1], 
                                  y   = s.ind.xy[,2])
    mu.x[[s]] <- tapply(all.ind.xy[[s]]$x,all.ind.xy[[s]]$ind,mean)
    mu.y[[s]] <- tapply(all.ind.xy[[s]]$y,all.ind.xy[[s]]$ind,mean)
  }
  par(oma=c(0,0,0,0))
  for(s in 1:length(scrFrame$caphist)){
    plot(scrFrame$traps[[s]][,c("X","Y")], asp=1, type="n",las=1)
    clr <- sample(colors(),nrow(tmp.ch))
    for(i in 1:nrow(tmp.ch)){
      to.x <- all.ind.xy[[s]]$x[all.ind.xy[[s]]$ind %in% i]
      to.y <- all.ind.xy[[s]]$y[all.ind.xy[[s]]$ind %in% i]
      segments(rep(mu.x[[s]][i],length(to.x)),
               rep(mu.y[[s]][i],length(to.y)), 
               to.x, to.y, col=clr[i], lwd=2)
    }
    points(scrFrame$traps[[s]][,c("X","Y")], pch="+",cex=1)
    points(mu.x[[s]],mu.y[[s]],pch=16,cex=1.5,col=clr)
  }  
  par(op)
}
