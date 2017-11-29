plot.scrFrame<-function(scrFrame, ax=TRUE, jit=1){
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
      if(dim(tmp.ch)[3]>1){
        pick <- apply(tmp.ch[i, , ], 1, sum) > 0
      }else{
        pick <- tmp.ch[i, , ] > 0
      }
      s.ind.xy <- rbind(s.ind.xy, tmp.tr[rep(1:nrow(tmp.tr),pick),c("X","Y")])
      s.ind <- c(s.ind,rep(i,sum(pick)))
      
    }
    all.ind.xy[[s]] <- data.frame(ind = s.ind, 
                                  x   = s.ind.xy[,1], 
                                  y   = s.ind.xy[,2])
    mu.x[[s]] <- jitter(tapply(all.ind.xy[[s]]$x,all.ind.xy[[s]]$ind,mean),factor=jit)
    mu.y[[s]] <- jitter(tapply(all.ind.xy[[s]]$y,all.ind.xy[[s]]$ind,mean),factor=jit)
  }
  for(t in 1:length(scrFrame$caphist)){
    plot(scrFrame$traps[[t]][,c("X","Y")], asp=1, type="n", las=1, 
         axes = ax, xlab = "", ylab = "")
    clr <- sample(colors(),nrow(tmp.ch))
    box(bty="o")
    for(j in 1:nrow(tmp.ch)){
      to.x <- all.ind.xy[[t]]$x[all.ind.xy[[t]]$ind %in% j]
      to.y <- all.ind.xy[[t]]$y[all.ind.xy[[t]]$ind %in% j]
      segments(rep(mu.x[[t]][j],length(to.x)),
               rep(mu.y[[t]][j],length(to.y)), 
               to.x, to.y, col=clr[j], lwd=2)
    }
    points(scrFrame$traps[[t]][,c("X","Y")], pch=3,cex=1)
    points(mu.x[[t]],mu.y[[t]],pch=21,cex=1.5,bg=clr)
  }  
}
