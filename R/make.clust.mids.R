make.clust.mids <-
function(C,S){ # could/should use cluster package for medians but doesn’t work for toy here
  mids <-
  data.frame( cl=unique(C[,3]),
              Xcoord=c(tapply(C[,1],C$cluster,mean)),
              Ycoord=c(tapply(C[,2],C$cluster,mean)))
  plot(S)
  points(C[,1:2],pch=21,bg=C[,3],cex=1.5)
  points(mids[,-1],pch="+",col=mids[,1],cex=1.5)
  return(mids)
}
