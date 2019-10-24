plot.scrdesign <- function(scrdesign, which=4){
  nzmin <- function(x) min(x[x!=0])
  if(which==1){
    plot(scrdesign$des)
  }
  if(which==2){
    plot(scrdesign$statespace, asp=1, pch=16, cex=0.1, 
         col="grey", axes=F, xlab="", ylab="")
    box(bty="o")
    points(scrdesign$alltraps,pch=3,col=1)
    points(scrdesign$optimaltraps,pch=3,lwd=3, col=4)
  }
  if(which==3){
    mm <- max(dist(scrdesign$alltraps))+2
    hist(apply(as.matrix(dist(scrdesign$optimaltraps)),1,nzmin),
     ylim=c(0,1.5),breaks = seq(0,mm,1), xlab="Distance", main="",
     col=adjustcolor(4,0.1),freq = FALSE)
    hist(as.matrix(dist(scrdesign$optimaltraps)),
      ylim=c(0,1.5),breaks = seq(0,mm,1), xlab="Distance", main="", add=T, 
      col=adjustcolor(2,0.1), freq = FALSE)
    abline(v=c(scrdesign$sigma,2*scrdesign$sigma), lty=c(1,2), col=2)  
  }
  if(which==4){
    mm <- max(dist(scrdesign$alltraps))
    plot(scrdesign$des)
    plot(scrdesign$statespace, asp=1, pch=16, cex=0.1, 
         col="grey", axes=F, xlab="", ylab="")
    box(bty="o")
    points(scrdesign$alltraps,pch=3,col=1)
    points(scrdesign$optimaltraps,pch=3,lwd=3, col=4)
    hist(apply(as.matrix(dist(scrdesign$optimaltraps)),1,nzmin),
     ylim=c(0,1.5),breaks = seq(0,mm,1), xlab="Distance", main="",
     col=adjustcolor(4,0.1),freq = FALSE)
    hist(as.matrix(dist(scrdesign$optimaltraps)),
      ylim=c(0,1.5),breaks = seq(0,mm,1), xlab="Distance", main="", add=T, 
      col=adjustcolor(2,0.1), freq = FALSE)
    abline(v=c(scrdesign$sigma,2*scrdesign$sigma), lty=c(1,2), col=2)  
  }
}

