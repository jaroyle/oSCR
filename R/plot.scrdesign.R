plot.scrdesign <- function(scrdesign, which=4){
  nzmin <- function(x) min(x[x!=0])
  if(which==1){
    plot(scrdesign$des)
  }
  if(which==2){
    plot(scrdesign$statespace[,1:2], asp=1, pch=16, cex=0.1, 
         col="grey", axes=F, xlab="", ylab="")
    box(bty="o")
    points(scrdesign$alltraps[,1:2],pch=3,col=1)
    points(scrdesign$optimaltraps[,1:2],pch=3,lwd=3, col="blue")
  }
  if(which==3){
    mm <- max(dist(scrdesign$alltraps[,1:2]))+2
    h1 <- max(hist(apply(as.matrix(dist(scrdesign$optimaltraps[,1:2])),1,nzmin),
                   breaks = seq(0,mm,mm/10),plot=FALSE)$density)
    h2 <- max(hist(as.matrix(dist(scrdesign$optimaltraps[,1:2])),
                   breaks = seq(0,mm,mm/10), plot=FALSE)$density)
    maxh <- max(h1,h2) + 0.1*max(h1,h2)
    hist(apply(as.matrix(dist(scrdesign$optimaltraps[,1:2])),1,nzmin),
         ylim=c(0,maxh),breaks = seq(0,mm,mm/10), xlab="Distance", main="",
         col=adjustcolor("navyblue",0.1),freq = FALSE)
    hist(as.matrix(dist(scrdesign$optimaltraps[,1:2])),
         ylim=c(0,1.5),breaks = seq(0,mm,mm/10), xlab="Distance", main="", 
         add=T, col=adjustcolor("maroon",0.1), freq = FALSE)
    abline(v=c(scrdesign$sigma,2*scrdesign$sigma), lty=c(1,2), col="red")
    legend("topright", 
           legend=c("Nearest neighbour","All distances","Sigma", "2 x Sigma"),
           pch = c(22,22,NA,NA), 
           col = c("black","black","red","red"),
           pt.bg = c(adjustcolor("navyblue",0.3),adjustcolor("maroon",0.3),NA,NA),
           pt.cex = c(1.75,1.75,1,1),
           lty = c(NA,NA,1,2),
           cex = 1.25)
  }
  if(which==4){
    mm <- max(dist(scrdesign$alltraps[,1:2]))+2
    plot(scrdesign$des)
    plot(scrdesign$statespace[,1:2], asp=1, pch=16, cex=0.1, 
         col="grey", axes=F, xlab="", ylab="")
    box(bty="o")
    points(scrdesign$alltraps[,1:2],pch=3,col=1)
    points(scrdesign$optimaltraps[,1:2],pch=3,lwd=3, col=4)
    h1 <- max(hist(apply(as.matrix(dist(scrdesign$optimaltraps[,1:2])),1,nzmin),
                   breaks = seq(0,mm,mm/10),plot=FALSE)$density)
    h2 <- max(hist(as.matrix(dist(scrdesign$optimaltraps[,1:2])),
                   breaks = seq(0,mm,mm/10), plot=FALSE)$density)
    maxh <- max(h1,h2) + 0.1*max(h1,h2)
    hist(apply(as.matrix(dist(scrdesign$optimaltraps[,1:2])),1,nzmin),
         ylim=c(0,maxh),breaks = seq(0,mm,mm/10), xlab="Distance", main="",
         col=adjustcolor("navyblue",0.1),freq = FALSE)
    hist(as.matrix(dist(scrdesign$optimaltraps[,1:2])),
         ylim=c(0,1.5),breaks = seq(0,mm,mm/10), xlab="Distance", main="", 
         add=T, col=adjustcolor("maroon",0.1), freq = FALSE)
    abline(v=c(scrdesign$sigma,2*scrdesign$sigma), lty=c(1,2), col="red")
    legend("topright", 
           legend=c("Nearest neighbour","All distances","Sigma", "2 x Sigma"),
           pch = c(22,22,NA,NA), 
           col = c("black","black","red","red"),
           pt.bg = c(adjustcolor("navyblue",0.3),adjustcolor("maroon",0.3),NA,NA),
           pt.cex = c(1.75,1.75,1,1),
           lty = c(NA,NA,1,2),
           cex = 1.25)
  }
}

