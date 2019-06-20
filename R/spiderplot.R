spiderplot <- function(scrFrame = NULL, session=1, y=NULL, 
                        traplocs=NULL, add=FALSE, return.stats=FALSE){
    
    if(!is.null(scrFrame)){
      traplocs <- scrFrame$trps[[session]]
      y <- scrFrame$caphist[[session]]
    }
      
    dither <- FALSE
    dx <- max(traplocs[, 1]) - min(traplocs[, 1])
    dy <- max(traplocs[, 2]) - min(traplocs[, 2])
    dx <- 0.01 * dx
    dy <- 0.01 * dy
    if(length(dim(y)) == 3) {
      if(dim(y)[2] == nrow(traplocs)) {
        nind <- dim(y)[1]
        ntraps <- dim(y)[2]
        nocc <- dim(y)[3]
        newy <- array(NA, dim = c(nind, nocc, ntraps))
        for (i in 1:nind) {
          newy[i, 1:nocc, 1:ntraps] <- t(y[i, , ])
        }
        y <- newy
      }
      y3d <- y
      J <- dim(y3d)[3]
      T <- dim(y3d)[2]
      nind <- dim(y3d)[1]
      if(add==FALSE){
        plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
      }else{
        points(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
      }
      avg.s <- matrix(NA, nrow = nind, ncol = 2)
      for(i in 1:nind){
        tmp <- NULL
        for(t in 1:T) {
          aa <- y3d[i, t, ]
          if(sum(aa) > 0){
            aa <- traplocs[aa > 0, ]
            tmp <- rbind(tmp, aa)
          }
        }
        avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
        delta <- c(runif(1, -dx, dx), runif(1, -dy, dy)) * ifelse(dither, 1, 0)
        points(avg.s[i, 1] + delta, avg.s[i, 2] + delta, pch = "S", 
               cex = 1, col = "red")
        for(m in 1:nrow(tmp)) {
          if(nrow(tmp) > 1){ 
            lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i, 2], tmp[m, 2]))
        }
      }
    }
    if(length(dim(y)) == 2) {
      y2d <- y
      J <- nrow(traplocs)
      T <- dim(y2d)[2]
      nind <- dim(y2d)[1]
      plot(traplocs, pch = 20, xlab = " ", ylab = " ", cex = 1.5)
      avg.s <- matrix(NA, nrow = nind, ncol = 2)
      for(i in 1:nind) {
        tmp <- NULL
        for (t in 1:T){
          aa <- y2d[i, t]
          if(aa <= J) {
            aa <- traplocs[aa, ]
            tmp <- rbind(tmp, aa)
          }
        }
        avg.s[i, ] <- c(mean(tmp[, 1]), mean(tmp[, 2]))
        points(avg.s[i, 1], avg.s[i, 2], pch = "S", cex = 1, col = "red")
        for(m in 1:nrow(tmp)){
          if(nrow(tmp) > 1){
            lines(c(avg.s[i, 1], tmp[m, 1]), c(avg.s[i, 2], tmp[m, 2]))
          }
        }
      }
    }
    points(traplocs, pch = 20)
    Cx <- mean(traplocs[, 1])
    Cy <- mean(traplocs[, 2])
    xcent <- sqrt((avg.s[, 1] - Cx)^2 + (avg.s[, 2] - Cy)^2)
    if(return.stats){
      return(list(xcent = xcent, avg.s = avg.s, center = c(Cx, Cy)))
    }  
  }
}
