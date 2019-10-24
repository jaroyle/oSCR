scrdesignOF <- function(v,          #subset/k index 
                        alltraps,   #all possible trap locations
                        statespace, #the statespace
                        N = 100,    #popsize
                        sigma,      #estimate of sigma
                        beta0,      #estimate of log(g0)
                        crit = 1    #which criterion
                        ){

  traps <- alltraps[v,]
  ntraps <- nrow(traps)
  p0 <- exp(beta0)
  pmat <- p0*exp(-e2dist(traps, statespace)^2/(2*sigma*sigma))

  # crit 1: pbar
  capt.pr <- 1-exp(-pmat)  
  notcapt.pr <- 1-capt.pr
  pbar <- mean(1 - exp(colSums(log(notcapt.pr))))
  crit1 <- pbar

  # crit 2: p2bar (1 - pr(no capts) - pr(capt in 1 trap))
  capt.pr <- 1-exp(-pmat)  
  notcapt.pr <- 1-capt.pr
  pbar <- mean(1 - exp(colSums(log(notcapt.pr))))
  p0times<- 1-pbar
  cap.trap.j <- 1-(notcapt.pr)
  bb <- cap.trap.j/(notcapt.pr)
  notcap <- exp(colSums(log(notcapt.pr)))
  bling <-  matrix(notcap, ncol=length(notcap), nrow=ntraps, byrow=TRUE)*bb
  p1time <- colSums(bling)
  crit2 <- mean(1 - p0times - p1time)

  # crit 3: pbar + p2bar
  crit3 <- crit1 + crit2 

  # return both 1-pbar and 1-p2bar which we will minimize using SCRdesign
  valueOF <- 1 - c(crit1, crit2, crit3)[crit]
  valueOF
}
