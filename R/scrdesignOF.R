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
  dmat <- e2dist(traps, statespace)
  pmat <- p0*exp(-dmat*dmat/(2*sigma*sigma))

  # crit 1: pbar
  pcaptmat <- 1-exp(-pmat)  
  pnotcaptmat <- 1-pcaptmat
  pbar <- mean(1 - exp(colSums(log(pnotcaptmat))))
  crit1 <- 1 - pbar

  
  # crit 2: p2bar (1 - pr(no capts) - pr(capt in 1 trap))
  p0times <- 1-pbar
  capattrapj <- pcaptmat
  bb <- pcaptmat/(pnotcaptmat)
  notcap <- exp(colSums(log(pnotcaptmat)))
  bling <- matrix(notcap, ncol=length(notcap), nrow=ntraps, byrow=TRUE)*bb
  p1time <- colSums(bling)
  p2bar <- mean(1 - p0times - p1time)
  crit2 <- 1-p2bar

  # crit 3: pbar + p2bar
  p12bar <- pbar + p2bar
  crit3 <- 1 - p12bar 

  # return both 1-pbar and 1-p2bar which we will minimize using SCRdesign
  valueOF <- c(crit1, crit2, crit3)[crit]
  valueOF
}

