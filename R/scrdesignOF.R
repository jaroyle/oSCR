scrdesignOF <- function(v,          #subset/k index 
                        alltraps,   #all possible trap locations
                        statespace, #the statespace
                        N = 100,    #popsize
                        sigma,      #estimate of sigma
                        beta0,      #estimate of log(g0)
                        crit = 1    #which criterion
                        ){

  if(!"density" %in% colnames(statespace)){
    stop("The statespace must have a named 'density' column containing expected pixel densities")
  }

  statespace$pr.density <- statespace$density/sum(statespace$density)
  traps <- alltraps[v,]
  ntraps <- nrow(traps)
  
  p0 <- exp(beta0)                                #lambda0
  dmat <- e2dist(traps, statespace[,1:2])         #distance matrix (add non-euc here)
  pmat <- p0*exp(-dmat*dmat/(2*sigma*sigma))      #x-by-s detection
  
  p_iscapt_mat <- 1-exp(-pmat)                      #Pr(cap) from Poisson (hazard)
  p_nocapt_mat <- 1-p_iscapt_mat                    #Pr(cap) from Poisson (hazard)
  p_nocapt <- exp(apply(log(p_nocapt_mat),2,sum)) #not capt in any trap
  p_iscapt <- 1 - p_nocapt                          #capt in any trap

  # crit 1: pbar
  pbar <- sum(p_iscapt * statespace$pr.density)
  crit1 <- -pbar

  # crit 2: p2bar
  ##WAS THIS:
  #p0times <- 1-pbar
  #capattrapj <- pcaptmat
  #bb <- pcaptmat/(pnotcaptmat)
  #notcap <- exp(apply(log(pnotcaptmat),2,sum))
  #bling <- matrix(notcap, ncol=length(notcap), nrow=ntraps, byrow=TRUE)*bb
  #p1time <- apply(bling,2,sum) * statespace$pr.density ##this right?
  #p2bar <- mean( 1- p0times - p1time)  # pr cap 2 or more times

  # crit 2: p2bar
  ##NOW THIS:
  bb <- p_iscapt_mat/p_nocapt_mat
  bling <- matrix(p_nocapt, ncol=length(p_nocapt), nrow=ntraps, byrow=TRUE)*bb
  p1time <- apply(bling,2,sum)
  
  p2bar <- sum((1 - p_nocapt - p1time) * statespace$pr.density)
  crit2 <- -p2bar

  # crit 3: pbar + p2bar (add weighting?)
  p12bar <- pbar + p2bar
  crit3 <- -p12bar 
  
  #crit 4: MC
  firstcaps <- sum((1-exp(-apply(pmat,2,sum))) * statespace$density)
  recaps <- sum(apply(pmat,2,sum) * statespace$density) - firstcaps
  crit4 <- -recaps

  # return both 1-pbar and 1-p2bar which we will minimize using SCRdesign
  valueOF <- c(crit1, crit2, crit3, crit4)[crit]
  valueOF
}

