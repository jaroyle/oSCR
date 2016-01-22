sim.SCR.ms <- function(sessions=1, sex=F, sex.ratio=0.5, N = 100, K = 5,
                       alpha0 = -2.5, sigma = 0.5, discard0 = TRUE,
                       array3d = FALSE, ssRes = 0.5){

    caphist <- list()
    traps <- list()
    indCovs <- list()
    ssDF <- list()
    if(sessions > 1){
      if(length(sex.ratio)==1) 
        sex.ratio <- rep(sex.ratio,sessions)
      if(length(N)==1)
        N <- rep(N,sessions)
      if(length(K)==1)         
        K <- rep(K,sessions)
      if(length(alpha0)==1) 
        alpha0 <- rep(alpha0,sessions)
      if(length(sigma)==1) 
        sigma <- rep(sigma,sessions)
    }
    if(length(sex.ratio)!=sessions) stop("check length of sex ratio vector")
    if(length(N)!=sessions) stop("check length of N vector")
    if(length(K)!=sessions) stop("check length of K vector")
    if(length(alpha0)!=sessions) stop("check length of alpha0 vector")
    if(length(sigma)!=sessions) stop("check length of sigma vector")

    for(s in 1:sessions){
     dat <- sim.scr(N=N[s],K=K[s],alpha0=alpha0[s], sigma=sigma[s], discard0=TRUE, array3d=TRUE)
     Y2d <- apply(dat$Y,c(1,2),sum)
     n0<- N-nrow(Y2d)
     caphist[[s]] <- dat$Y
     traps[[s]] <- dat$traplocs; colnames(traps[[s]]) <- c("X","Y")
     indCovs[[s]] <- data.frame(sex = rbinom(nrow(dat$Y),1,sex.ratio))
     ssDF[[s]] <- dat$ss; colnames(ssDF[[s]]) <- c("X","Y")
   }
    sf <- list(caphist=caphist,
               traps=traps,
               indCovs=NULL,
               trapCovs=NULL,
               trapOperation=NULL,
               type="scr",
               nocc=K)
    if(sex==T) sf$indCovs <- indCovs
    class(sf) <- "scrFrame"
    sf <-sf
    ssDF <- ssDF
    dat <- dat
    return(list(sf=sf,ssDF=ssDF,dat=dat))
   }
