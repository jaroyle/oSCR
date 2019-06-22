p2bar.fn <-
function(X, S, N = 100, sigma, beta0 = log( K*p0 ) ){
    # compute pbar

    p0<- exp(beta0)

    gr<- S
 traps<- X
ntraps<- nrow(X)
 pp<- e2dist(traps, gr)
 pp<- p0*exp(- pp*pp/(2*sigma*sigma))  # ntraps x nrow(statespace)
# hazard
 pp<- 1-exp(-pp)
    notpp<- 1-pp  # same :  ntraps x nrow(statespace)

# pbar<-   mean(  1- exp(colSums(log(notpp^K)))    )

# pbar =  1- Pr(not captured in any trap) = 1  - ( prod_{j} 1-p[j] )
 notcap <-          exp( colSums(log(notpp)))  # not captured in any trap
 pbar<-   mean(  1- exp(colSums(log(notpp)))    )  # scalar
 p0times<- 1-pbar

 #notcap <-  exp( colSums(log(notpp^K)) )
 #cap.trap.j<- 1-(notpp)^K
 #bb<- cap.trap.j/((1-pp)^K)

 # I forget what all of this is, I think the probability of being captured only in trap j
 cap.trap.j<- 1-(notpp)
 bb<- cap.trap.j/((1-pp))
 bling<-  matrix(notcap, ncol=length(notcap), nrow=ntraps, byrow=TRUE)*bb

 # Probability of being captured in exactly 1 trap, any trap
 p1time<- colSums(bling)

 # p2bar = 1 minuts probability of being captured no times minus probability of being caught in exactly 1 trap
    p2bar <- mean( 1- p0times - p1time)  # pr cap 2 or more times

    # return both 1-pbar and 1-p2bar which we will minimize using SCRdesign
c(1-pbar,1-p2bar)
}
