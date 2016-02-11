otis <-
function(y2d,mhst=NULL, cov=NULL){

K<- ncol(y2d)
y<- apply(y2d,1,sum)
nx<- rep(0, K)
names(nx)<- 1:K
yt<- table(y)
nx[names(yt)]<- yt
nind<- sum(nx)


expit<- function(x){
 1/(1+exp(-x))
}


prevcap<- matrix(0,nrow=nrow(y2d),ncol=K)
first<-rep(NA,nind)
for(i in 1:nind){
first[i]<- min(  (1:K)[y2d[i,]>0]   )
if(first[i]<K)
    prevcap[(first[i]+1):K,]<- 1
}






Mblik<- function(parms, behave=FALSE, cov=NULL){

if(!is.null(cov))
  covmat<- matrix(cov, ncol=K,nrow=nrow(prevcap) , byrow=TRUE)

if(behave == TRUE & is.null(cov)){
p<-  expit(parms[1]+ parms[2]*prevcap)
n0<-  exp(parms[3])
p<- rbind(p, rep(expit(parms[1]), K) )
}
if(behave==TRUE & !is.null(cov)){
  p<-  expit(parms[1]+ parms[2]*prevcap + parms[3]*covmat)
  n0<- exp(parms[4])
  p<- rbind(p,  expit(parms[1] + parms[3]*cov)  )
}
if(behave==FALSE & !is.null(cov)){
  p<-  expit(parms[1]+ parms[2]*covmat)
  n0<- exp(parms[3])
  p<- rbind(p,  expit(parms[1] + parms[2]*cov)  )
}


N<-nind + n0
y0<- rbind(y2d, rep(0,K))
llik<- rep(NA,nrow(y2d)+1    ) 
for(i in 1:nrow(y0)){
   llik[i]<- sum(dbinom(y0[i,],1,p[i,],log=TRUE))
}
-1*(lgamma(N+1) - lgamma(n0+1) + sum( c(rep(1,nind),n0)*llik )  )
}

M0lik<- function(parms){
 p<-  matrix(expit( parms[1] ), nrow=nind,ncol=K)
 n0<-  exp(parms[2])
 N<-nind + n0
 y0<- rbind(y2d, rep(0,K))
 p<- rbind(p, rep(expit(parms[1]), K) )
 llik<- rep(NA,nrow(y2d)+1    ) 
 for(i in 1:nrow(y0)){
     llik[i]<- sum(dbinom(y0[i,],1,p[i,],log=TRUE))
  }
  -1*(lgamma(N+1) - lgamma(n0+1) + sum( c(rep(1,nind),n0)*llik )  )
}
# old
lik0<-function(parms){
p<-  expit(parms[1])
n0<-  exp(parms[2])
N<-nind + n0
cpvec<- dbinom(0:K,K,p)
-1*(lgamma(N+1) - lgamma(n0+1) + sum(c(n0,nx)*log(cpvec) ))
}

Mhlikxx<-function(parms){
mu<- parms[1]
p0<- exp(mu)/(1+exp(mu))
sigma<-exp(parms[2])
n0<-exp(parms[3])

il<-rep(NA,K+1)
for(k in 0:K){
il[k+1]<-integrate(
function(x){
 dbinom(k,K,expit(mu+x) )*dnorm(x,0,sigma)
},
lower=-Inf,upper=Inf)$value
}

-1*(    lgamma(n0+nind+1) - lgamma(n0+1) + sum(c(n0,nx)*log(il)))
}



# Likelihood function with sigma fixed (near 0) result should be close to Model M0 result
Mhlik <-function(parms, y2d ){

# Just two parameters
mu<- parms[1]
sigma<-exp(parms[2])
n0<-exp(parms[3])
K<- ncol(y2d)
 
nind<- nrow(y2d)

# add a row of zeros onto the matrix
 y2d<- rbind(y2d, rep(0,K))
# a vector to hold the result
 il<-rep(NA,nrow(y2d))

# loop over rows of y2d and compute the marginal likelihood by evaluating the integral on a fine grid
# (note: the "integrate" function doesn't work either..... )
for(i in 1:nrow(y2d) ){
 ym <- as.vector(y2d[i,])
 eta.gr<- seq(-8,8,.1) 
 f <-  function(eta, sig = sigma){
           #dbinom( sum(ym),K, expit(mu+eta), log=FALSE)*dnorm(eta,0,sig)
          apply(sapply(expit(mu+eta), dbinom, x = ym, size = 1),2,prod) * dnorm(eta,0,sig)        
     }
  f.gr<- f(eta.gr)
  # Marginal likelihood right here:
  il[i]<- sum(f.gr*0.1)
}
# Sum this stuff up
-1*(    lgamma(n0+nind+1) - lgamma(n0+1) + sum(c(rep(1,nind),n0)*log(il)))
}




Mh2likxxx<-function(parms){
mu1<- parms[1]
mu2<- parms[2]
p1<- exp(mu1)/(1+exp(mu1))
p2<- exp(mu2)/(1+exp(mu2))

psi<-exp(parms[3])/(1+exp(parms[3]))
n0<-exp(parms[4])

il<-rep(NA,K+1)
for(k in 0:K){
   il[k+1]<-  dbinom(k,K,p1)*psi + dbinom(k,K,p2)*(1-psi)
}

-1*(    lgamma(n0+nind+1) - lgamma(n0+1) + sum(c(n0,nx)*log(il)))
}

Mh2lik<-function(parms, y2d, cov=NULL){
if(is.null(cov)){
mu1<- parms[1]
mu2<- parms[2]
p1<- exp(mu1)/(1+exp(mu1))
p2<- exp(mu2)/(1+exp(mu2))
psi<-exp(parms[3])/(1+exp(parms[3]))
n0<-exp(parms[4])
}


nind<- nrow(y2d); K<- ncol(y2d)
y2d<- rbind(y2d, rep(0, K))

il<-rep(NA,nrow(y2d) )
for(i in 1:nrow(y2d)){
   il[i]<-  exp(sum(dbinom(y2d[i,],1,p1,log=TRUE)))*psi +   exp(sum(dbinom(y2d[i,],1,p2,log=TRUE)))*(1-psi)
}

-1*(    lgamma(n0+nind+1) - lgamma(n0+1) + sum(c(rep(1,nind), n0)*log(il)))
}
 

if(is.null(mhst)){
mug<- seq(-4, 1, 1)
sig<- seq(-3, 2, 0.5)
n0<- seq(0.1, 4.5, 0.6)
gr<-as.matrix(expand.grid(mug, sig, n0))

out<- matrix(NA,nrow=nrow(gr),ncol=4)
for(s in 1:nrow(gr)){
  cat("iter: ",s,fill=TRUE)
  out[s,]<- c(gr[s,],Mhlik(gr[s,], y2d=y2d))

}

mhst<- out[out[,4]==min(out[,4]), 1:3]
}

 

# want a table of M0 Mb Mt Mh1 Mh2
# N SE D SE AIC  (D computed if trap grid is provied)
# need to compute MMDM and so forth....





tmp<-(nlm(M0lik,c(0,0),hessian=TRUE))
Nhat<- sum(nx) + exp( tmp$estimate[2]  )
phat<- exp(tmp$estimate[1])/(1+exp(tmp$estimate[1]))
SE<- exp(tmp$estimate[2]) * sqrt(diag(solve(tmp$hessian) ))[2]
AIC.M0<- tmp$minimum*2 + 2*2
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )
M0.out<- c(tmp$minimum, length(tmp$estimate), Nhat, SE, AIC.M0)
M0.parms<- tmp$estimate


tmp<- nlm(Mhlik,mhst, y2d=y2d,hessian=TRUE)
n0hat<- tmp$estimate[3]
Nhat<- sum(nx) + exp(n0hat)
muhat<- tmp$estimate[1]
sighat<- exp(tmp$estimate[2])
SE<- sqrt(diag(solve(tmp$hessian) ))[3]
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )

AIC.Mh<- 2*tmp$minimum + 2*3
Mh.out<- c(tmp$minimum, length(tmp$estimate), Nhat, SE, AIC.Mh)
Mh.parms<- tmp$estimate

 

mh2st<-c(log(phat/(1-phat)), log(phat/(1-phat))-1, 0, log(n0hat))
tmp<- nlm(Mh2lik,mh2st,y2d=y2d, hessian=TRUE)
Nhat<- sum(nx) + exp(tmp$estimate[4])
SE<- sqrt(diag(solve(tmp$hessian) ))[4]
mu1hat<- tmp$estimate[1]
mu2hat<- tmp$estimate[2]
psihat<- tmp$estimate[3]
psihat<- exp(psihat)/(1+exp(psihat))
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )
AIC.Mh2<- 2*tmp$minimum + 2*4
Mh2.out<- c(tmp$minimum, length(tmp$estimate), Nhat, SE, AIC.Mh2)
Mh2.parms<- tmp$estimate

## model Mcov only
if(!is.null(cov)){
st<-c(log(phat/(1-phat)), 0, log(n0hat)-1 )

tmp<- nlm(Mblik,st,behave=FALSE, cov=cov, hessian=TRUE)
Nhat<- sum(nx) + exp(tmp$estimate[3])
SE<- sqrt(diag(solve(tmp$hessian) ))[3]
beta0hat<- tmp$estimate[1]
beta1hat<- tmp$estimate[2]
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )
AIC.Mcov<- 2*tmp$minimum + 2*3
Mcov.out<- c(tmp$minimum, length(tmp$estimate), Nhat, SE, AIC.Mcov)
Mcov.parms<- tmp$estimate
}else{
Mbcov.out<- c(NULL, NULL, NULL, NULL, NULL)
Mbcov.parms<- rep(NULL, 4)
}

## model Mb
mbst<-c(log(phat/(1-phat)), 0, log(n0hat)-1 )
tmp<- nlm(Mblik,mbst,behave=TRUE, hessian=TRUE)
Nhat<- sum(nx) + exp(tmp$estimate[3])
SE<- sqrt(diag(solve(tmp$hessian) ))[3]
beta0hat<- tmp$estimate[1]
beta1hat<- tmp$estimate[2]
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )
AIC.Mb<- 2*tmp$minimum + 2*3
Mb.out<- c(tmp$minimum, length(tmp$estimate), Nhat, SE, AIC.Mb)
Mb.parms<- tmp$estimate

if(!is.null(cov)){
## model Mb+cov
mbst<-c(log(phat/(1-phat)), 0, 0, log(n0hat)-1 )
tmp<- nlm(Mblik,mbst,hessian=TRUE, behave=TRUE, cov=cov)
Nhat<- sum(nx) + exp(tmp$estimate[4])
SE<- sqrt(diag(solve(tmp$hessian) ))[4]
beta0hat<- tmp$estimate[1]
beta1hat<- tmp$estimate[2]
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )

AIC.Mbcov<- 2*tmp$minimum + 2*4
Mbcov.out<- c(tmp$minimum, length(tmp$estimate), Nhat, SE, AIC.Mbcov)
Mbcov.parms<- tmp$estimate
}else{
Mbcov.out<- c(NULL, NULL, NULL, NULL, NULL)
Mbcov.parms<- rep(NULL, 4)
}

M0.parms
Mh.parms
Mh2.parms
Mcov.parms
Mb.parms
Mbcov.parms
names(M0.parms) <-   c("mu","n0")
names(Mh.parms) <-   c("mu","logsigma","n0")
names(Mh2.parms)<-   c("mu1","mu2","logit(psi)","n0")
names(Mcov.parms)<-  c("mu","beta","n0")
names(Mb.parms)<-    c("mu","behave","n0")
names(Mbcov.parms)<- c("mu","behave","beta","n0")
unam<- unique(c(names(M0.parms),names(Mh.parms),names(Mh2.parms),names(Mcov.parms),names(Mb.parms),names(Mbcov.parms)))
coef<- matrix(NA,nrow=6,ncol=length(unam))
colnames(coef)<- unam
mods<- c("p(.)","p(h; ln)", "p(h; fm)","p(b)","p(b+cov)","p(cov)")
coef[1,names(M0.parms)]<- M0.parms
coef[2,names(Mh.parms)]<- Mh.parms
coef[3,names(Mh2.parms)]<- Mh2.parms
coef[4, names(Mb.parms)]<- Mb.parms
coef[5,names(Mbcov.parms)]<- Mbcov.parms
coef[6,names(Mcov.parms)]<- Mcov.parms

tmp<- rbind(M0.out,Mh.out, Mh2.out,Mb.out, Mbcov.out,Mcov.out)
colnames(tmp)<- c("logL","K", "Nhat","SE","AIC")
tmp<- data.frame(tmp)

aic.tab<- cbind("model" = mods, tmp)
coef<- data.frame(coef)
coef<- cbind("model" = mods, coef)
coef<- coef[order(aic.tab[,"AIC"]),]


aic.tab<- aic.tab[order(aic.tab[,"AIC"]),]

list(aic.tab=aic.tab, coef.tab=coef )

}

















