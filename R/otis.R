otis <-
function(y2d,mhst=NULL){

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









Mblik<- function(parms){

}

lik0<-function(parms){
p<-  expit(parms[1])
n0<-  exp(parms[2])
N<-nind + n0
cpvec<- dbinom(0:K,K,p)
-1*(lgamma(N+1) - lgamma(n0+1) + sum(c(n0,nx)*log(cpvec) ))
}

Mhlik<-function(parms){
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


Mh2lik<-function(parms){
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



if(is.null(mhst)){
mug<- seq(-4, 1, 1)
sig<- seq(-3, 2, 0.5)
n0<- seq(0.1, 4.5, 0.6)
gr<-as.matrix(expand.grid(mug, sig, n0))

out<- matrix(NA,nrow=nrow(gr),ncol=4)
for(s in 1:nrow(gr)){
  cat("iter: ",s,fill=TRUE)
  out[s,]<- c(gr[s,],Mhlik(gr[s,]))

}

mhst<- out[out[,4]==min(out[,4]), 1:3]
}



# want a table of M0 Mb Mt Mh1 Mh2
# N SE D SE AIC  (D computed if trap grid is provied)
# need to compute MMDM and so forth....


m0parms<- c("p0","n0")
mhparms<- c("mu","sigma","n0")
mh2parms<- c("mu1","mu2","psi","n0")
mbparms<- c("mu","beta","n0")


tmp<-(nlm(lik0,c(0,0),hessian=TRUE))
Nhat<- sum(nx) + exp( tmp$estimate[2]  )
phat<- exp(tmp$estimate[1])/(1+exp(tmp$estimate[1]))
SE<- exp(tmp$estimate[2]) * sqrt(diag(solve(tmp$hessian) ))[2]
AIC.M0<- tmp$minimum*2 + 2*2
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )
M0.out<- c(Nhat, SE, AIC.M0)


tmp<- nlm(Mhlik,mhst, hessian=TRUE)
n0hat<- tmp$estimate[3]
Nhat<- sum(nx) + exp(n0hat)
muhat<- tmp$estimate[1]
sighat<- exp(tmp$estimate[2])
SE<- sqrt(diag(solve(tmp$hessian) ))[3]
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )

AIC.Mh<- 2*tmp$minimum + 2*3
Mh.out<- c(Nhat, SE, AIC.Mh)


mh2st<-c(log(phat/(1-phat)), log(phat/(1-phat))-1, 0, 1)
tmp<- nlm(Mh2lik,c(-2,0,0,0 ),hessian=TRUE)
Nhat<- sum(nx) + exp(tmp$estimate[4])
SE<- sqrt(diag(solve(tmp$hessian) ))[4]
mu1hat<- tmp$estimate[1]
mu2hat<- tmp$estimate[2]
psihat<- tmp$estimate[3]
psihat<- exp(psihat)/(1+exp(psihat))
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )

AIC.Mh2<- 2*tmp$minimum + 2*4
Mh2.out<- c(Nhat, SE, AIC.Mh2)


tmp<- rbind(M0.out,Mh.out, Mh2.out)
colnames(tmp)<- c("Nhat","SE","AIC")
tmp
}
