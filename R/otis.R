otis <-
function(y2d){

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
# want a table of M0 Mb Mt Mh1 Mh2
# N SE D SE AIC  (D computed if trap grid is provied)
# need to compute MMDM and so forth....


tmp<- nlm(Mhlik,c(-2,0,0 ),hessian=TRUE)
Nhat<- sum(nx) + exp(tmp$estimate[3])
SE<- sqrt(diag(solve(tmp$hessian) ))[3]
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )

AIC.Mh<- 2*tmp$minimum + 2*3

Mh.out<- c(Nhat, SE, AIC.Mh)

tmp<- nlm(Mh2lik,c(-2,0,0,0 ),hessian=TRUE)
Nhat<- sum(nx) + exp(tmp$estimate[4])
SE<- sqrt(diag(solve(tmp$hessian) ))[4]
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )

AIC.Mh2<- 2*tmp$minimum + 2*3

Mh2.out<- c(Nhat, SE, AIC.Mh2)



tmp<-(nlm(lik0,c(0,0),hessian=TRUE))
Nhat<- sum(nx) + exp( tmp$estimate[2]  )
SE<- exp(tmp$estimate[2]) * sqrt(diag(solve(tmp$hessian) ))[2]
AIC.M0<- tmp$minimum*2 + 2*2
#Dhat<- Nhat/buffers
#SE.Dhat<- sqrt( ( (1/buffers)^2 )*SE*SE  )

M0.out<- c(Nhat, SE, AIC.M0)

tmp<- rbind(M0.out,Mh.out, Mh2.out)
colnames(tmp)<- c("Nhat","SE","AIC")
tmp
}
