make.trapCovs <-
function(lol){
#lol = list of lists
# lol$cov1 = covariate 1 which is a list of Nsessions where lol$cov1$sess1 = ntraps x nocc etc..
# lol$cov2 = covariate 2 same format
ncovs<- length(lol)
nsessions<- length(lol[[1]])
cnames<- names(lol)
 

outlist<- list(rep(NULL, nsessions))


for(s in 1:nsessions){
ntrps<- nrow(lol[[1]][[s]])
tmp<- matrix(NA,nrow=ntrps, ncol=ncovs)
for(c in 1:ncovs){
   tmp[,c] <-   lol[[c]][[s]][,1]   # occasion 1 
}
  tmp<-  data.frame(tmp)
  colnames(tmp)<- cnames
  outlist[[s]] <- tmp
 }
 

if(ncovs>1){
for(s in 1:nsessions){
ntrps<- nrow(lol[[1]][[s]])
nocc<- ncol(lol[[1]][[s]])   # just use covariate 1 to check number of occs
for(k in 1:nocc){
tmp<- matrix(NA,nrow=ntrps, ncol=ncovs)
for(c in 1:ncovs){
   tmp[,c] <-   lol[[c]][[s]][,k]   # occasion k 
}
  tmp<-  data.frame(tmp)
  colnames(tmp)<- cnames
  outlist[[s]][[k]] <- tmp
 } # close nocc
}

} # close ncovs>1

outlist
}
