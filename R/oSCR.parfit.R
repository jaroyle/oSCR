oSCR.parfit <-
function(mods, data, ncores=2){
library(doParallel)

sf<- data$sf
ss<- data$ss
cs<- data$cs
 
wrapperX<- function(model.number, mods, sf, ss, cs=NULL){
 
    mod <- list(mods[[1]][[model.number]], #density
              mods[[2]][[model.number]], #detection
              mods[[3]][[model.number]], #sigma
              mods[[4]][[model.number]]) #asu
 
  fm <- oSCR.fit(scrFrame=sf, ssDF=ss, costDF=cs, model=mod, trimS=2.5,
                 se=TRUE, distmet="euc",sexmod="constant")
  save(fm, file=paste("model",model.number,".RData",sep=""))
  return(fm)
}
nmods<- nrow(mods)
cl<-makeCluster(ncores)  
registerDoParallel(cl)
out <-foreach(i=1:nmods) %dopar% {
library(oSCR)
tmp<- wrapperX(i, mods, data$sf, data$ssDF)
  return(tmp)
}
stopCluster(cl)

save("out",file="models.RData")

tmp<- list()
for(i in 1:length(out)){
    tmp[[i]]<- out[[i]]
    tmp[[i]]$call$model <- list("list",paste(mods[i,1]),paste(mods[i,2]),paste(mods[i,3]),paste(mods[i,4]))
}
out<- fitList.oSCR(tmp,rename=TRUE )
return(out)

}
