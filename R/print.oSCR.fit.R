print.oSCR.fit <-
function(x, burn=NULL, ...) {
mod <- x$call[["model"]]
tmpFit <- x$outStats
for(i in c(3,4,5)){
 tmpFit[,i] <- round(x$outStats[,i],3)
}
cat(" Model: ", paste(mod)[-1],fill=TRUE)
cat(" Run time: ", x$proctime," minutes",fill=TRUE)
cat(" ",fill=TRUE)
cat("-------------------------------------------------",fill=TRUE)
cat(" Summary table of model parameters:","\n")

print(tmpFit,...)
cat("-------------------------------------------------",fill=TRUE)
cat(" AIC: ", x$AIC,fill=TRUE)

}
