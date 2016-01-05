print.oSCR.fit <-
function(x, burn=NULL, ...) {
mod <- x$call[["model"]]
cat(" Model: ", paste(mod)[-1],fill=TRUE)
cat(" Run time: ", x$proctime," seconds",fill=TRUE)
cat(" ",fill=TRUE)
cat("-------------------------------------------------",fill=TRUE)
cat(" Summary table of model parameters:","\n")

print(x$outStats,...)


}
