print.oSCR.fit <-
function(x, burn=NULL, ...) {

cat(" Model: ", x$call,fill=TRUE)
cat(" Run time: ", x$proctime,fill=TRUE)
cat(" ",fill=TRUE)
cat("-------------------------------------------------",fill=TRUE)
cat("Summaries table of model parameters:","\n")

print(x$outStats,...)


}
