print.oSCR.fit <-
function(x, burn=NULL, ...) {
    
    
    # processing
    
    
cat(" iterations (total): ", NULL,fill=TRUE)
cat(" ",fill=TRUE)
cat(" burn-in: ", NULL, fill=TRUE)
cat(" total posterior samples: ", NULL,  fill=TRUE)
cat("-------------------------------------------------------",fill=TRUE)
cat("Posterior summaries of model parameters:","\n")

print(x$outStats,...)


}
