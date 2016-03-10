predict.oSCR <-
function(scrFrame, scr.fit, ssDF, costDF=NULL){
 
  mles <- scr.fit$rawOutput$estimate
  call <- scr.fit$call
  
  oSCR.fit2 <- oSCR.fit
  # modify the arguments to match the original call (skip the data)
  call.fix <- names(call)[!names(call) %in% c("","scrFrame","ssDF","costDF")]
  for(i in 1:length(call.fix)){
    formals(oSCR.fit2)[[call.fix[i]]]<- call[[call.fix[i]]]
  }
  prediction <- oSCR.fit2(scrFrame=scrFrame, ssDF=ssDF, costDF=costDF,
                          #critical arguments that differ for prediction
                          start.vals=mles, predict=TRUE)
  return(prediction)
 
}
