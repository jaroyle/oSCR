fitList.oSCR <- function(x, rename=F, names = NULL, drop=NULL){
  
  if(rename==FALSE & is.null(names(x)))
    fl.names <- 1:length(x)
  
  if(rename==FALSE & !is.null(names(x)))
    fl.names <- names(x)
  
  if(rename==TRUE){
#    n1 <- paste(lapply(x,function(z) z$call$model[2]))
#    n1 <- gsub("\\(","",n1)
#    n1 <- gsub("\\)","",n1)
#    n1 <- gsub(" ","",n1)
    n1 <- paste(lapply(x,function(z) update.formula(formula(z$call$model[[2]]),NULL ~ .)))    
    n1 <- ifelse(n1 %in% "NULL",".",n1)
#    n2 <- paste(lapply(x,function(z) z$call$model[3]))
#    n2 <- gsub("\\(","",n2)
#    n2 <- gsub("\\)","",n2)
#    n2 <- gsub(" ","",n2)
    n2 <- paste(lapply(x,function(z) update.formula(formula(z$call$model[[3]]),NULL ~ .)))    
    n2 <- ifelse(n2 %in% "NULL",".",n2)
#    n3 <- paste(lapply(x,function(z) z$call$model[4]))
#    n3 <- gsub("\\(","",n3)
#    n3 <- gsub("\\)","",n3)
#    n3 <- gsub(" ","",n3)
    n3 <- paste(lapply(x,function(z) update.formula(formula(z$call$model[[4]]),NULL ~ .)))    
    n3 <- ifelse(n3 %in% "NULL",".",n3)
#    n4 <- paste(lapply(x,function(z) z$call$model[5]))
#    n4 <- gsub("\\(","",n4)
#    n4 <- gsub("\\)","",n4)
#    n4 <- gsub(" ","",n4)
    count.mods <- lapply(x, function(z)length(z$call$model))
    if(any(count.mods>4)){
      extract.fn <- function(z){
        if(length(z$call$model)>4)
          tmp <- unlist(update.formula(formula(z$call$model[[5]]),NULL ~ .))
        if(length(z$call$model)==4)
          tmp <- "."
        return(tmp)
      }
     n4 <- lapply(x,extract.fn)    
    }
    fl.names <- paste("D(",n1,") ","p(",n2,") ",
                      "sig(",n3,") ","asu(",n4,")", sep="")
    
  }
  names(x) <- fl.names
  class(x) <- "oSCR.fitList"
  return(x)
}

