
################################################################################
## A function that sub sets an scrFrame object

 subFrame <- function(scrFrame,subs){

    nms <- names(scrFrame)
    nms <- setdiff(nms,c("nocc","type"))
    miss.nms <- setdiff(c("caphist","traps","trapCovs","indCovs","trapOperation"),nms)
    scrFrame.sub <- list()
    for(i in nms){
      scrFrame.sub[[i]] <- list()
     if(is.null(scrFrame[[i]])){
       scrFrame.sub[[i]] <- NULL
     }else{
      for(j in 1:length(subs)){
        scrFrame.sub[[i]][[j]] <- scrFrame[[i]][[subs[j]]]
      }
     }
    }
    for(i in miss.nms){
      scrFrame.sub[[i]] <- NULL
    }
    new.nocc <- max(unlist(lapply(scrFrame.sub$caphist,function(x) dim(x)[3])))
    scrFrame.sub$nocc <- scrFrame$nocc
    scrFrame.sub$type <- scrFrame$type
    class(scrFrame.sub) <- "scrFrame"
    return(scrFrame.sub)
 }

################################################################################
