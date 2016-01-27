

#modSel - a function for generating an ordered (by deltaAIC) model selection
#         table from an 'oSCR.fitList'

# x: an 'oSCR.fitList object

modSel.oSCR <- function(x){

  if(class(x)=="oSCR.fitList"){
    ms <- list()
    #AIC table
    df.out <- data.frame(model = names(x),
                         logL = unlist(lapply(x, function(y) y$rawOutput$minimum)),
                         K = unlist(lapply(x, function(y) length(y$rawOutput$estimate))),
                         AIC = unlist(lapply(x, function(y) y$AIC)))

    df.out$dAIC <- df.out$AIC - min(df.out$AIC)
    df.out2 <- df.out[order(df.out$dAIC),]
    rownames(df.out2) <- NULL
    ms[["aic.tab"]] <- df.out2

    #coefficient table
    coef.df <- NULL
    aic.df <- NULL
    for(i in 1:length(x)){
      tmp.df <- x[[i]]$coef.mle
      tmp.df$model <- names(x)[i]
      coef.df <- rbind(coef.df,tmp.df)
    }
    coef.tab <- data.frame(tapply(coef.df$mle,list(coef.df$model,coef.df$param),unique))
    coef.tab$model <- rownames(coef.tab)
    coef.tab <- coef.tab[,c("model",setdiff(colnames(coef.tab),"model"))]
    rownames(coef.tab) <- NULL
    coef.out <- merge(coef.tab,df.out[,c("model","AIC")],by="model")

    ms[["coef.tab"]] <- coef.out
    class(ms) <- "oSCR.modSel"
    return(ms)
  }else{
    print("Object is not of class oSCR.fit or oSCR.fitList")
  }
}


