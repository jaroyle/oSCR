test.oSCR <- function(test=c(1,2,3)[1]){
  data(test)

  if(test %in% c(1,3)){
    message("Testing oSCR.fit() on the redbacked salamander data (no sex).")

    cat("\n")
  
    data(rbs)
    rbs1.sf <- subFrame(rbs.sf,2)
    rbs1.ss <- make.ssDF(rbs1.sf,res=0.5,buffer=3)
    rbs.test <- oSCR.fit(scrFrame=rbs1.sf,ssDF=rbs1.ss)
  
    cat("Output table should look like this:\n")
  
    print(rbs.known)

    cat("\n")
  
    cat("Current output table looks like this:\n")
  
    print(rbs.test)

    cat("\n")
  }
  if(test %in% c(2,3)){
    message("Testing oSCR.fit() on the ocelot data (has sex).")

    cat("\n")
  
    data(ocelot)
    ocelot.test <- oSCR.fit(scrFrame=ocelot.sf,ssDF=ocelot.ss)
  
    cat("Output table should look like this:\n")

    cat("\n")
  
    print(ocelot.known)
  
    cat("\n")
  
    cat("Current output table looks like this:\n")
  
    print(ocelot.test)
  }
}
