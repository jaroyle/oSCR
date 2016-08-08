test.oSCR <- function(){
  data(test)

  message("Testing oSCR.fit() on the redbacked salamander data (no sex).")
  data(rbs)
  rbs.test <- oSCR.fit(scrFrame=rbs.sf,ssDF=rbs.ss)
  
  cat("Output table should look like this:")
  rbs.known
  
  cat("Current output table looks like this:")
  rbs.test
  
  message("Testing oSCR.fit() on the ocelot data (has sex).")
  data(ocelot)
  ocelot.test <- oSCR.fit(scrFrame=ocelot.sf,ssDF=ocelot.ss)
  
  cat("Output table should look like this:")
  ocelot.known
  
  cat("Current output table looks like this:")
  ocelot.test
}
