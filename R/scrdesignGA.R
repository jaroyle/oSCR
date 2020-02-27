# ststespace: state space points
# all.traps:  set of trap locations
# fix:        a set of points that are fixed (i.e., always in the design)
# clust.mids: centroid of a cluster (for selction purposes?)
# ntraps:     number of location we WANT
# ndesigns:   number of different 'optimal' designs to return
# nn:         number of neighbours for the swapping algorithm
# beta0:      estimated bld
# sigma:      estimated sigma
# crit: which criterion to use

scrdesignGA <- function(statespace = NULL,
                        alltraps = NULL, 
                        ntraps = 9, 
                        beta0 = -0.6, 
                        sigma = 2, 
                        crit = 1, 
                        N = 100, #why?
                        verbose=1,
                        ...){

  if(is.null(statespace)) stop("Must supply a 'statespace' object (coords of the study area)")
  if(is.null(alltraps))   stop("Must supply a 'alltraps' object (coords of all possible trap locations)")
  
  statespace <- data.frame(statespace)
  
  if(ncol(statespace)==2){
    statespace$density <- 1 / nrow(statespace)
  } 
  
  des <- kofnGA(n = nrow(alltraps), 
                k = ntraps, 
                OF = scrdesignOF,
                verbose = verbose,
                ...,
                alltraps = alltraps,
                statespace = statespace,
                beta0 = beta0,
                sigma = sigma,
                crit=crit)
  optimaltraps <- alltraps[des$bestsol,]
  scrdesign <- list(des=des, statespace=statespace, alltraps=alltraps, optimaltraps=optimaltraps, 
                    sigma=sigma, beta0=beta0)
  class(scrdesign) <- "scrdesign"
  return(scrdesign)
}
