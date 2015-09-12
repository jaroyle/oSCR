make.S <-
function(scrFrame,buffer,res,xy=NULL){
                 
 #scrFrame is an object of class 'scrFrame'
                  
 if(class(scrFrame)!='scrFrame')stop("Need an object of class scrFrame to create the state space")
                   #xy are column nummbers to locate x and y coordinates
                  
 if(is.null(xy)) xy <- c(1,2)
                   #res is the resolution for your state space
                 
  if(res==NULL) stop("You didnt provide a resolution value!")
                   #buffer is the required buffer 
                   
if(res==NULL)stop("You didnt provide a buffer value!")
                 
   S <- list()
                 
   trpls <- scrFrame$traps
                  
 for(i in 1:length(trpls)){
                     bl <- apply(trpls[[i]][,xy],2,min)
                     tr <- apply(trpls[[i]][,xy],2,max)
                     sxy <- expand.grid(seq(bl[1]-buffer,tr[1]+buffer,res),
                                        seq(bl[2]-buffer,tr[2]+buffer,res))
                     dd <- apply(e2dist(sxy,trpls[[i]][,xy]),1,min)
                     S[[i]] <- sxy[dd<=buffer,]
                   }
                    return(S)
               }
