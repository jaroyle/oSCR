SCR23darray <-
function(edf, tdf){
### Returns 3-d array "ind x trap x occasion"

# Check for dups
    uq<- paste(edf[,2],edf[,3],edf[,4])
    uq<- unique(uq)
    if(any(table(uq)>1)) cat("Duplicate captures (same individual, trap, occasion) present in data set, these are not used",fill=TRUE)

nind<-max(edf[,2])
ntraps<-nrow(tdf)
nperiods<-ncol(tdf)-3
per.id<- as.numeric(dimnames(tdf)[[2]][4:ncol(tdf)])

ind.id<- edf[,2]
trap.id<- edf[,4]

if( length(per.id) != length(min(per.id):max(per.id)) ){
 x<- 1:nperiods
 names(x)<-as.character(per.id)
 per.id <- x[as.character(edf[,3])]
}
else{
per.id<-edf[,3]
}

y<-array(0,c(nind,ntraps, nperiods))

tmp<-cbind(ind.id,trap.id,per.id)
y[tmp]<-1
y
}
