trackclustersizes<- function(clustersize, gn, Nmaxsizes){
  allclustersizes= matrix(0, nrow=gn, ncol=gn)
  allclustersizes[1,1]<- gn
  maxelements<- matrix(0, nrow=gn, ncol=Nmaxsizes)
  maxelements[1,]<- 1
  for (i in 1:(gn-1)){
    allclustersizes[i+1,]<- allclustersizes[i,]
    allclustersizes[i+1, clustersize[i,1]]<- allclustersizes[i+1, clustersize[i,1]]-1
    allclustersizes[i+1,clustersize[i,2]]<-allclustersizes[i+1,clustersize[i,2]]-1
    allclustersizes[i+1,clustersize[i,3]]<-allclustersizes[i+1,clustersize[i,3]]+1
    myallclustersizes<- allclustersizes[i+1,1:(i+1)]
    allpossiblesizes<-seq(1:(i+1))
    allpossiblesizes<-allpossiblesizes[myallclustersizes>0]
    myallclustersizes<-myallclustersizes[myallclustersizes>0]
    temp<- c()
    macssize<-length(myallclustersizes)
    for (j in 1: min(macssize,Nmaxsizes)){
      temp<- c(temp, rep(allpossiblesizes[macssize-j+1],myallclustersizes[macssize-j+1]))
    }
    if (length(temp)>Nmaxsizes){
      temp<-temp[1:Nmaxsizes]}
    maxelements[i+1,1:length(temp)]<-temp
    
  }
  return_statement<- list(acs=allclustersizes, me=maxelements) #list of outputs
  return(return_statement)}