computeclustersize<- function(lm_m){ #linkmatmerge
  add<- matrix(data=0, nrow= nrow(lm_m), ncol=3)
  lm_m<- cbind(lm_m, add)
  for (i in 1:nrow(lm_m)){
    if (lm_m[i,1]<0 & lm_m[i,2]<0){
      lm_m[i,3]<-1
      lm_m[i,4]<-1
      lm_m[i,5]<-2
    }
    else if (lm_m[i,1]<0 & lm_m[i,2]>0){
      lm_m[i,3]<-1
      lm_m[i,4]<- lm_m[lm_m[i,2],5]
      lm_m[i,5]<- 1+lm_m[lm_m[i,2],5]
    }
    else if (lm_m[i,1]>0 & lm_m[i,2]<0){
      lm_m[i,3]<- lm_m[lm_m[i,1],5]
      lm_m[i,4]<- 1
      lm_m[i,5]<- 1+ lm_m[lm_m[i,1],5]
    }
    else{
      lm_m[i,3]<-lm_m[lm_m[i,1],5]
      lm_m[i,4]<-lm_m[lm_m[i,2],5]
      lm_m[i,5]<-lm_m[lm_m[i,2],5]+lm_m[lm_m[i,1],5]
    }
  }
  return(lm_m[,3:5])
}


computeclusternumber<- function(clustersize, minsize){
  cl_num<- rep(0, nrow(clustersize))
  for (i in 1: length(cl_num)){
    if (clustersize[i,1]>=minsize & clustersize[i,2]>=minsize){
      cl_num[i+1]<- cl_num[i] -1}
    else if(clustersize[i,1]>=minsize & clustersize[i,2]<minsize){
      cl_num[i+1]<-cl_num[i]}
    else if (clustersize[i,1]<minsize & clustersize[i,2]>=minsize){
      cl_num[i+1] <- cl_num[i]
    }
    else if (clustersize[i,3]>=minsize){
      cl_num[i+1]<- cl_num[i]+1
    }
    else{cl_num[i+1]<-cl_num[i]}}
  return(cl_num[-1])
}


trackclustersizes<- function(clustersize, gn, Nmaxsizes){
  allclustersizes= matrix(0, nrow=gn, ncol=gn)
  allclustersizes[1,1]<- gn
  maxelements= matrix(0, nrow=gn, ncol=Nmaxsizes)
  maxelements[1,]<- 1
  for (i in 1:nrow(clustersize)){
    allclustersizes[i+1,]<- allclustersizes[i,]
    allclustersizes[i+1, clustersize[i,1]]<- allclustersizes[i+1, clustersize[i,1]]-1
    allclustersizes[i+1,clustersize[i,2]]<-allclustersizes[i+1,clustersize[i,2]]-1
    allclustersizes[i+1,clustersize[i,3]]<-allclustersizes[i+1,clustersize[i,3]]+1
    myallclustersizes<- allclustersizes[i+1,1:i+1]
    allpossiblesizes<-seq(1:i+1)
    allpossiblesizes<-allpossiblesizes[myallclustersizes>0]
    myallclustersizes<-myallclustersizes[myallclustersizes>0]
    temp<- c()
    macssize<-length(myallclustersizes)
    for (j in 1: min(macssize,Nmaxsizes)-1){
      temp<- c(temp, rep(allpossiblesizes[macssize-j],myallclustersizes[macssize-j]))
    }
    if (length(temp)>Nmaxsizes){
      temp<-temp[1:Nmaxsizes]}
    maxelements[i+1,1:length(temp)]<-temp
    
  }
  
  return_statement<- list(acs=allclustersizes, me=maxelements) #list of outputs
  return(return_statement)}