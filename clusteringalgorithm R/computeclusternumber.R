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
