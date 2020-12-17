#input: cluster nummer en dataset
complement<- function(cluster, CDS){
  sam<-pData(CDS)
  r<-sam%>%
    mutate(pos=c(1:nrow(sam)))%>%
    filter(Cluster==cluster)%>%
    select(pos) #This splits the data set into two groups, one of the cluster and the remaining clusters
  r<- as.vector(r$pos)
  setx<- CDS[,r]
  setx<- as.data.frame(exprs(setx))
  setx<- rowMeans(setx)
  setx<-setx[setx>1]
  setx<- names(setx)# these are the names of expressed modules in cluster x
  #The rest of the clusters
  sety<-CDS[,-r]
  sety<- as.data.frame(exprs(sety))
  sety<- rowMeans(sety)
  sety<- sety[sety>1]
  sety<- names(sety)
  #determine the complement
  complement<- setdiff(setx, sety)
  complement<- setdiff(complement, sety)
  return(complement)
}
