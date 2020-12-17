require(qlcMatrix)
require(pracma)
require(zipfR)
computenoisesizes<- function(datasize1, datasize2, repeats){
  setwd("D:/Internship AMOLF/R analysis")
  sizethresh<- matrix(NaN, nrow=(datasize1-1), ncol=repeats)
  for (j in 1:repeats){
    myrand<- randn(datasize1, datasize2)
    myrandnorm<- sqrt(rowSums(myrand^2))
    myrand<- sweep(myrand, 1,STATS=myrandnorm)
    mydistmat<- acos(corSparse(t(myrand)))/pi
    mydistmat[is.nan(mydistmat)]<-0
    dist_matrix<- as.dist(mydistmat)
    
    mylinkmat=hclust(dist_matrix,"single")
    mylinkmatmerge<-mylinkmat$merge
    mylinkmat<- cbind(mylinkmatmerge, mylinkmat$height)
    
    source("functions/computeclustersize.R")
    clustersize= computeclustersize(mylinkmat)
    clustersize= as.data.frame(clustersize)
    
    
    source("functions/computeclusternumber.R")
    clustern<- computeclusternumber(clustersize, 3)  
    maxelements<-matrix(NaN, nrow=(datasize1-1),ncol=2)
    maxelement=1
    for (i in 1:(datasize1-1)){
      if (as.numeric(clustersize[i,3])>maxelement){
        maxelement= clustersize[i,3]
        maxelements[i,]<- c(clustersize[i,3],mylinkmat[i,3])
      }
    }
    maxelements<- maxelements[!is.nan(maxelements[,1]),]
    
    sizethresh[,j]<- approx(maxelements[,1], maxelements[,2],c(2:datasize1), method="linear")$y  #ask steffen

    if (j%%10==0){
      print(j)
    }

  }
  sizethreshmean= apply(sizethresh, MARGIN=1, FUN=mean)
  sizethreshstd=apply(sizethresh,MARGIN=1,FUN=sd)
  sizes<-list(sizethresh=sizethresh, sizethreshmean=sizethreshmean, sizethreshstd=sizethreshstd, maxelement=maxelements)
  return(sizes)
}



