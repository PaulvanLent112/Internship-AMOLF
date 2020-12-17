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

