computingclusterflag<- function(mergepoint, myclusterIDreal){
  clusterflag<- rep(0, 2)
  if (!is.nan(myclusterIDreal[1])){
    if (myclusterIDreal[1]<mergepoint){
      clusterflag[1]=2
    }
    else {clusterflag[1]=1}
  }
  if (!is.nan(myclusterIDreal[2])){
    if (myclusterIDreal[2]<mergepoint) {clusterflag[2]=2}
    else {clusterflag[2]=1}
  }
  clusterflag<- sum(clusterflag)
  return(clusterflag)
  }