
extractallclusters<- function(gn, clustersize, mylinkmat, mydim, Nref, falsecluster, falsemembers, output){
  require(pracma)
  require(zipfR)
  
  #Make supportfolder
  supportfolder<- paste0(getwd(),"/PercolationReference")
  if (dir.exists(supportfolder)){
  }
  else{dir.create("~/PercolationReference")}
  
  #load reference for percolation transition
  if (file.exists(paste0(supportfolder,'/PercolationDimensionReference.mat'))){
    percol<- R.matlab::readMat("PercolationDimensionReference.mat")
    dimpercolsmooth<- percol$dimpercolsmooth
    mydimlist<- percol$mydimlist
  }
  else{
    print("Warning: Reference file is missing")
    print("Assuming percolation at degree 1, not valid for D<20.")
    dimpercolsmooth<- c(1,1)
    mydimlist<- c(4,100)
  }
  inflectionpoint<- which.min(-diff(dimpercolsmooth) %/% diff(mydimlist))
  dimpercolsmooth<- approx(c(mydimlist[1:inflectionpoint],10^12), c(dimpercolsmooth[1:inflectionpoint],1),seq(4,mydim))
   
  #load reference for largest sizes
  
  stl_file<- paste0(getwd(),"/sizethreshlog",as.character(gn),"-",as.character(Nref),".csv")
  if (file.exists(stl_file)){
    sizethreshlog<- read.csv(stl_folder)
    if (nrow(stlsize)<mydim-3){
      sizethreshlog[(stlsize+1):(mydim-3),]<-matrix(list(),nrow(mydim-3-stlsize),ncol=2)}
  sizethreshlist= matrix(0, nrow=1, ncol=nrow(sizethreshlog))
  for (i in 1:nrow(sizethreshlog)){
    sizethreshlist[i]<- !isempty(sizethreshlog[[i,1]]) #????
  }
  }
  else{
    sizethreshlog<- matrix(list(),nrow=mydim-3,ncol=2)
    sizethreshlist<- rep(0, mydim-3)
  }

  
  
  #trackclusterorigin
  
  #REMOVED THE MYLINKMAT-GN and changed the if-statements to as.numeric(x) to avoid errors 
  #Now this part works
  trackclustersize= seq(NaN, gn-1)
  icount=0
  for (i in 1:(gn-1)){
    if (as.numeric(clustersize[i,3])==2){
      icount=icount+1
      trackclustersize[i]=icount
    }
    else if (as.numeric(clustersize[i,1])<as.numeric(clustersize[i,2])){
      trackclustersize[i]=trackclustersize[mylinkmat[i,2]]
    }
    else if (as.numeric(clustersize[i,1])>as.numeric(clustersize[i,2])){
      trackclustersize[i]= trackclustersize[mylinkmat[i,1]]
    }
    else if (as.numeric(trackclustersize[mylinkmat[i,1]])>as.numeric(trackclustersize[mylinkmat[i,2]])){
      trackclustersize[i]=trackclustersize[mylinkmat[i,2]]
    }
    else{
      trackclustersize[i]=trackclustersize[mylinkmat[i,2]]
    }
  }
  

  
  
  #actualclustering
  maxpercolthreshold=fzero(function(x) gn* Rbeta(sin(x*pi)^2, (mydim-2)/2, 0.5)/2- dimpercolsmooth[mydim-3],
                           mylinkmat[nrow(mylinkmat),3])  #Is this correct?
  
  maxpercolind<- which.min(abs(maxpercolthreshold-mylinkmat[,3]))
  clusterIDreal<- seq(NaN, icount)
  
  #scan from last merging to first merging
  
  for (i in maxpercolind:3){
    mergepoints<- mylinkmat[i, 1:2]
    if (sum(mergepoints>0)==2){ #Only if clusters merge this can be a percolation point for the two branches, since then they are real clusters 
      myclusterIDreal<- clusterIDreal[trackclustersize[mergepoints-gn]]
      #checks whether ID of both branches is associated with real cluster
     # clusterflag<- computingclusterflag(i, myclusterIDreal)
      computingclusterflag<-source("functions/computingclusterflag")
      clusterflag<- computingclusterflag(i, myclusterIDreal)
      if (clusterflag<2){
        if (trackclustersize[mergepoints[1]-gn]==trackclustersize[i]){
          smind=2
          largind=1
          }
        else{
          smind=1
          largind=2
        }
        if (clusterflag==1){
          clusterIDreal<- clusterIDreal[trackclustersize[mergepoints[largind]-gn]]
        }
        mydim=which.min(abs(gn*Rbeta(sin(mylinkmat[i,3]*pi)^2,(c(4:mydim)-2)/2,1/2)/2-dimpercolsmooth[1:mydim-3]))
          #line 118
        mydim=mydim+3
        
        if (sizethreshlist[mydimlist-3]==0){
          size_compute<-computenoisesizes5(gn,mydim,Nref)
          sizethreshstd<- size_compute$sizethreshstd
          sizethreshlog[[(mydim-3),1]]<-size_compute$sizethreshmean
          sizethreshlog[[(mydim-3),2]]<-size_compute$sizethreshstd
          sizethreshlist[mydim-3]<- 1
        }
        else{
          sizethreshmean<- sizethreshlog[[(mydim-3),1]]
          sizethreshstd<- sizethreshlog[[(mydim-3),1]]
        }
        sizethresh<- sizethreshmean -(falsecluster *sizethreshstd)
        sizethresh[1]<-0 #Dont use clusters of size 2
        
        #check whether smaller cluster includes one or more real clusters
        checkIDlist<- mergepoints[smind]-gn
        checkflag=0
        
        while (!isempty(checkIDlist)){
          #Check whether smaller cluster includes on or more real clusters
          abovethresh1<- as.matrix(0, nrow=0, ncol=checkIDlist[1])
          #find all instances with the same ID
          temp<- trackclustersize[1:checkIDlist[1]==trackclustersize[checkIDlist[1]]] #this probably doesnt work
          mymergings<- mylinkmat[1:checkIDlist[1],3]
          myclustersize<- clustersize[1:checkIDlist[1],3]
          abovethresh1[temp]<- cumsum(mymergings[temp]<sizethresh[myclustersize[temp]-1])>0
          #line 157
          if (sum(abovethresh1)>0){
            sizethreshpos<- round(falsemembers*myclustersize[temp])-1
            sizethreshpos[sizethreshpos<1]=1
            sizethreshpos[sizethreshpos>=gn]= gn-1
            abovethresh2= as.matrix(0, nrow=1, ncol=checkIDlist[1])
            abovethresh2[temp]= mymergings[temp]<sizethreshmean[sizethreshpos]
            realclustpos= finds((abovethresh1*abovethresh2)==1,1) 
            realclustpos= realclustpos[length(realclustpos)] #R
            if (!isempty(realclustpos)){
              clusterIDreal[trackclustersize[checkIDlist[1]]]<-realclustpos
              checkflag=1
            }
            #if criterion is not met, set realclusterpos to beginning
            #(this allows to identify small real clusters on long
            #branches (very tight) that however merge into a larger
            #a non-real cluster before the percolation threshold
            else{
              realclustpos<- finds(trackclustersize==trackclustersize[checkIDlist[1]])
              realclustpos<- realclustpos[1]
              
            }}
          #If criterion is not met, set realclusterpos to beginning of cluster
          else{
            realclustpos<- finds(trackclustersize==trackclustersize[checkIDlist[1]])
            realclustpos<- realclustpos[1]
          }
          #line 181
          if (realclustpos<checkIDlist[1]){ #there are branchings above
            #find all branchings above
            tempnew=temp * matrix(0, nrow=realclustpos, ncol=1)# not sure what is going on here
            if (sum(tempnew)==0){
              print('error tempnew')
            }
            mymergingspos<-mylinkmat[1:checkIDlist[1],1:2]
            mymergingspos<- mymergingspos[tempnew>0,]
            mymerginspos<- mymergingspos[colSums(mymergingspos>gn)==2,]
            if (!isempty(mymergingspos)){ #cluster branching in
              checkIDlisttemp<- Reshape(mymergingspos-gn,c(),1) #???
              #Remove branches that have been just considered
              checkIDlisttemp[trackclustersize[checkIDlisttemp]==trackclustersize[checkIDlist[1]]]=c() #???
              checkIDlist=c(checkIDlist, checkIDlisttemp) #?????
              
            }
          }
            checkIDlist[1]=c()
          } #line 204
          if (checkflag==1){ #small cluster was real
            checkIDlist= mergepoints[largind]-gn
            while (!isempty(checkIDlist)){
              #Check whether larger cluster was above size threshold before merging
              abovethresh1<- matrix(0, nrow=1, ncol=checkIDlist)
              temp<- trackclustersize[1:checkIDlist[1]==trackclustersize[checkIDlist[1]]] #????
              mymergings<- mylinkmat[1:checkIDlist[1],3]
              myclustersize<- clustersize[1:checkIDlist[1],3]
              abovethresh1[temp]<-cumsum(mymergings[temp]<sizethresh[myclustersize[temp]-1])>0
              if (sum(abovethresh)>0){#only if the first criterion is met
                sizethreshpos<- round(falsemembers*myclustersize[temp])-1
                sizethreshpos[sizethreshpos<1]=1
                sizethreshpos[sizethreshpos>=gn]= gn-1
                abovethresh2= as.matrix(0, nrow=1, ncol=checkIDlist[1])
                abovethresh2[temp]= mymergings[temp]<sizethreshmean[sizethreshpos]
                realclustpos= finds((abovethresh1*abovethresh2)==1,1) 
                realclustpos= realclustpos[length(realclustpos)] #R
                if (!isempty(realclustpos)){
                  clusterIDreal[trackclustersize[checkIDlist[1]]]<-realclustpos
                  checkflag=1
                }
                else{
                  realclustpos<- Finds(trackclustersize==trackclustersize[checkIDlist[1]])
                  realclustpos<- realclustpos[1]
                }
              }
              else{
                realclustpos<- Finds(trackclustersize==trackclustersize[checkIDlist[1]])
                realclustpos<- realclustpos[1] #line 233
              }
              
              
              if (realclustpos<checkIDlist[1]){#there are branchings above
                #find all branchings above
                tempnew<- temp * matrix(0, nrow=realclustpos, ncol=1)#ask
                if (sum(tempnew)==0){
                  print("Error. Tempnew")
                }
                mymergingspos<- mylinkmat[1:checkIDlist[1],1:2]
                mymergingspos<- mymergingspos[tempnew>0,]
                mymergingspos<- mymergingspos[colSums(mymergingspos>gn)==2,]
                if (!isempty(mymergingspos)){#clusters branching in
                  checkIDlisttemp<- Reshape() #Ask?
                  #remove branches that have been just considered
                  checkIDlisttemp[trackclustersize[checkIDlisttemp]==trackclustersize[checkIDlist[1]]]=c() #???
                  checkIDlist<- c(checkIDlist,checkIDlisttemp) #???
                }}
              checkIDlist[1]<- c()
            }}
          #if there was a subclustering but the large cluster did not
          #change, remove the large cluster event above
          if (checkflag>0 & clusterflag==1){
            if (clusterIDrealref==clusterIDreal[trackclustersize[mergepoints[largind]-gn]]){
              clusterIDreal[trackclustersize[mergepoints[largind]-gn]]<-NaN  } 
          }
        
      }}}
  
    allclusters<- clusterIDreal[!is.nan(clusterIDreal)]
    
    #Remove later
    myclusterIDs<- trackclustersize[allclusters]
    myclusterIDsUnique<- unique(myclusterIDs)
    myclusterIDsUniquesize<- length(myclusterIDsUniqie)
    if (myclusterIDsUniquesize<sum(!is.nan(clusterIDreal))){
      print("Error. Clusters not distinct!")
    }
    
    #Save reference for largest cluster sizes (287)
    
    if (output>0){
      print(paste0("We have found ",sum(!is.nan(clusterIDreal)), " clusters."))
    }
    #We have to return something
}


  
    

    
        
              
                
              
              
            
            
          
          
  

