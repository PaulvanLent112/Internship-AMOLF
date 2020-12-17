
#Creates a distance matrix, corresponding neurons dataset and normalized count matrix for a subset of the data

#Function input: cds, neuron_dataset, num_cells_expressed
#output: a distance matrix of the subsetted 

CDS_subset<- function(cds, neuron_dataset,num_cells_expressed){
  library(monocle3)
  pData(cds)$position<- c(1:ncol(cds))
  cds_subset<- monocle3::choose_cells(cds)
  position<- pData(cds_subset)$position
  detach(package:monocle3, unload=TRUE)
  neurons2<- neuron_dataset[,position]
  #Normalization on the genes > num_cells_expressed
  library(monocle)
  expressed_genes<-which(rowSums(exprs(neurons2)>0)>num_cells_expressed)
  neurons2<- neurons2[expressed_genes,]
  
  norm<- sctransform::vst(exprs(neurons2), method="nb_fast", min_cells=num_cells_expressed,n_genes=500)$y
  dist_matrix<- bigcor(t(norm), fun="cor")
  
  dist_matrix1<- matrix(NaN, nrow=nrow(norm), ncol=nrow(norm))
  temp=1
  for (i in 1:nrow(norm)){
    temp2= temp+ nrow(norm)-1
    rowi= dist_matrix[temp:temp2]
    dist_matrix1[i,]<- rowi
    temp=temp2+1
    
  }
  dist_matrix1<- acos(dist_matrix1)/pi
  dist_matrix1[is.nan(dist_matrix1)]<-0
  dist_matrix1<-round(dist_matrix1, 6) #rond af tot 5 decimalen
  
  mylist<- list("Norm"=norm, "dist_matrix"= dist_matrix1, "neurons2"=neurons2)
  
  return(mylist)
}
