Cluster_expression<- function(cds, clusterindex){
  #The expression of gene modules for a certain cluster
  clusters<- as.numeric(clusters(cds))
  ind<- which(clusters==clusterindex)
  sub<- cds[,ind]
  expr_sub<- rowMeans(exprs(sub))
  group<-data.frame(sheets=rownames(sub), expr=expr_sub)
  group<- group[order(group$expr, decreasing=TRUE),]
  group<- group[1:10,]
  
  plot<- ggplot(group)+
    geom_col(mapping=aes(x=sheets, y=expr), fill="gold1")+
    xlab("Modules")+
    ylab("Normalized expression")+
    ggtitle(paste0("Cluster ",clusterindex))+
    theme_classic()
  
  return(plot)
}

Neuron_expression<- function(cds,Neurontype){
  neurontypes<- cds$Neuron.type
  ind<- which(neurontypes==Neurontype)
  sub<- cds[,ind]
  expr_sub<- rowMeans(exprs(sub))
  group<-data.frame(sheets=rownames(sub), expr=expr_sub)
  group<- group[order(group$expr, decreasing=TRUE),]
  group<- group[1:10,]
  
  plot<- ggplot(group)+
    geom_col(mapping=aes(x=sheets, y=expr), fill="gold1")+
    xlab("Modules")+
    ylab("Normalized expression")+
    ggtitle(paste0("Neuron type ",Neurontype))+
    theme_classic()
  
  return(plot)}

module_cluster_expression<- function(cds, module_number){
  #The expression of a gene module for all cluster
  M_expression<- exprs(cds[module_number,])
  sheets<-rownames(cds)
  cl<- as.numeric(clusters(cds))
  
  cluster_e<- rep(0, max(cl))
  for (i in 1:max(cl)){
    position<-which(cl==i)
    mean_per_c<- mean(M_expression[,position])
    
    cluster_e[i]<-mean_per_c 
  }
  
  cl2<- as.character(seq(1,max(cl)))
  clustername<- cl2
  
  
  plot<- ggplot()+
    geom_col(mapping=aes(x=clustername, y=cluster_e), fill="gold1")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title.align=0.5)+
    xlab("Cluster")+
    ylab("Normalized expression")+
    ggtitle(paste0("Module ",module_number))
    
  
  return(plot)}

Module_Neuron_expression<- function(cds, module_number){
  #Module expression for all neurons
  M_expression<- exprs(cds[module_number,])
  N<- pData(cds)$Neuron.type
  levels<- levels(as.factor(N))
  
  neuron_e<- rep(0, length(levels))
  for (i in 1:length(levels)){
    position<- which(N==levels[i])
    mean_per_n<- mean(M_expression[,position])
    neuron_e[i]<- mean_per_n
  }
  plot<- ggplot()+
    geom_col(mapping=aes(x=levels, y=neuron_e),fill="gold1")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title.align=0.5, axis.text.y=element_text(size=6))+
    xlab("Neuron type")+
    ylab("Normalized expression")+
    ggtitle(paste0("Module", module_number))
  
  return(plot)
  
}


Module_expressed_in<- function(cds, module_number, expr_thresh){
  #Module expression for all neurons
  M_expression<- exprs(cds[module_number,])
  N<- pData(cds)$Neuron.type
  levels<- levels(as.factor(N))
  
  neuron_e<- rep(0, length(levels))
  for (i in 1:length(levels)){
    position<- which(N==levels[i])
    mean_per_n<- mean(M_expression[,position])
    neuron_e[i]<- mean_per_n
  }
  plot<-levels[which(neuron_e>expr_thresh)]
 
  return(plot)
  
}

