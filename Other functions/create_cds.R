create_cds<- function(output, filename){
  filename<- filename
  source("functions/other_functions/R_functions.R") 
  feature_vector<-feature_vector(filename) #loads vector of excel sheet indexes of genes
  list_of_sheets<- excel_sheets(filename)
  vec_list<- vec_list(feature_vector) #all genes from the excel sheet
  index<- index(feature_vector) #Indices of the genes for subsetting the Average Gene Module Dataframe.
  
  norm<- output$Norm[vec_list,]
  
  AGM_expression<- matrix(NaN ,nrow=length(list_of_sheets), ncol=ncol(norm))
  for (i in 1: length(index)){
    sub<-norm[index[[i]],]
    avg<- colMeans(sub)
    AGM_expression[i,]<-avg
  }
  rownames(AGM_expression)<- list_of_sheets
  
  expr_matrix<- AGM_expression
  colnames(expr_matrix)<- colnames(output$neurons2)
  rownames(expr_matrix)<- list_of_sheets
  
  size=c()
  for (i in 1:length(feature_vector)){
    size=c(size, length(index[[i]]))
  }
  
  module_data<- data.frame(size=size)
  rownames(module_data)<-rownames(expr_matrix)
  module_data$gene_short_name<- rownames(module_data)
  
  sample_data<- data.frame(Sample=output$neurons2$Sample)
  rownames(sample_data)<- output$neurons2$Barcode
  sample_data$Neuron.type<- output$neurons2$Neuron.type
  sample_data$Experiment<- output$neurons2$Experiment
  library(monocle3)
  cds <- monocle3::new_cell_data_set(expr_matrix,
                                     cell_metadata= sample_data, gene_metadata=module_data) 
  
  return(cds)
}
