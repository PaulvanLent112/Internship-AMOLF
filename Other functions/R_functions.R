feature_vector<- function(file_name){
  require(readxl)
  list_of_sheets<- excel_sheets(file_name)
  feature_vector<- c()
  for (i in 1:length(list_of_sheets)){
    sheet<-read_xlsx(file_name, sheet = i, col_names=FALSE)
    feature_vector<- c(feature_vector, sheet)}
  return(feature_vector)
  }


vec_list<- function(feature_vector){
  #This stores all genes in 1 vector, so that we can easily subset the CDS dataset
  vec_list<- c()
  for (i in 1:length(feature_vector)){
    vec_list<- c(vec_list, feature_vector[i]$...1)
  }
  return(vec_list)
}

index<- function(feature_vector){
  #Gives us the integers that represent each module of the vec_list
  index<- list()
  count<-1
  for (i in 1:length(feature_vector)){
    module_index<-  seq(from=count, to=count-1+length(feature_vector[[i]]))
    index<- append(index, list(module_index))
    count<- count+ length(feature_vector[[i]])
  }
  return(index)
}



