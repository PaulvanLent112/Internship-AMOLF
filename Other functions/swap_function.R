swappy = function(v,a,b){  # where v is a dataframe, a and b are the columns indexes to swap
  
  name = deparse(substitute(v))
  
  helpy = v[,a]
  v[,a] = v[,b]
  v[,b] = helpy
  
  
  name1 = colnames(v)[a] 
  name2 = colnames(v)[b] 
  
  colnames(v)[a] = name2
  colnames(v)[b] = name1
  
  assign(name,value = v , envir =.GlobalEnv)
}


