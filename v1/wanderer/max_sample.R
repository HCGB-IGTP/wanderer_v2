
max_sample <- function(sample_size, DataType, tissue){
  
  sample_size<-sample_size[sample_size[,3]==tissue,]
  
  if(DataType=="methylation") sample_size <- as.numeric(sample_size[,4:5])
  if(DataType=="expression") sample_size <- as.numeric(sample_size[,6:7])
  
  return(sample_size)
}
  