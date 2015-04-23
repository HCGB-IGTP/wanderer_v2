
RNAseq_data_table <- function(ddGene, ddExons){
  
  order_columns <- match(colnames(ddExons)[-1], colnames(ddGene)[-1])
  dataGene <- ddGene[,-1]
  dataGene <- data.frame(GENE=ddGene[,1],dataGene[,order_columns])
  for(i in 1:6) colnames(dataGene) <- sub("\\.","-",colnames(dataGene))
  dataGene <- t(dataGene)
  dataGene <- data.frame(rownames(dataGene), dataGene[,1])
  colnames(dataGene) <- c("PATIENT",paste0(dataGene[1,2]))
  dataGene <- dataGene[-1,]
  return(dataGene)
}





