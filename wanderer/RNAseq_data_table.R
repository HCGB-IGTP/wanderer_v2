
RNAseq_data_table <- function(ddGene, ddExons){
  
  order_columns <- match(colnames(ddExons)[-1], colnames(ddGene)[-1])
  nams<-colnames(ddGene)[-1]
  nams<-nams[order_columns]
  dataGene <- as.data.frame(ddGene[,-1])
  dataGene <- data.frame(ddGene[,1],dataGene[,order_columns])
  colnames(dataGene)<-c("GENE",nams)
  for(i in 1:6) colnames(dataGene) <- sub("\\.","-",colnames(dataGene))
  dataGene <- t(dataGene)
  dataGene <- data.frame(rownames(dataGene), dataGene[,1])
  colnames(dataGene) <- c("PATIENT",paste0(dataGene[1,2]))
  dataGene <- dataGene[-1,]
  return(dataGene)
}

