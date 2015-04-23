
data_RNAseqGene <- function(con, dataexpr, tissue){
  
  dd <- dataexpr[['exons2']]
  geneid <- unique(dd$rnaseqgeneid)
    
  if(!is.na(geneid)){
    
    #RNAseq at gene level tumor data download
    ddT <- try(dbSendQuery(con, paste0("select * from ", tissue, "_tumor.illuminahiseq_rnaseqv2_by_gene
                                     where gene = '", geneid, "';")), silent=TRUE)
    if(class(ddT)=="try-error"){
      ddT<-NULL
    } else{
      ddT <- fetch(ddT, n = -1)
      colnames(ddT) <- toupper(colnames(ddT))
      for(i in 1:6) colnames(ddT) <- sub("_","-",colnames(ddT))
    }
    
    
    #RNAseq at gene level normal data download
    ddN <- try(dbSendQuery(con, paste0("select * from ", tissue, "_normal.illuminahiseq_rnaseqv2_by_gene
                                     where gene = '", geneid, "';")), silent=TRUE)
    if(class(ddN)=="try-error"){
      ddN<-NULL
    } else{
      ddN <- fetch(ddN, n = -1)
      colnames(ddN) <- toupper(colnames(ddN))
      for(i in 1:6) colnames(ddN) <- sub("_","-",colnames(ddN))
    }
    results <- list(ddN, ddT)
    names(results) <- c("Normal", "Tumor")
    
  } else {
    results <- NULL
  }
  return(results)
  
}


