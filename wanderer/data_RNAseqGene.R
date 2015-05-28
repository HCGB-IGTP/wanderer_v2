
data_RNAseqGene <- function(con, geneNamesType, geneName, tissue){
  
  exonsaux <- dbSendQuery(con, statement = paste0("select * from annotations.exons_annot where ", geneNamesType, " = '", geneName, "'"))
  exonsaux <- fetch(exonsaux, n = -1)
  exonsaux <- exonsaux[,-c(1,2,4,6,7)]
  exonsaux <- exonsaux[!duplicated(exonsaux),]
  geneid <- unique(exonsaux$rnaseqgeneid)
  
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
    ddT[,-1] <- log2(ddT[,-1] + 1)
    
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
    ddN[,-1] <- log2(ddN[,-1] + 1)
    
    
    results <- list(ddN, ddT)
    names(results) <- c("Normal", "Tumor")
    
  } else {
    results <- NULL
  }
  return(results)
  
  }

