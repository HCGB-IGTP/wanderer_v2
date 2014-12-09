#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

######################################
#function to query expression data of a gene from db

expression_data <- function(con, geneName, geneNamesType, tissue){
  
  
  if (geneNamesType == "genename") geneNamesType_label <- "Gene Name"
  if (geneNamesType == "emsemblgeneid") geneNamesType_label <- "Ensembl Gene ID"
  
  ################
  #RNAseq annotation download
  probes <- dbSendQuery(con, statement = paste0("select * from annotations.exons_annot
            where chr = (select chr from annotations.exons_annot where ", geneNamesType, " = '", geneName, "' limit 1)
            and exon_start > (select genestart as sstart from annotations.exons_annot
            where ", geneNamesType, " = '", geneName, "' limit 1)
            and exon_end < (select geneend as send from annotations.exons_annot
            where ", geneNamesType, " = '", geneName, "' limit 1);"))
  
  probes <- fetch(probes, n = -1)
  if(dim(probes)[1] == 0) stop(paste0("There are not exons in this region"))
  if(dim(probes)[1] != 0){
    ordering <- order(probes$exon_start)
    probes <- probes[ordering,]
    
    ##############
    probes_collapsed <- paste0("('", paste(probes$exon, collapse="','"), "')")
    
    #methylation array tumor data download
    ddT <- dbSendQuery(con, paste0("select * from ", tissue, "_tumor.illuminahiseq_rnaseqv2
          where exon in ", probes_collapsed, ";"))
    ddT <- fetch(ddT, n = -1)
    ddT <- ddT[ordering,]
    ddT[,2:dim(ddT)[2]] <- log2(ddT[,2:dim(ddT)[2]]+1)
    
    #methylation array normal data download
    ddN <- dbSendQuery(con, paste0("select * from ", tissue, "_normal.illuminahiseq_rnaseqv2
            where exon in ", probes_collapsed, ";"))
    ddN <- fetch(ddN, n = -1)
    ddN <- ddN[ordering,]
    ddN[,2:dim(ddN)[2]] <- log2(ddN[,2:dim(ddN)[2]]+1)
    
    #removing NA's
    naT <- which(apply(is.na(ddT), 1, sum) == (dim(ddT)[2] - 1))
    if(length(naT)>0){
      probes <- probes[-naT,]
      ddN <- ddN[-naT,]
      ddT <- ddT[-naT,]
    }
    
    naN <- which(apply(is.na(ddN), 1, sum) == (dim(ddN)[2] - 1))
    if(length(naN)>0){
      probes <- probes[-naN,]
      ddT <- ddT[-naN,]
      ddN <- ddN[-naN,]
    }
    if(dim(probes)[1] == 0) stop(paste0("There are not exons in this region"))
    
  }
  
  results <- list(ddN, ddT, probes)
  return(results)
}


