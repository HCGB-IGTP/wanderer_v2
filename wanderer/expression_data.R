#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

######################################
#function to query expression data of a gene from db

expression_data <- function(con, geneName, geneNamesType, tissue){
  
  if (geneNamesType == "genename") geneNamesType_label <- "Gene Name"
  if (geneNamesType == "emsemblgeneid") geneNamesType_label <- "Ensembl Gene ID"
  
  exonsaux <- dbSendQuery(con, statement = paste0("select * from annotations.exons_annot where ", geneNamesType, " = '", geneName, "' limit 1"))
  exonsaux <- fetch(exonsaux, n = -1)
  
  ############################################
  #RNAseq annotation download
  if(dim(exonsaux)[1] > 0){
    
    exons <- dbSendQuery(con, statement = paste0("select * from annotations.exons_annot where ", geneNamesType, " = '", geneName,"';"))
    exons <- fetch(exons, n = -1)
    if(dim(exons)[1]>0){
      ordering <- order(exons$exon_start)
      exons <- exons[ordering,]
      
      ##############
      exons_collapsed <- paste0("('", paste(exons$exon, collapse="','"), "')")
      
      #methylation array tumor data download
      ddT <- try(dbSendQuery(con, paste0("select * from ", tissue, "_tumor.illuminahiseq_rnaseqv2
          where exon in ", exons_collapsed, ";")),silent=TRUE)
      if(class(ddT)=="try-error"){
        ddT<-NULL
      } else{
        ddT <- fetch(ddT, n = -1)
        ddT <- ddT[ordering,]
        ddT[,2:dim(ddT)[2]] <- log2(ddT[,2:dim(ddT)[2]]+1)
        colnames(ddT) <- toupper(colnames(ddT))
        for(i in 1:6) colnames(ddT) <- sub("_","-",colnames(ddT))
      }
      
      
      #methylation array normal data download
      ddN <- try(dbSendQuery(con, paste0("select * from ", tissue, "_normal.illuminahiseq_rnaseqv2
            where exon in ", exons_collapsed, ";")),silent=TRUE)
      if(class(ddN)=="try-error"){
        ddN<-NULL
      } else{
        ddN <- fetch(ddN, n = -1)
        ddN <- ddN[ordering,]
        ddN[,2:dim(ddN)[2]] <- log2(ddN[,2:dim(ddN)[2]]+1)
        colnames(ddN) <- toupper(colnames(ddN))
        for(i in 1:6) colnames(ddN) <- sub("_","-",colnames(ddN))
      }
      
      #removing NA's
      if(!is.null(ddT)){
        naT <- which(apply(is.na(ddT), 1, sum) == (dim(ddT)[2] - 1))
        if(length(naT)>0){
          exons <- exons[-naT,]
          ddN <- ddN[-naT,]
          ddT <- ddT[-naT,]
        }
      }
      
      if(!is.null(ddN)){
        naN <- which(apply(is.na(ddN), 1, sum) == (dim(ddN)[2] - 1))
        if(length(naN)>0){
          exons <- exons[-naN,]
          ddT <- ddT[-naN,]
          ddN <- ddN[-naN,]
        }
      }
      
      results <- list(ddN, ddT, exons, empty = FALSE, geneNamesType_label)
    } else{
      results <- list(empty = TRUE, geneNamesType_label)
      stop(paste0("There are not exons in this region"))
    }
  } else{
    results <- list(empty = TRUE, geneNamesType_label)
    stop(paste0("There are not exons in this region"))
  }
  
  return(results)
}