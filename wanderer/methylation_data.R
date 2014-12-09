#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

######################################
#function to query methylation data of a gene from db

methylation_data <- function(con, geneName, geneNamesType, tissue){
  
  
  
  if (geneNamesType == "genename") geneNamesType_label <- "Gene Name"
  if (geneNamesType == "emsemblgeneid") geneNamesType_label <- "Ensembl Gene ID"
  
  
  ################
  #methylation array annotation download
  probes <- dbSendQuery(con, statement = paste0("select * from annotations.humanmethylation450kannot
            where chr = (select chr from annotations.humanmethylation450kannot where ", geneNamesType, " = '", geneName, "' limit 1)
            and probestart > (select genestart-", 50000, " as sstart from annotations.humanmethylation450kannot
            where ", geneNamesType, " = '", geneName, "' limit 1)
            and probeend < (select geneend+", 50000, " as send from annotations.humanmethylation450kannot
            where ", geneNamesType, " = '", geneName, "' limit 1);"))
  
  probes <- fetch(probes, n = -1)
  print(dim(probes))
  if(dim(probes)[1] == 0) stop(paste0("There are not probes in this region"))
  if(dim(probes)[1] != 0){
    probes <- probes[order(probes$probe),]
    ordering <- order(probes$probestart)
    probes <- probes[ordering,]
    ##############
    probes_collapsed <- paste0("('", paste(probes$probe, collapse="','"), "')")
    
    #methylation array tumor data download
    ddT <- dbSendQuery(con, paste0("select * from ", tissue, "_tumor.humanmethylation450
          where probe in ", probes_collapsed, ";"))
    ddT <- fetch(ddT, n = -1)
    print(dim(ddT))
    ddT <- ddT[order(ddT$probe),]
    ddT <- ddT[ordering,]
    
    #methylation array normal data download
    ddN <- dbSendQuery(con, paste0("select * from ", tissue, "_normal.humanmethylation450
          where probe in ", probes_collapsed, ";"))
    ddN <- fetch(ddN, n = -1)
    print(dim(ddN))
    ddN <- ddN[order(ddN$probe),]
    ddN <- ddN[ordering,]
    
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
  }
  if(dim(probes)[1] == 0) stop(paste0("There are not probes in this region"))
  
  results <- list(ddN, ddT, probes)
  return(results)
}