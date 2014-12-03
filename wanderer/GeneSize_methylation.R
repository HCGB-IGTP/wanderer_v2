
GeneSize_methylation <- function(con, geneName, geneNamesType){
  
  if (geneNamesType == "genename") geneNamesType_label <- "Gene Name"
  if (geneNamesType == "ensemblgeneid") geneNamesType_label <- "Ensembl Gene ID"
  

  ################
  #methylation array annotation download
  probesaux <- dbSendQuery(con, statement = paste0("select * from annotations.humanmethylation450kannot where ", geneNamesType, " = '", geneName, "' limit 1"))
  probesaux <- fetch(probesaux, n = -1)
  if(dim(probesaux)[1] == 0) stop(paste0("Your gene ", geneName, " does not correspon to any ", geneNamesType_label, " or is not in the probes annotation"))
  
  lengthGene <- probesaux$geneend - probesaux$genestart + 1 
  return(lengthGene)
}


