#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

GeneSize <- function(con, geneName, geneNamesType){
  
  if (geneNamesType == "genename") geneNamesType_label <- "Gene Name"
  if (geneNamesType == "emsemblgeneid") geneNamesType_label <- "Ensembl Gene ID"
  

  gene <- dbSendQuery(con, statement = paste0("select * from annotations.reduced_genome where ", geneNamesType, " = '", geneName, "'"))
  gene <- fetch(gene, n = -1)
  
  if(dim(gene)[1] == 0) return(list(0))
  
  if(dim(gene)[1] > 1) return(list(1))
  
  if(dim(gene)[1] == 1){
  #lengthGene <- gene$gene_end - gene$gene_start + 1 
  return(list(gene$gene_start, gene$gene_end))
  }
}


