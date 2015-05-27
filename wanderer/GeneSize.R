#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

GeneSize <- function(con, geneName, geneNamesType){
  
  gene <- dbSendQuery(con, statement = paste0("select * from annotations.reduced_genome where ", geneNamesType, " = '", geneName, "'"))
  gene <- fetch(gene, n = -1)
  
  if(dim(gene)[1] == 0) return(list(0))
  if(dim(gene)[1] > 1) return(list(1))
  
  if(dim(gene)[1] == 1) return(list(gene$gene_start, gene$gene_end, gene$chrom))
}


