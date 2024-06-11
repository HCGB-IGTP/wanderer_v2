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


##' Extracts the genesymbol of a ensembl id (ENSG identifier)
##'
##' Based on static, local-db pairing
##' @title 
##' @param con open connection to the database
##' @param ensembl_id the ENSG identifier
##' @return a string with the genesymbol
##' @author Izaskun Mallona
ensembl_to_gene_symbol <- function(con, ensembl_id) {
    stmt <- "SELECT genename FROM annotations.reduced_genome WHERE emsemblgeneid = '%s'"
    gene_symbol <- fetch(dbSendQuery(con, statement = sprintf(stmt, ensembl_id)), n = -1)

    return(as.character(gene_symbol))                              
}
