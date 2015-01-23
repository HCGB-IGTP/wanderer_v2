
#data filtering by zoom

data_expr_filtering <- function(results, sampleSize, tissue, zoom){
  
  tissue_label <- sampleSize[sampleSize[,3]==tissue,1]
  
  ddN<-results[[1]]
  ddT<-results[[2]]
  exons<-results[[3]]
  
  xmin <- zoom[1]
  xmax <- zoom[2]
  
  exons2 <- exons[exons$exon_start+1 >= xmin & exons$exon_start <= xmax,]
  ddN2 <- ddN[ddN$EXON%in%exons2$exon,]
  ddT2 <- ddT[ddT$EXON%in%exons2$exon,]
  
  if(dim(exons2)[1]>0){
    results_filt <- list(ddN2 = ddN2, ddT2 = ddT2, exons2 = exons2, xmax = xmax, xmin = xmin, tissue_label = tissue_label)
  } else{
    results_filt <- list(exons2 = data.frame())
  }
  
  return(results_filt)  
  
}
