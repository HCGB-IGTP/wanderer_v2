
#data filtering by zoom

data_meth_filtering <- function(results, sampleSize, tissue, zoom){
  
  tissue_label <- sampleSize[sampleSize[,3]==tissue,1]
  
  ddN<-results[[1]]
  ddT<-results[[2]]
  probes<-results[[3]]
  
  xmin <- zoom[1]
  xmax <- zoom[2]
  
  probes2 <- probes[probes$cg_start >= xmin & probes$cg_start <= xmax,]
  if(!is.null(ddN)){    
    ddN2 <- ddN[ddN$PROBE%in%probes2$probe,]
  } else{
    ddN2 <- NULL
  }
  if(!is.null(ddT)){  
    ddT2 <- ddT[ddT$PROBE%in%probes2$probe,]
  } else{
    ddT2 <- NULL
  }
  
  #if(dim(probes2)[1]>0){
    results_filt <- list(ddN2 = ddN2, ddT2 = ddT2, probes2 = probes2, xmax = xmax, xmin = xmin, tissue_label = tissue_label)
#   } else{
#     results_filt <- list(probes2 = data.frame())
#   }
  return(results_filt)  
  
}
