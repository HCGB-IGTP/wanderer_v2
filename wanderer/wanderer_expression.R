#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

######################################
#function to plot a region with the exons and the gene


wanderer_expression <- function(results, geneName, geneNamesType, sampleSize, tissue, zoom, walk, npointsN, npointsT, plotmean, plotting, geneLine){
  
  tissue_label <- sampleSize[sampleSize[,3]==tissue,1]
  data_label <- "Illumina HiSeq RNAseq"
  
  ddN<-results[[1]]
  ddT<-results[[2]]
  probes<-results[[3]]
  
  #axis limits
  gmin <- unique(probes$genestart[probes[,paste0(geneNamesType)] == geneName])
  gmax <- unique(probes$geneend[probes[,paste0(geneNamesType)] == geneName])
  gstrand <- unique(probes$strand[probes[,paste0(geneNamesType)] == geneName])
  gchr <- unique(probes$chr[probes[,paste0(geneNamesType)] == geneName])
  
  
  if(is.null(walk)) walk = 0 
  if(is.null(zoom)) zoom = 0 
  xmin <- gmin + zoom + walk
  xmax <- gmax - zoom + walk
  
  probes2 <- probes[probes$exon_start >= xmin & probes$exon_start <= xmax,]
  ddN2 <- ddN[ddN$exon%in%probes2$exon,]
  ddT2 <- ddT[ddT$exon%in%probes2$exon,]
  
  if(plotting){
    if(!plotmean){
      #sample of patients
      if(npointsN <= (dim(ddN)[2]-1)){
        set.seed(11)
        sampN <- sample(2:dim(ddN)[2],npointsN)
      } else stop(paste0("The maximum number of Normal samples in ", tissue_label, " for ", data_label, " is ", (dim(ddN)[2]-1)))
      if(npointsT <= (dim(ddT)[2]-1)){
        set.seed(11)
        sampT <- sample(2:dim(ddT)[2],npointsT)
      } else stop(paste0("The maximum number of Tumoral samples in ", tissue_label, " for ", data_label, " is ", (dim(ddT)[2]-1)))
      
      
      #y axis  
      ymax <- max(ddT[,-1], ddN[,-1])
      
      
      ###############
      #plot
      par(mfrow = c(2,1))
      
      #plot for Normals
      par(mai = par()$mai + c(1.8,0,0,0))
      plot(probes$exon_start, ddN[,sampN[1]], type = "b", xlim = c(xmin, xmax), 
           pch = 1, cex = 0.7, xaxt = "n", col = "black", ylim = c(-0.2, ymax), 
           ylab = "log2(rpkm + 1)", xlab = "", las = 1)
      for(pl in sampN[-1]) lines(probes$exon_start, ddN[,pl], type = "b", pch = 1, cex = 0.7, col = "black")
      title(paste0("Expression of ", geneName, " ", "Normal ", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$exon, at = probes$exon_start, las=3)
      if(geneLine){
        if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      #plot for Tumors
      plot(probes$exon_start, ddT[,sampT[1]], type = "b", xlim = c(xmin, xmax), 
           pch = 1, cex = 0.7, xaxt = "n", col = "black", ylim = c(-0.2, ymax),
           ylab = "log2(rpkm + 1)", xlab = "", las = 1)
      for(pl in sampT[-1]) lines(probes$exon_start, ddT[,pl], type = "b", pch = 1, cex = 0.7, col = "black")
      title(paste0("Expression of ", geneName, " ", "Tumor ", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$exon, at = probes$exon_start, las=3) 
      if(geneLine){
        if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
    }
    
    if(plotmean){
      #y axis  
      ymax <- max(ddT[,-1], ddN[,-1])
      
      
      ###############
      #plot
      par(mfrow = c(2,1))
      
      #plot for Normals
      par(mai = par()$mai + c(2,0,0,0))
      boxplot(t(ddN[,2:dim(ddN)[2]]), at = probes$exon_start, names = NULL, xlim = c(xmin, xmax), 
              pch = 20, cex = 0.5, xaxt = "n", ylim = c(-0.2, ymax), 
              ylab = "log2(rpkm + 1)", xlab = "", las = 1, boxwex = 2000, varwidth = FALSE)
      title(paste0("Expression of ", geneName, " ", "Normal ", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$exon, at = probes$exon_start, las=3)
      if(geneLine){
        if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      #plot for Tumors
      boxplot(t(ddT[,2:dim(ddT)[2]]), at = probes$exon_start, names = NULL, xlim = c(xmin, xmax), 
              pch = 20, cex = 0.5, xaxt = "n", ylim = c(-0.2, ymax), 
              ylab = "log2(rpkm + 1)", xlab = "", las = 1, boxwex = 2000, varwidth = FALSE)#, border = "#dd4814")
      title(paste0("Expression of ", geneName, " ", "Tumor ", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$exon, at = probes$exon_start, las=3)
      if(geneLine){
        if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
    }
  }
  
  plotting_results <- list(ddN2 = ddN2, ddT2 = ddT2, probes2 = probes2)
  
  return(plotting_results) 
}
