#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

######################################
#function to plot a region with the probes and the gene

wanderer_methylation <- function(results, geneName, geneNamesType, sampleSize, tissue, zoom, walk, npointsN, npointsT, CpGislands, plotmean, plotting, geneLine){

  
  tissue_label <- sampleSize[sampleSize[,3]==tissue,1]
  data_label <- "450k Methylation Array"
  
  
  ddN<-results[[1]]
  ddT<-results[[2]]
  probes<-results[[3]]
  
  #axis limits
  gmin <- unique(probes$genestart[probes[,paste0(geneNamesType)] == geneName])
  gmax <- unique(probes$geneend[probes[,paste0(geneNamesType)] == geneName])
  gstrand <- unique(probes$genestrand[probes[,paste0(geneNamesType)] == geneName])
  gchr <- unique(probes$chr[probes[,paste0(geneNamesType)] == geneName])
  
  if(is.null(walk)) walk = 0 
  if(is.null(zoom)) zoom = 0 
  xmin <- gmin + zoom + walk
  xmax <- gmax - zoom + walk
  postep <- xmax - xmin
  positions <- round(seq(xmin, xmax, postep/5),0)
  
  probes2 <- probes[probes$cg_start >= xmin & probes$cg_start <= xmax,]
  ddN2 <- ddN[ddN$probe%in%probes2$probe,]
  ddT2 <- ddT[ddT$probe%in%probes2$probe,]
  
  
  
  #profile plot
  if(plotting){
    #color to mark CpGislands
    colort <- rep("black", dim(probes)[1])
    if(CpGislands) colort[!is.na(probes$cpgiid)] <- "forestgreen"
    
    
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
      
      
      ###############
      #plot
      par(mfrow = c(2,1))
      
      #plot for Normals
      par(mai = par()$mai + c(0.3,0,1,0))

      plot(probes$cg_start, ddN[,sampN[1]], type = "b", xlim = c(xmin, xmax), 
           pch =  1, cex = 0.7, xaxt = "n", col = "black", ylim = c(-0.2, 1), 
           ylab = "beta value", xlab = "", las = 1)
      for(pl in sampN[-1]) lines(probes$cg_start, ddN[,pl], type = "b", pch =  1, cex = 0.7, col = "black")
      title(paste0("Methylation of ", geneName, " ", "Normal ", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$probe, at = probes$cg_start, col.axis = FALSE) 
      axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      mtext(side = 1, text = probes$probe, at = probes$cg_start, las = 3, col = colort, line = 1) 
      mtext(side = 3, text = positions, at = positions, las = 1, line = 1, cex = 0.7) 
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)

      
      #plot for Tumors
      plot(probes$cg_start, ddT[,sampT[1]], type = "b", xlim = c(xmin, xmax), 
           pch = 1, cex = 0.7, xaxt = "n", col = "black", ylim = c(-0.2, 1),
           ylab = "beta value", xlab = "", las = 1)
      for(pl in sampT[-1]) lines(probes$cg_start, ddT[,pl], type = "b", pch = 1, cex = 0.7, col = "black")
      title(paste0("Methylation of ", geneName, " ", "Tumor ", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$probe, at = probes$cg_start, col.axis = FALSE) 
      axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      mtext(side = 1, text = probes$probe, at = probes$cg_start, las = 3, col = colort, line = 1)
      mtext(side = 3, text = positions, at = positions, las = 1, line = 1, cex = 0.7) 
      
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      par(xpd = TRUE)
      arrows(xmax - 10000, -1, xmax, -1, lwd = 2,  length=0.05, angle=90, code=3)
      text(x = xmax - 5000, y = -0.9, labels = "10Kb")
      par(xpd = FALSE)
    }
    
    if(plotmean){
      ###########
      #plot
      par(mfrow = c(2,1))
      
      #plot for Normals
      par(mai = par()$mai + c(0.3,0,1,0))
      boxplot(t(ddN[,2:dim(ddN)[2]]), at = probes$cg_start, names = NULL, xlim = c(xmin, xmax), 
              pch = 20, cex = 0.5, xaxt = "n", ylim = c(-0.2, 1), 
              ylab = "beta value", xlab = "", las = 1, boxwex = 2000, varwidth = FALSE)
      title(paste0("Methylation of ", geneName, " ", "Normal ", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$probe, at = probes$cg_start, col.axis = FALSE) 
      axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      mtext(side = 1, text = probes$probe, at = probes$cg_start, las = 3, col = colort, line = 1)
      mtext(side = 3, text = positions, at = positions, las = 1, line = 1, cex = 0.7) 
      
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      #plot for Tumors
      boxplot(t(ddT[,2:dim(ddT)[2]]), at = probes$cg_start, names = NULL, xlim = c(xmin, xmax), 
              pch = 20, cex = 0.5, xaxt = "n", ylim = c(-0.2, 1), 
              ylab = "beta value", xlab = "", las = 1, boxwex = 2000, varwidth = FALSE)
      title(paste0("Methylation of ", geneName, " ", "Tumor ", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$probe, at = probes$cg_start, col.axis = FALSE)
      axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      mtext(side = 1, text = probes$probe, at = probes$cg_start, las = 3, col = colort, line = 1)
      mtext(side = 3, text = positions, at = positions, las = 1, line = 1, cex = 0.7) 
      
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      par(xpd = TRUE)
      arrows(xmax - 10000, -1, xmax, -1, lwd = 2,  length=0.05, angle=90, code=3)
      text(x = xmax - 5000, y = -0.9, labels = "10Kb")
      par(xpd = FALSE)
    }
  }
  results <- list(ddN2, ddT2, probes2)
  return(results)
  
}



