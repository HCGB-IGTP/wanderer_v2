#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

######################################
#function to plot a region with the probes and the gene

wanderer_methylation <- function(results, geneName, geneNamesType, sampleSize, tissue, zoom, npointsN, npointsT, CpGislands, plotmean, plotting, geneLine){

  options(scipen=20)
  
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
  
  
  if(is.null(zoom)){
    xmin <- gmin
    xmax <- gmax
  }
  if(!is.null(zoom)){
    xmin <- zoom[1]
    xmax <- zoom[2]
  }
  postep <- xmax - xmin
  positions <- round(seq(xmin, xmax, postep/5),0)
  
  
  probes2 <- probes[probes$cg_start >= xmin & probes$cg_start <= xmax,]
  ddN2 <- ddN[ddN$probe%in%probes2$probe,]
  ddT2 <- ddT[ddT$probe%in%probes2$probe,]
  
  if(dim(probes2)[1] == 0) stop(paste0("There are not probes in this region"))
  
  
  if(positions[2]-positions[1]>=100000) escala <- 100000
  if(positions[2]-positions[1]<100000 & positions[2]-positions[1]>=10000) escala <- 10000
  if(positions[2]-positions[1]<10000 & positions[2]-positions[1]>=5000) escala <- 5000
  if(positions[2]-positions[1]<5000 & positions[2]-positions[1]>=1000) escala <- 1000
  if(positions[2]-positions[1]<1000) escala <- 0
  
  
  
  #profile plot
  if(plotting){
    #color to mark CpGislands
    colort <- rep("black", dim(probes2)[1])
    if(CpGislands) colort[!is.na(probes2$cpgiid)] <- "forestgreen"
    
    
    if(!plotmean){
      
      #sample of patients
      if(npointsN <= (dim(ddN2)[2]-1)){
        set.seed(11)
        sampN <- sample(2:dim(ddN2)[2],npointsN)
      } else stop(paste0("The maximum number of Normal samples in ", tissue_label, " for ", data_label, " is ", (dim(ddN2)[2]-1)))
      if(npointsT <= (dim(ddT)[2]-1)){
        set.seed(11)
        sampT <- sample(2:dim(ddT2)[2],npointsT)
      } else stop(paste0("The maximum number of Tumoral samples in ", tissue_label, " for ", data_label, " is ", (dim(ddT2)[2]-1)))
      
      
      ###############
      #plot
      par(mfrow = c(2,1))
      
      #plot for Normals
      par(mai = par()$mai + c(0.3,0,1,0))

      plot(probes2$cg_start, ddN2[,sampN[1]], type = "b", xlim = c(xmin, xmax), 
           pch =  1, cex = 0.7, axes = FALSE, col = "black", ylim = c(-0.2, 1), 
           ylab = "Methylation (beta value)", xlab = "", las = 1)
      for(pl in sampN[-1]) lines(probes2$cg_start, ddN2[,pl], type = "b", pch =  1, cex = 0.7, col = "black")
      title(paste0("Methylation of ", geneName, " in Normal\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes2$probe, at = probes2$cg_start, col.axis = FALSE)
      axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
      axis(3, labels = positions, at = positions, col.axis = FALSE) 
      par(xpd = FALSE)
      mtext(side = 1, text = probes2$probe, at = probes2$cg_start, las = 3, col = colort, line = 1) 
      mtext(side = 3, text = positions, at = positions, las = 1, line = 1, cex = 1) 
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      if(escala!=0){
        par(xpd = TRUE)
        arrows(xmin, 1.25, xmin + escala, 1.25, lwd = 2,  length=0.05, angle=90, code=3)
        text(x = xmin + (escala/2), y = 1.3, labels = paste0((escala/1000),"Kb"))
      }
      par(xpd = FALSE)

      
      #plot for Tumors
      plot(probes2$cg_start, ddT2[,sampT[1]], type = "b", xlim = c(xmin, xmax), 
           pch = 1, cex = 0.7, axes = FALSE, col = "black", ylim = c(-0.2, 1),
           ylab = "Methylation (beta value)", xlab = "", las = 1)
      for(pl in sampT[-1]) lines(probes2$cg_start, ddT2[,pl], type = "b", pch = 1, cex = 0.7, col = "black")
      title(paste0("Methylation of ", geneName, " in Tumor\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes2$probe, at = probes2$cg_start, col.axis = FALSE)
      axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
      axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      par(xpd = FALSE)
      mtext(side = 1, text = probes2$probe, at = probes2$cg_start, las = 3, col = colort, line = 1)
      mtext(side = 3, text = positions, at = positions, las = 1, line = 1, cex = 1) 
      
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)

    }
    
    if(plotmean){
      ###########
      #plot
      par(mfrow = c(2,1))
      
      #plot for Normals
      par(mai = par()$mai + c(0.3,0,1,0))
      boxplot(t(ddN2[,2:dim(ddN)[2]]), at = probes2$cg_start, names = NULL, xlim = c(xmin, xmax), 
              pch = 20, cex = 0.5, axes = FALSE, ylim = c(-0.2, 1), 
              ylab = "Methylation (beta value)", xlab = "", las = 1, boxwex = (xmax-xmin)*0.03, varwidth = FALSE)
      title(paste0("Methylation of ", geneName, " in Normal\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes2$probe, at = probes2$cg_start, col.axis = FALSE)
      axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
      axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      par(xpd = FALSE)
      mtext(side = 1, text = probes2$probe, at = probes2$cg_start, las = 3, col = colort, line = 1)
      mtext(side = 3, text = positions, at = positions, las = 1, line = 1, cex = 1) 
      
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      par(xpd = TRUE)
      arrows(xmin, 1.25, xmin + 10000, 1.25, lwd = 2,  length=0.05, angle=90, code=3)
      text(x = xmin + 5000, y = 1.3, labels = "10Kb")
      par(xpd = FALSE)
      
      #plot for Tumors
      boxplot(t(ddT2[,2:dim(ddT)[2]]), at = probes2$cg_start, names = NULL, xlim = c(xmin, xmax), 
              pch = 20, cex = 0.5, axes = FALSE, ylim = c(-0.2, 1), 
              ylab = "Methylation (beta value)", xlab = "", las = 1, boxwex = (xmax-xmin)*0.03, varwidth = FALSE)
      title(paste0("Methylation of ", geneName, " in Tumor\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes2$probe, at = probes2$cg_start, col.axis = FALSE)
      axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
      axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      par(xpd = FALSE)
      mtext(side = 1, text = probes2$probe, at = probes2$cg_start, las = 3, col = colort, line = 1)
      mtext(side = 3, text = positions, at = positions, las = 1, line = 1, cex = 1) 
      
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
    }
  }
  
  plotting_results <- list(ddN2 = ddN2, ddT2 = ddT2, probes2 = probes2)
  
  options(scipen=0)
  
  return(plotting_results)
}



