#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

######################################
#function to plot a region with the exons and the gene


wanderer_expression <- function(results_filt, geneName, geneNamesType, npointsN, npointsT, plotmean, plotting, geneLine){
  
  options(scipen=20)
  
  #data_label <- "Illumina HiSeq RNAseq"
  
  
  ddN2 <- results_filt$ddN2
  ddT2 <- results_filt$ddT2
  exons2 <- results_filt$exons2
  tissue_label <- results_filt$tissue_label
  xmin <- results_filt$xmin
  xmax <- results_filt$xmax
  
  #axis limits
  gmin <- unique(exons2$genestart[exons2[,paste0(geneNamesType)] == geneName])
  gmax <- unique(exons2$geneend[exons2[,paste0(geneNamesType)] == geneName])
  gstrand <- unique(exons2$strand[exons2[,paste0(geneNamesType)] == geneName])
  gchr <- unique(exons2$chr[exons2[,paste0(geneNamesType)] == geneName])
  
  
  posmin <- ((xmin%/%1000)-1)*1000
  posmax <- ((xmax%/%1000)+1)*1000
  postep <- posmax - posmin
  positions <- round(seq(posmin, posmax, postep/5),0)
  positions <- (positions%/%1000)*1000
  
  
  
  if(dim(exons2)[1] > 0){
    
    if(positions[2]-positions[1]>=100000) escala <- 100000
    if(positions[2]-positions[1]<100000 & positions[2]-positions[1]>=10000) escala <- 10000
    if(positions[2]-positions[1]<10000 & positions[2]-positions[1]>=5000) escala <- 5000
    if(positions[2]-positions[1]<5000 & positions[2]-positions[1]>=1000) escala <- 1000
    if(positions[2]-positions[1]<1000) escala <- 0
    
    
    sampN <- c(2:dim(ddN2)[2])
    sampT <- c(2:dim(ddT2)[2])
    
    
    if(plotting){
      if(!plotmean){
        #sample of patients
        if(npointsN < (dim(ddN2)[2]-1)){
          set.seed(1234)
          sampN <- sample(x=c(2:dim(ddN2)[2]),size=npointsN)
        } 
        
        if(npointsT < (dim(ddT2)[2]-1)){
          set.seed(1234)
          sampT <- sample(x=c(2:dim(ddT2)[2]),size=npointsT)
        } 
        
        #y axis  
        ymax <- max(ddT2[,-1], ddN2[,-1])
        
        
        ###############
        #plot
        par(mfrow = c(2,1))
        
        #plot for Normals
        par(mai = par()$mai + c(2,0,1,0))
        
        plot(exons2$exon_start, ddN2[,sampN[1]], type = "b", xlim = c(xmin, xmax), 
             pch = 1, cex = 0.7, axes = FALSE, col = "black", ylim = c(-0.2, ymax), 
             ylab = "Expression log2(rpkm + 1)", xlab = "", las = 1)
        for(pl in sampN[-1]) lines(exons2$exon_start, ddN2[,pl], type = "b", pch = 1, cex = 0.7, col = "black")
        title(paste0("Expression of ", geneName, " in Normal (n=", length(sampN), ")\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
        axis(1, labels = exons2$exon, at = exons2$exon_start, col.axis = FALSE)
        axis(2, labels = seq(0, ymax, 1), at = seq(0, ymax, 1), las = 1)
        axis(3, labels = positions, at = positions, col.axis = FALSE) 
        
        par(xpd = FALSE)
        mtext(side = 1, text = exons2$exon, at = exons2$exon_start, las = 3, line = 1) 
        mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
        
        
        if(geneLine){
          if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
          if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        }
        box(lwd = 1.5)
        
        if(escala!=0){
          par(xpd = TRUE)
          arrows(xmin, ymax + (ymax/3), xmin + escala, ymax + (ymax/3), lwd = 2,  length=0.05, angle=90, code=3)
          text(x = xmin + (escala/2), y = ymax + (ymax/3) + 0.5, labels = paste0((escala/1000),"Kb"))
        }
        par(xpd = FALSE)
        
        #plot for Tumors
        plot(exons2$exon_start, ddT2[,sampT[1]], type = "b", xlim = c(xmin, xmax), 
             pch = 1, cex = 0.7, axes = FALSE, col = "black", ylim = c(-0.2, ymax),
             ylab = "Expression log2(rpkm + 1)", xlab = "", las = 1)
        for(pl in sampT[-1]) lines(exons2$exon_start, ddT2[,pl], type = "b", pch = 1, cex = 0.7, col = "black")
        title(paste0("Expression of ", geneName, " in Tumor (n=", length(sampT), ")\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
        axis(1, labels = exons2$exon, at = exons2$exon_start, col.axis = FALSE)
        axis(2, labels = seq(0, ymax, 1), at = seq(0, ymax, 1), las = 1)
        axis(3, labels = positions, at = positions, col.axis = FALSE) 
        
        par(xpd = FALSE)
        mtext(side = 1, text = exons2$exon, at = exons2$exon_start, las = 3, line = 1) 
        mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
        
        if(geneLine){
          if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
          if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        }
        box(lwd = 1.5)
      }
      
      
      if(plotmean){
        #y axis  
        ymax <- max(ddT2[,-1], ddN2[,-1])
        
        
        ###############
        #plot
        par(mfrow = c(2,1))
        
        #plot for Normals
        par(mai = par()$mai + c(2,0,1,0))
        boxplot(t(ddN2[,2:dim(ddN2)[2]]), at = exons2$exon_start, names = NULL, xlim = c(xmin, xmax), 
                pch = 20, cex = 0.5, axes = FALSE, ylim = c(-0.2, ymax), 
                ylab = "Expression log2(rpkm + 1)", xlab = "", las = 1, boxwex = (xmax-xmin)*0.03, varwidth = FALSE)
        title(paste0("Expression of ", geneName, " in Normal (n=", (dim(ddN2)[2]-1), ")\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
        axis(1, labels = exons2$exon, at = exons2$exon_start, col.axis = FALSE)
        axis(2, labels = seq(0, ymax, 1), at = seq(0, ymax, 1), las = 1)
        axis(3, labels = positions, at = positions, col.axis = FALSE) 
        
        par(xpd = FALSE)
        mtext(side = 1, text = exons2$exon, at = exons2$exon_start, las = 3, line = 1) 
        mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
        
        
        if(geneLine){
          if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
          if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        }
        box(lwd = 1.5)
        
        if(escala!=0){
          par(xpd = TRUE)
          arrows(xmin, ymax + (ymax/3), xmin + escala, ymax + (ymax/3), lwd = 2,  length=0.05, angle=90, code=3)
          text(x = xmin + (escala/2), y = ymax + (ymax/3) + 0.5, labels = paste0((escala/1000),"Kb"))
        }
        par(xpd = FALSE)
        
        #plot for Tumors
        boxplot(t(ddT2[,2:dim(ddT2)[2]]), at = exons2$exon_start, names = NULL, xlim = c(xmin, xmax), 
                pch = 20, cex = 0.5, axes = FALSE, ylim = c(-0.2, ymax), 
                ylab = "Expression log2(rpkm + 1)", xlab = "", las = 1, boxwex = (xmax-xmin)*0.03, varwidth = FALSE)#, border = "#dd4814")
        title(paste0("Expression of ", geneName, " in Tumor (n=", (dim(ddT2)[2]-1), ")\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
        axis(1, labels = exons2$exon, at = exons2$exon_start, col.axis = FALSE)
        axis(2, labels = seq(0, ymax, 1), at = seq(0, ymax, 1), las = 1)
        axis(3, labels = positions, at = positions, col.axis = FALSE) 
        
        par(xpd = FALSE)
        mtext(side = 1, text = exons2$exon, at = exons2$exon_start, las = 3, line = 1) 
        mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
        
        
        if(geneLine){
          if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
          if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        }
        box(lwd = 1.5)
        
      }
    }
    
    plotting_results <- list(ddN2 = ddN2, ddT2 = ddT2, exons2 = exons2)
    
    options(scipen=0)
    
    return(plotting_results)
    
  } else{
    
    stop(paste0("There are not exons in this region"))
    
  }
  
}
