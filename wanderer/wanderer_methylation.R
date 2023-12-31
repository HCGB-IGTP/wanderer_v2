#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

######################################
#function to plot a region with the probes and the gene

wanderer_methylation <- function(results_filt, geneName, geneNamesType, npointsN, npointsT, CpGislands, plotmean, plotting, geneLine, proportional){
  
  options(scipen=20)
  
  #data_label <- "450k Methylation Array"
  
  ddN2 <- results_filt$ddN2
  ddT2 <- results_filt$ddT2
  probes2 <- results_filt$probes2
  tissue_label <- results_filt$tissue_label
  
  if(proportional){
    xmin <- results_filt$xmin
    xmax <- results_filt$xmax
    
    #axis limits
    gmin <- unique(probes2$genestart[probes2[,paste0(geneNamesType)] == geneName])
    gmax <- unique(probes2$geneend[probes2[,paste0(geneNamesType)] == geneName])
    gstrand <- unique(probes2$genestrand[probes2[,paste0(geneNamesType)] == geneName])
    gchr <- unique(probes2$chr)
    if(length(gmin)==0 & length(gmax)>0) gmin <- xmin
    if(length(gmin)>0 & length(gmax)==0) gmax <- xmax
    if(length(gmin)==0 & length(gmax)==0) gmin <- gmax <- NULL
    
    
    posmin <- ((xmin%/%1000)-1)*1000
    posmax <- ((xmax%/%1000)+1)*1000
    postep <- posmax - posmin
    positions <- round(seq(posmin, posmax, postep/5),0)
    positions <- (positions%/%1000)*1000
    positions <- positions[positions>=xmin & positions<=xmax]
    
    if(positions[2]-positions[1]>=100000) escala <- 100000
    if(positions[2]-positions[1]<100000 & positions[2]-positions[1]>=10000) escala <- 10000
    if(positions[2]-positions[1]<10000 & positions[2]-positions[1]>=5000) escala <- 5000
    if(positions[2]-positions[1]<5000 & positions[2]-positions[1]>=1000) escala <- 1000
    if(positions[2]-positions[1]<1000) escala <- 0
  }
  
  if(!proportional){
    xmin <- 1
    xmax <- dim(probes2)[1]
    probes2$cg_start <- seq(xmin,xmax,1)
    
    gmin <- gmax <- gchr <- NULL
    escala <- 0
  }
  
  if(!is.null(ddN2)){
    if(dim(ddN2)[2]>2) sampN <- c(2:dim(ddN2)[2])
    if(dim(ddN2)[2]==2) sampN <- 2
  } else{
    sampN <- 0
  }
  if(!is.null(ddT2)){
    if(dim(ddT2)[2]>2) sampT <- c(2:dim(ddT2)[2])
    if(dim(ddT2)[2]==2) sampT <- 2
  } else{
    sampT <- 0
  }
  
  
  #profile plot
  if(plotting){
    #color to mark CpGislands
    colort <- rep("black", dim(probes2)[1])
    if(CpGislands) colort[!is.na(probes2$cpgiid)] <- "forestgreen"
    
    
    if(!plotmean){
      
      #sample of patients
      if(!is.null(ddN2)){
        if(npointsN==0) npointsN <- 30
        if(npointsN<=(dim(ddN2)[2]-1)){
          set.seed(1234)
          sampN <- sample(x=c(2:dim(ddN2)[2]),size=npointsN)
        }
      }
      
      if(!is.null(ddT2)){
        if(npointsT==0) npointsT <- 30
        if(npointsT<=(dim(ddT2)[2]-1)){
          set.seed(1234)
          sampT <- sample(x=c(2:dim(ddT2)[2]),size=npointsT)
        }
      }
      
      ###############
      #plot
      par(mfrow = c(2,1))
      
      #plot for Normals
      par(mai = par()$mai + c(0.3,0,1,0))
      if(!is.null(ddN2)){
        if(dim(ddN2)[2]>2){
          plot(probes2$cg_start, ddN2[,sampN[1]], type = "b", xlim = c(xmin, xmax), 
               lwd = 1.2, cex = 0.7, axes = FALSE, col = "dodgerblue", ylim = c(-0.2, 1), 
               ylab = "Methylation (beta value)", xlab = "", las = 1)
          for(pl in sampN[-1]) lines(probes2$cg_start, ddN2[,pl], type = "b", lwd = 1.2, cex = 0.7, col = "dodgerblue")
          if(proportional) title(paste0("Methylation of ", geneName, " in Normal (n=", length(sampN), ")\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
          if(!proportional) title(paste0("Methylation of ", geneName, " in Normal (n=", length(sampN), ")\n", tissue_label))
          
        }
        if(dim(ddN2)[2]==2){
          plot(probes2$cg_start, ddN2[,2], type = "b", xlim = c(xmin, xmax), 
               lwd = 1.2, cex = 0.7, axes = FALSE, col = "dodgerblue", ylim = c(-0.2, 1), 
               ylab = "Methylation (beta value)", xlab = "", las = 1)
          if(proportional) title(paste0("Methylation of ", geneName, " in Normal (n=1)\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
          if(!proportional) title(paste0("Methylation of ", geneName, " in Normal (n=1)\n", tissue_label))
          
        }
      } else{
        plot((xmax+xmin)/2, 0.5, type = "n", xlim = c(xmin, xmax), 
             axes = FALSE, ylim = c(-0.2, 1), 
             ylab = "Methylation (beta value)", xlab = "", las = 1)
        text((xmax+xmin)/2, 0.5, "No data available", col = "dodgerblue")
        if(proportional) title(paste0("Methylation of ", geneName, " in Normal (n=0)\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
        if(!proportional) title(paste0("Methylation of ", geneName, " in Normal (n=0)\n", tissue_label))
        
      }
      axis(1, labels = probes2$probe, at = probes2$cg_start, col.axis = FALSE)
      axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
      if(proportional) axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      par(xpd = FALSE)
      mtext(side = 1, text = probes2$probe, at = probes2$cg_start, las = 3, col = colort, line = 1) 
      if(proportional) mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
      if(proportional & geneLine & !is.null(gmin)){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      if(escala!=0 & proportional){
        par(xpd = TRUE)
        arrows(xmin, 1.25, xmin + escala, 1.25, lwd = 2,  length=0.05, angle=90, code=3)
        text(x = xmin + (escala/2), y = 1.3, labels = paste0((escala/1000),"Kb"))
      }
      par(xpd = FALSE)
      
      
      #plot for Tumors
      if(!is.null(ddT2)){
        if(dim(ddT2)[2]>2){
          plot(probes2$cg_start, ddT2[,sampT[1]], type = "b", xlim = c(xmin, xmax), 
               lwd = 1.2, cex = 0.7, axes = FALSE, col = "darkred", ylim = c(-0.2, 1),
               ylab = "Methylation (beta value)", xlab = "", las = 1)
          for(pl in sampT[-1]) lines(probes2$cg_start, ddT2[,pl], type = "b", lwd = 1.2, cex = 0.7, col = "darkred")
          if(proportional) title(paste0("Methylation of ", geneName, " in Tumor (n=", length(sampT), ")\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
          if(!proportional) title(paste0("Methylation of ", geneName, " in Tumor (n=", length(sampT), ")\n", tissue_label)) 
        }
        if(dim(ddT2)[2]==2){
          plot(probes2$cg_start, ddT2[,2], type = "b", xlim = c(xmin, xmax), 
               lwd = 1.2, cex = 0.7, axes = FALSE, col = "darkred", ylim = c(-0.2, 1),
               ylab = "Methylation (beta value)", xlab = "", las = 1)
          if(proportional) title(paste0("Methylation of ", geneName, " in Tumor (n=1)\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
          if(!proportional) title(paste0("Methylation of ", geneName, " in Tumor (n=1)\n", tissue_label))   
        }
      } else{
        plot((xmax+xmin)/2, 0.5, type = "n", xlim = c(xmin, xmax), 
             axes = FALSE, ylim = c(-0.2, 1), 
             ylab = "Methylation (beta value)", xlab = "", las = 1)
        text((xmax+xmin)/2, 0.5, "No data available", col = "darkred")
        if(proportional) title(paste0("Methylation of ", geneName, " in Tumor (n=0)\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
        if(!proportional) title(paste0("Methylation of ", geneName, " in Tumor (n=0)\n", tissue_label))
        
      }
      axis(1, labels = probes2$probe, at = probes2$cg_start, col.axis = FALSE)
      axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
      if(proportional) axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      par(xpd = FALSE)
      mtext(side = 1, text = probes2$probe, at = probes2$cg_start, las = 3, col = colort, line = 1)
      if(proportional) mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
      
      if(proportional & geneLine & !is.null(gmin)){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
    }
    
    if(plotmean){
      ###########
      #plot
      par(mfrow = c(2,1))
      
      #plot for Normals
      par(mai = par()$mai + c(0.3,0,1,0))
      if(!is.null(ddN2)){
        boxplot(t(ddN2[,2:dim(ddN2)[2]]), at = probes2$cg_start, names = NULL, xlim = c(xmin, xmax), 
                pch = 20, cex = 0.5, axes = FALSE, ylim = c(-0.2, 1), border = "dodgerblue",
                ylab = "Methylation (beta value)", xlab = "", las = 1, boxwex = (xmax-xmin)*0.03, varwidth = FALSE)
        if(proportional) title(paste0("Methylation of ", geneName, " in Normal (n=", (dim(ddN2)[2]-1), ")\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
        if(!proportional) title(paste0("Methylation of ", geneName, " in Normal (n=", (dim(ddN2)[2]-1), ")\n", tissue_label))
        
      } else{
        plot((xmax+xmin)/2, 0.5, type = "n", xlim = c(xmin, xmax), 
             axes = FALSE, ylim = c(-0.2, 1), 
             ylab = "Methylation (beta value)", xlab = "", las = 1)
        text((xmax+xmin)/2, 0.5, "No data available", col = "dodgerblue")
        if(proportional) title(paste0("Methylation of ", geneName, " in Normal (n=0)\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
        if(!proportional) title(paste0("Methylation of ", geneName, " in Normal (n=0)\n", tissue_label))
        
      }
      axis(1, labels = probes2$probe, at = probes2$cg_start, col.axis = FALSE)
      axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
      if(proportional) axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      par(xpd = FALSE)
      mtext(side = 1, text = probes2$probe, at = probes2$cg_start, las = 3, col = colort, line = 1)
      if(proportional) mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
      
      if(proportional & geneLine & !is.null(gmin)){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      if(escala!=0 & proportional){
        par(xpd = TRUE)
        arrows(xmin, 1.25, xmin + escala, 1.25, lwd = 2,  length=0.05, angle=90, code=3)
        text(x = xmin + (escala/2), y = 1.3, labels = paste0((escala/1000),"Kb"))
      }
      
      par(xpd = FALSE)
      
      #plot for Tumors
      if(!is.null(ddT2)){
        boxplot(t(ddT2[,2:dim(ddT2)[2]]), at = probes2$cg_start, names = NULL, xlim = c(xmin, xmax), 
                pch = 20, cex = 0.5, axes = FALSE, ylim = c(-0.2, 1), border="darkred",
                ylab = "Methylation (beta value)", xlab = "", las = 1, boxwex = (xmax-xmin)*0.03, varwidth = FALSE)
        if(proportional) title(paste0("Methylation of ", geneName, " in Tumor (n=", (dim(ddT2)[2]-1), ")\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
        if(!proportional) title(paste0("Methylation of ", geneName, " in Tumor (n=", (dim(ddT2)[2]-1), ")\n", tissue_label))
        
      } else{
        plot((xmax+xmin)/2, 0.5, type = "n", xlim = c(xmin, xmax), 
             axes = FALSE, ylim = c(-0.2, 1), 
             ylab = "Methylation (beta value)", xlab = "", las = 1)
        text((xmax+xmin)/2, 0.5, "No data available", col = "darkred")
        if(proportional) title(paste0("Methylation of ", geneName, " in Tumor (n=0)\n", tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
        if(!proportional) title(paste0("Methylation of ", geneName, " in Tumor (n=0)\n", tissue_label))
        
      }
      axis(1, labels = probes2$probe, at = probes2$cg_start, col.axis = FALSE)
      axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
      if(proportional) axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      par(xpd = FALSE)
      mtext(side = 1, text = probes2$probe, at = probes2$cg_start, las = 3, col = colort, line = 1)
      if(proportional) mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
      
      if(proportional & geneLine & !is.null(gmin)){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
    }
  }
  
  plotting_results <- list(ddN2 = ddN2, ddT2 = ddT2, probes2 = probes2)
  
  options(scipen=0)
  
  return(plotting_results)
  
}



