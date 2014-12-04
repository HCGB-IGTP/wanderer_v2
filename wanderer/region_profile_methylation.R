#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

######################################
#function to plot a region with the probes and the gene

region_profile_methylation <- function(con, geneName, geneNamesType, sampleSize, tissue, zoom, walk, npointsN, npointsT, CpGislands, plotmean, plotting, geneLine){

  tissue_label <- sampleSize[sampleSize[,3]==tissue,1]
  data_label <- "450k Methylation Array"
  
  if (geneNamesType == "genename") geneNamesType_label <- "Gene Name"
  if (geneNamesType == "emsemblgeneid") geneNamesType_label <- "Ensembl Gene ID"
  
  

  ################
  #methylation array annotation download
  probes <- dbSendQuery(con, statement = paste0("select * from annotations.humanmethylation450kannot
  where chr = (select chr from annotations.humanmethylation450kannot where ", geneNamesType, " = '", geneName, "' limit 1)
  and probestart > (select genestart-", 10000, " as sstart from annotations.humanmethylation450kannot
     where ", geneNamesType, " = '", geneName, "' limit 1)
  and probeend < (select geneend+", 10000, " as send from annotations.humanmethylation450kannot
      where ", geneNamesType, " = '", geneName, "' limit 1);"))
  
  probes <- fetch(probes, n = -1)
  print(dim(probes))
  if(dim(probes)[1] == 0) stop(paste0("There are not probes in this region"))
  if(dim(probes)[1] != 0){
    probes <- probes[order(probes$probe),]
    ordering <- order(probes$probestart)
    probes <- probes[ordering,]
    ##############
    probes_collapsed <- paste0("('", paste(probes$probe, collapse="','"), "')")
    
    #methylation array tumor data download
    ddT <- dbSendQuery(con, paste0("select * from ", tissue, "_tumor.humanmethylation450
  where probe in ", probes_collapsed, ";"))
    ddT <- fetch(ddT, n = -1)
    print(dim(ddT))
    ddT <- ddT[order(ddT$probe),]
    ddT <- ddT[ordering,]
    
    #methylation array normal data download
    ddN <- dbSendQuery(con, paste0("select * from ", tissue, "_normal.humanmethylation450
  where probe in ", probes_collapsed, ";"))
    ddN <- fetch(ddN, n = -1)
    print(dim(ddN))
    ddN <- ddN[order(ddN$probe),]
    ddN <- ddN[ordering,]
    
    #removing NA's
    naT <- which(apply(is.na(ddT), 1, sum) == (dim(ddT)[2] - 1))
    if(length(naT)>0){
      probes <- probes[-naT,]
      ddN <- ddN[-naT,]
      ddT <- ddT[-naT,]
    }
    
    naN <- which(apply(is.na(ddN), 1, sum) == (dim(ddN)[2] - 1))
    if(length(naN)>0){
      probes <- probes[-naN,]
      ddT <- ddT[-naN,]
      ddN <- ddN[-naN,]
    }
  }
  if(dim(probes)[1] == 0) stop(paste0("There are not probes in this region"))
  
  
  #axis limits
  gmin <- unique(probes$genestart[probes[,paste0(geneNamesType)] == geneName])
  gmax <- unique(probes$geneend[probes[,paste0(geneNamesType)] == geneName])
  gstrand <- unique(probes$genestrand[probes[,paste0(geneNamesType)] == geneName])
  gchr <- unique(probes$chr[probes[,paste0(geneNamesType)] == geneName])
  
  if(is.null(walk)) walk = 0 
  if(is.null(zoom)) zoom = 0 
  xmin <- gmin + zoom + walk
  xmax <- gmax - zoom + walk
  
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
      par(mai = par()$mai + c(0.3,0,0,0))
      plot(probes$cg_start, ddN[,sampN[1]], type = "b", xlim = c(xmin, xmax), 
           pch = 20, xaxt = "n", col = "black", ylim = c(-0.2, 1), 
           ylab = "beta value", xlab = "", las = 1)
      for(pl in sampN[-1]) lines(probes$cg_start, ddN[,pl], type = "b", pch = 20, col = "black")
      title(paste0("Methylation of ", geneName, " ", "Normal ", tissue, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$probe, at = probes$cg_start, col.axis = FALSE) 
      mtext(side = 1, text = probes$probe, at = probes$cg_start, las = 3, col = colort, line = 1) 
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      #plot for Tumors
      plot(probes$cg_start, ddT[,sampT[1]], type = "b", xlim = c(xmin, xmax), 
           pch = 20, xaxt = "n", col = "black", ylim = c(-0.2, 1),
           ylab = "beta value", xlab = "", las = 1)
      for(pl in sampT[-1]) lines(probes$cg_start, ddT[,pl], type = "b", pch = 20, col = "black")
      title(paste0("Methylation of ", geneName, " ", "Tumor ", tissue, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$probe, at = probes$cg_start, col.axis = FALSE) 
      mtext(side = 1, text = probes$probe, at = probes$cg_start, las = 3, col = colort, line = 1)
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
      par(mai = par()$mai + c(2,0,0,0))
      boxplot(t(ddN[,2:dim(ddN)[2]]), at = probes$cg_start, names = NULL, xlim = c(xmin, xmax), 
              pch = 20, cex = 0.5, xaxt = "n", ylim = c(-0.2, 1), 
              ylab = "beta value", xlab = "", las = 1, boxwex = 2000)
      title(paste0("Methylation of ", geneName, " ", "Normal ", tissue, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$probe, at = probes$cg_start, col.axis = FALSE) 
      mtext(side = 1, text = probes$probe, at = probes$cg_start, las = 3, col = colort, line = 1)
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      #plot for Tumors
      boxplot(t(ddT[,2:dim(ddT)[2]]), at = probes$cg_start, names = NULL, xlim = c(xmin, xmax), 
              pch = 20, cex = 0.5, xaxt = "n", ylim = c(-0.2, 1), 
              ylab = "beta value", xlab = "", las = 1, boxwex = 2000)
      title(paste0("Methylation of ", geneName, " ", "Tumor ", tissue, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$probe, at = probes$cg_start, col.axis = FALSE) 
      mtext(side = 1, text = probes$probe, at = probes$cg_start, las = 3, col = colort, line = 1)
      if(geneLine){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
    }
  }
  results <- list(ddN2, ddT2, probes2)
  return(results)
  
}



