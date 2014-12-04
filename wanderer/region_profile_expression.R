#!/usr/bin/env R
#
# @package wanderer
# @author Anna Diez

######################################
#function to plot a region with the exons and the gene

region_profile_expression <- function(con, geneName, geneNamesType, sampleSize, tissue, zoom, walk, npointsN, npointsT, plotmean, plotting, geneLine){
  
  tissue_label <- sampleSize[sampleSize[,3]==tissue,1]
  data_label <- "Illumina HiSeq RNAseq"
  
  if (geneNamesType == "genename") geneNamesType_label <- "Gene Name"
  if (geneNamesType == "emsemblgeneid") geneNamesType_label <- "Ensembl Gene ID"
  
  ################
  #RNAseq annotation download
  probes <- dbSendQuery(con, statement = paste0("select * from annotations.exons_annot
  where chr = (select chr from annotations.exons_annot where ", geneNamesType, " = '", geneName, "' limit 1)
  and exon_start > (select genestart as sstart from annotations.exons_annot
     where ", geneNamesType, " = '", geneName, "' limit 1)
  and exon_end < (select geneend as send from annotations.exons_annot
      where ", geneNamesType, " = '", geneName, "' limit 1);"))
  
  probes <- fetch(probes, n = -1)
  if(dim(probes)[1] == 0) stop(paste0("There are not exons in this region"))
  if(dim(probes)[1] != 0){
    ordering <- order(probes$exon_start)
    probes <- probes[ordering,]
    
    ##############
    probes_collapsed <- paste0("('", paste(probes$exon, collapse="','"), "')")
    
    #methylation array tumor data download
    ddT <- dbSendQuery(con, paste0("select * from ", tissue, "_tumor.illuminahiseq_rnaseqv2
                                   where exon in ", probes_collapsed, ";"))
    ddT <- fetch(ddT, n = -1)
    ddT <- ddT[ordering,]
    ddT[,2:dim(ddT)[2]] <- log2(ddT[,2:dim(ddT)[2]]+1)
    
    #methylation array normal data download
    ddN <- dbSendQuery(con, paste0("select * from ", tissue, "_normal.illuminahiseq_rnaseqv2
                                   where exon in ", probes_collapsed, ";"))
    ddN <- fetch(ddN, n = -1)
    ddN <- ddN[ordering,]
    ddN[,2:dim(ddN)[2]] <- log2(ddN[,2:dim(ddN)[2]]+1)
    
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
    if(dim(probes)[1] == 0) stop(paste0("There are not exons in this region"))
    
  }
  
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
           pch = 20, xaxt = "n", col = "black", ylim = c(-0.2, ymax), 
           ylab = "log2(rpkm + 1)", xlab = "", las = 1)
      for(pl in sampN[-1]) lines(probes$exon_start, ddN[,pl], type = "b", pch = 20, col = "black")
      title(paste0("Expression of ", geneName, " ", "Normal ", tissue, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$exon, at = probes$exon_start, las=3)
      if(geneLine){
        if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      #plot for Tumors
      plot(probes$exon_start, ddT[,sampT[1]], type = "b", xlim = c(xmin, xmax), 
           pch = 20, xaxt = "n", col = "black", ylim = c(-0.2, ymax),
           ylab = "log2(rpkm + 1)", xlab = "", las = 1)
      for(pl in sampT[-1]) lines(probes$exon_start, ddT[,pl], type = "b", pch = 20, col = "black")
      title(paste0("Expression of ", geneName, " ", "Tumor ", tissue, "\n", gchr, ": ", xmin, " - ", xmax))
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
              ylab = "log2(rpkm + 1)", xlab = "", las = 1, boxwex = 2000)
      title(paste0("Expression of ", geneName, " ", "Normal ", tissue, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$exon, at = probes$exon_start, las=3)
      if(geneLine){
        if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
      #plot for Tumors
      boxplot(t(ddT[,2:dim(ddT)[2]]), at = probes$exon_start, names = NULL, xlim = c(xmin, xmax), 
              pch = 20, cex = 0.5, xaxt = "n", ylim = c(-0.2, ymax), 
              ylab = "log2(rpkm + 1)", xlab = "", las = 1, boxwex = 2000)#, border = "#dd4814")
      title(paste0("Expression of ", geneName, " ", "Tumor ", tissue, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = probes$exon, at = probes$exon_start, las=3)
      if(geneLine){
        if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      
    }
  }
  
  results <- list(ddN2, ddT2, probes2)
  return(results)
}


