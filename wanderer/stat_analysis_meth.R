#############################################3
#statistical analysis for each probe
###############################################


stat_analysis_meth <- function(results_filt, geneName, geneNamesType, CpGislands, geneLine, plotting){
  
  ddN2 <- results_filt$ddN2
  row.names(ddN2) <- ddN2[,1]
  ddN <- ddN2[,-1]
  
  ddT2 <- results_filt$ddT2
  row.names(ddT2) <- ddT2[,1]
  ddT <- ddT2[,-1]
  
  probes2 <- results_filt$probes2
  tissue_label <- results_filt$tissue_label
  xmin <- results_filt$xmin
  xmax <- results_filt$xmax
  
  
  #wilcoxon test
  results<-data.frame(probe=0,wilcox_stat=0,pval=0)
  for(j in 1:dim(ddN)[1]){
    test<-wilcox.test(x=as.numeric(ddN[j,]),y=as.numeric(ddT[j,]))
    results[j,]<-c(rownames(ddN)[j],test$statistic,test$p.value)
  }
    
    results$adj.pval<-p.adjust(results$pval, method="BH")
    results_stats<-merge(probes2,results,by="probe")
    results_stats<-results_stats[order(results_stats$cg_start),]
    aux1<-which(colnames(results_stats)=="emsemblgeneid")
    colnames(results_stats)[aux1]<-"ENSEMBL_geneID"
  
    ########################################################################################################
    #plot
    if(plotting){
      
      asterisc<-results_stats$adj.pval<0.05
      pasterisc<-probes2$probe
      pasterisc[asterisc]<-paste0("* ",probes2$probe[asterisc])
      
      
      #axis limits
      gmin <- unique(probes2$genestart[probes2[,paste0(geneNamesType)] == geneName])
      gmax <- unique(probes2$geneend[probes2[,paste0(geneNamesType)] == geneName])
      gstrand <- unique(probes2$genestrand[probes2[,paste0(geneNamesType)] == geneName])
      gchr <- unique(probes2$chr[probes2[,paste0(geneNamesType)] == geneName])
      if(length(gmin)==0 & length(gmax)>0) gmin <- xmin
      if(length(gmin)>0 & length(gmax)==0) gmax <- xmax
      if(length(gmin)==0 & length(gmax)==0) gmin <- gmax <- NULL
      
      posmin <- ((xmin%/%1000)-1)*1000
      posmax <- ((xmax%/%1000)+1)*1000
      postep <- posmax - posmin
      positions <- round(seq(posmin, posmax, postep/5),0)
      positions <- (positions%/%1000)*1000
      
      
      mddN <- apply(ddN,1,mean,na.rm=TRUE)
      mddT <- apply(ddT,1,mean,na.rm=TRUE)
      
      
      #color to mark CpGislands
      colort <- rep("black", dim(probes2)[1])
      if(CpGislands) colort[!is.na(probes2$cpgiid)] <- "forestgreen"
      
      par(mai = par()$mai + c(0.7,0,1,0))
      layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(3.2,0.8))
      
      plot(probes2$cg_start, mddN, type="b",xlim = c(xmin, xmax), 
           pch =  1, cex = 0.7, axes = FALSE, col = "forestgreen", ylim = c(-0.2, 1), 
           ylab = "Mean Methylation (beta value)", xlab = "", las = 1, lwd=1.2) 
      lines(probes2$cg_start, mddT, col="darkred", lwd=1.2)
      points(probes2$cg_start, mddT, col="darkred", pch=1, cex = 0.7)
      title(paste0("Mean Methylation of ", geneName, "\n" ,tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      axis(1, labels = pasterisc, at = probes2$cg_start, col.axis = FALSE)
      axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
      axis(3, labels = positions, at = positions, col.axis = FALSE) 
      mtext(side = 1, text = pasterisc, at = probes2$cg_start, las = 3, col = colort, line = 1) 
      mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
      if(geneLine & !is.null(gmin)){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "#dd4814", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      plot(1, type = "n", xlim = c(0.5,1.5), ylim = c(0.5, 1.5), axes = FALSE, col="white", xlab = "", ylab = "")
      par(xpd=TRUE)
      legend(0.2, 1.5, c(paste0("Normal (n=", dim(ddN)[2], ")"), paste0("Tumor (n=", dim(ddT)[2], ")")), lty=1, lwd=1.2, col=c("forestgreen","darkred"))
      par(xpd=FALSE)
      
    }
    return(results_stats)
  }
  