#############################################3
#statistical analysis for each probe
###############################################


stat_analysis_meth <- function(results_filt, geneName, geneNamesType, CpGislands, geneLine, plotting, proportional){
  
  ddN2 <- results_filt$ddN2
  row.names(ddN2) <- ddN2[,1]
  ddN <- ddN2[,-1]
  
  ddT2 <- results_filt$ddT2
  row.names(ddT2) <- ddT2[,1]
  ddT <- ddT2[,-1]
  probes2 <- results_filt$probes2
  
  if(is.null(dim(ddN)) | is.null(dim(ddT))){
    results_stats <- probes2
  } else{
    
    tissue_label <- results_filt$tissue_label
    
    #wilcoxon test
    results<-data.frame(probe=0,wilcox_stat=0,pval=0)
    for(j in 1:dim(ddN)[1]){
      test<-wilcox.test(x=as.numeric(ddN[j,]),y=as.numeric(ddT[j,]))
      results[j,]<-c(rownames(ddN)[j],test$statistic,test$p.value)
    }
    
    results$adj.pval<-p.adjust(results$pval, method="BH")
    results$Norm_nsamples<-rep(dim(ddN)[2],dim(results)[1])
    results$Norm_mean<-apply(ddN,1,mean,na.rm=TRUE)
    results$Norm_sd<-apply(ddN,1,sd,na.rm=TRUE)
    results$Tum_nsamples<-rep(dim(ddT)[2],dim(results)[1])
    results$Tum_mean<-apply(ddT,1,mean,na.rm=TRUE)
    results$Tum_sd<-apply(ddT,1,sd,na.rm=TRUE)
    results<-results[,c(1,5:10,2:4)]
    results_stats<-merge(probes2,results,by="probe")
    results_stats<-results_stats[order(results_stats$cg_start),]
    aux1<-which(colnames(results_stats)=="emsemblgeneid")
    colnames(results_stats)[aux1]<-"ENSEMBL_geneID"
    
    ########################################################################################################
    #plot
    if(plotting){
      
      asterisc<-results_stats$adj.pval<0.05
      asterisc[is.na(asterisc)]<-FALSE
      pasterisc<-probes2$probe
      if(sum(asterisc)>0) pasterisc[asterisc]<-paste0("* ",probes2$probe[asterisc])
      
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
      }
      
      if(!proportional){
        xmin <- 1
        xmax <- dim(probes2)[1]
        probes2$cg_start <- seq(xmin,xmax,1)
        
        gmin <- gmax <- gchr <- NULL
      }
      
      
      
      mddN <- apply(ddN,1,mean,na.rm=TRUE)
      mddT <- apply(ddT,1,mean,na.rm=TRUE)
      
      
      #color to mark CpGislands
      colort <- rep("black", dim(probes2)[1])
      if(CpGislands) colort[!is.na(probes2$cpgiid)] <- "forestgreen"
      
      par(mai = par()$mai + c(0.7,0,1,0))
      plot(probes2$cg_start, mddN, type="b",xlim = c(xmin, xmax), 
           pch =  1, cex = 0.7, axes = FALSE, col = "dodgerblue", ylim = c(-0.2, 1), 
           ylab = "Mean Methylation (beta value)", xlab = "", las = 1, lwd=1.2) 
      lines(probes2$cg_start, mddT, col="darkred", lwd=1.2)
      points(probes2$cg_start, mddT, col="darkred", pch=5, cex = 0.7)
      if(proportional) title(paste0("Mean Methylation of ", geneName, "\n" ,tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      if(!proportional) title(paste0("Mean Methylation of ", geneName, "\n" ,tissue_label))
      
      axis(1, labels = pasterisc, at = probes2$cg_start, col.axis = FALSE)
      axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
      if(proportional) axis(3, labels = positions, at = positions, col.axis = FALSE) 
      mtext(side = 1, text = pasterisc, at = probes2$cg_start, las = 3, col = colort, line = 1) 
      if(proportional) mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
      if(proportional & geneLine & !is.null(gmin)){
        if(gstrand == -1) arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
        if(gstrand == 1) arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      par(xpd=TRUE)
      legend(xmin, 1.5, c(paste0("Normal (n=", dim(ddN)[2], ")"), paste0("Tumor (n=", dim(ddT)[2], ")")), lty=1, lwd=1.2, col=c("dodgerblue","darkred"))
      par(xpd=FALSE)
      
    }
  }
  
  return(results_stats)
  
}
