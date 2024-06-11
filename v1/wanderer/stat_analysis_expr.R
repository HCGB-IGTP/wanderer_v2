#############################################3
#statistical analysis for each expr
###############################################


stat_analysis_expr <- function(results_filt, geneName, geneNamesType, pvalThres, geneLine, plotting, proportional){
  
  ddN2 <- results_filt$ddN2
  row.names(ddN2) <- ddN2[,1]
  ddN <- ddN2[,-1]
  
  ddT2 <- results_filt$ddT2
  row.names(ddT2) <- ddT2[,1]
  ddT <- ddT2[,-1]
  
  exons2 <- results_filt$exons2
  
  if(is.null(dim(ddN)) | is.null(dim(ddT))){
    results_stats <- exons2
  } else{
    
    tissue_label <- results_filt$tissue_label
    
    
    #wilcoxon test
    results<-data.frame(exon=0,wilcox_stat=0,pval=0)
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
    
    results_stats<-merge(exons2,results,by="exon")
    results_stats<-results_stats[order(results_stats$exon_start),]
    aux1<-which(colnames(results_stats)=="emsemblgeneid")
    colnames(results_stats)[aux1]<-"ENSEMBL_geneID"
    aux2<-which(colnames(results_stats)=="emsembltransid")  
    colnames(results_stats)[aux2]<-"ENSEMBL_transcriptID"
    ########################################################################################################
    #plot
    if(plotting){
      old.scipen <- getOption("scipen")
      options("scipen"=999)
      asterisc<-results_stats$adj.pval<pvalThres
      asterisc[is.na(asterisc)]<-FALSE
      pasterisc<-results_stats$exon
      if(sum(asterisc)>0) pasterisc[asterisc]<-paste0("* ",results_stats$exon[asterisc])
      options("scipen"=old.scipen)
      
          
      if(proportional){
        xmin <- results_filt$xmin
        xmax <- results_filt$xmax
        
        #axis limits
        gmin <- unique(exons2$genestart[exons2[,paste0(geneNamesType)] == geneName])
        gmax <- unique(exons2$geneend[exons2[,paste0(geneNamesType)] == geneName])
        gstrand <- unique(exons2$strand[exons2[,paste0(geneNamesType)] == geneName])
        gchr <- unique(exons2$chr)
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
        xmax <- dim(exons2)[1]
        exons2$exon_start <- seq(xmin,xmax,1)
        
        gmin <- gmax <- gchr <- NULL
      }
      
      mddN <- apply(ddN,1,mean,na.rm=TRUE)
      mddT <- apply(ddT,1,mean,na.rm=TRUE)
      
      #y axis  
      ymax <- max(mddN, mddT)+1
      
      par(mai = par()$mai + c(2,0,1,0))
      plot(exons2$exon_start, mddN, type="b",xlim = c(xmin, xmax), 
           pch =  1, cex = 0.7, axes = FALSE, col = "dodgerblue", ylim = c(-0.2, ymax), 
           ylab = "Mean Expression log2(rpkm + 1)", xlab = "", las = 1, lwd=1.2)

      if(proportional) title(paste0("Mean Expression of ", geneName, "\n" ,tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
      if(!proportional) title(paste0("Mean Expression of ", geneName, "\n" ,tissue_label))
      
      axis(1, labels = pasterisc, at = exons2$exon_start, col.axis = FALSE)
      axis(2, labels = seq(0, ymax, 1), at = seq(0, ymax, 1), las = 1)
      if(proportional) axis(3, labels = positions, at = positions, col.axis = FALSE) 
      
      par(xpd = FALSE)
      mtext(side = 1, text = pasterisc, at = exons2$exon_start, las = 3, line = 1) 
      if(proportional) mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
      
      
      if(proportional & geneLine & !is.null(gmin)){
        if(gstrand == "-") arrows(gmax, -0.2, gmin, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
        if(gstrand == "+") arrows(gmin, -0.2, gmax, -0.2, cex = 2, col = "black", lwd = 2, length = 0.1)
      }
      box(lwd = 1.5)
      par(xpd=TRUE)
      legend(xmin, ymax + (ymax/3), c(paste0("Normal (n=", dim(ddN)[2], ")"), paste0("Tumor (n=", dim(ddT)[2], ")"), paste0("adj. pval<", pvalThres)), pch=c("","","*"), lty=c(1,1,0), lwd=c(1.2, 1.2, 0), col=c("dodgerblue","darkred","black"), yjust=0)
      par(xpd=FALSE)
      
      par(new=TRUE)
      plot(exons2$exon_start, mddT, type="b", xlim = c(xmin, xmax), ylim = c(-0.2, ymax), col="darkred", pch=1, axes=FALSE, cex = 0.7, lwd=1.2, xlab="", ylab="")
      
    }
  }
  
  return(results_stats)
  
}
