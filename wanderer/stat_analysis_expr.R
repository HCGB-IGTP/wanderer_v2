#############################################3
#statistical analysis for each expr
###############################################


stat_analysis_expr <- function(results_filt, geneName, geneNamesType, geneLine, plotting){
  
  ddN2 <- results_filt$ddN2
  row.names(ddN2) <- ddN2[,1]
  ddN <- ddN2[,-1]
  
  ddT2 <- results_filt$ddT2
  row.names(ddT2) <- ddT2[,1]
  ddT <- ddT2[,-1]
  
  exons2 <- results_filt$exons2
  tissue_label <- results_filt$tissue_label
  xmin <- results_filt$xmin
  xmax <- results_filt$xmax
  
  
  grupsf <- c(rep("Normal",dim(ddN)[2]), rep("Tumor",dim(ddT)[2]))
  nom_sample <- c(colnames(ddN), colnames(ddT))
  
  dd <- cbind(ddN, ddT)
  colnames(dd) <- paste0(grupsf, "_", nom_sample)
  
  
  noms<-unique(grupsf) 
  
  #ajustem el model
  fit<-lmFit(dd,design=model.matrix(~grupsf),na.rm=TRUE)
  colnames(fit$coefficients)<-row.names(contrasts(as.factor(grupsf)))
  #moderated-t
  fit2<-contrasts.fit(fit,contrasts=contrasts(as.factor(grupsf)))
  qwe<-eBayes(fit2)
  
  
  #taula final
  results<-topTable(qwe,number=dim(dd)[1],adjust="BH")
  results<-data.frame(id=row.names(results),results,stringsAsFactors=FALSE)
  results_stats<-merge(exons2,results,by.y="id",by.x="exon")
  
  ########################################################################################################
  #plot
  if(plotting){
    
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
    
    
    mddN <- apply(ddN,1,mean,na.rm=TRUE)
    mddT <- apply(ddT,1,mean,na.rm=TRUE)
    
    #y axis  
    ymax <- max(ddT[,-1], ddN[,-1])
    
    par(mai = par()$mai + c(0.3,0,1,0))
    layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(3.2,0.8))
    
    plot(exons2$exon_start, mddN, type="b",xlim = c(xmin, xmax), 
         pch =  1, cex = 0.7, axes = FALSE, col = "forestgreen", ylim = c(-0.2, ymax), 
         ylab = "Mean Expression log2(rpkm + 1)", xlab = "", las = 1, lwd=1.2) 
    lines(exons2$exon_start, mddT, col="darkred", lwd=1.2)
    points(exons2$exon_start, mddT, col="darkred", pch=1)
    
    
    title(paste0("Mean Expression of ", geneName, "\n" ,tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
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

    plot(1, type = "n", xlim = c(0.5,1.5), ylim = c(0.5, 1.5), axes = FALSE, col="white", xlab = "", ylab = "")
    par(xpd=TRUE)
    legend(0.2, 1.5, c(paste0("Normal (n=", dim(ddN)[2], ")"), paste0("Tumor (n=", dim(ddT)[2], ")")), lty=1, lwd=1.2, col=c("forestgreen","darkred"))
    par(xpd=FALSE)
    
  }
    return(results_stats)
}
