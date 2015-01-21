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
  
  
  grups <- c(rep("Normal",dim(ddN)[2]), rep("Tumor",dim(ddT)[2]))
  nom_sample <- c(colnames(ddN), colnames(ddT))
  
  dd <- cbind(ddN, ddT)
  colnames(dd) <- paste0(grups, "_", nom_sample)
  
  
  #per poder calcular l'M-value no podem tenir CM=0 o CM=1
  for (i in 1:dim(dd)[2]){
    dd[dd[,i]==0 & !is.na(dd[,i]),i]<-1/(100)
    dd[dd[,i]==1 & !is.na(dd[,i]),i]<-1-(1/(100))
  }
  
  #M-values
  mvals<-function(x){
    mvals<-log2(x/(1-x))
  }
  mdd<-mvals(dd)
  
  noms<-unique(grups) 
  #taula_resum<-noms_taula<-0
  
  dd1<-mdd[,grups==noms[1]]
  dd2<-mdd[,grups==noms[2]]
  grupsf<-c(rep(noms[1],dim(dd1)[2]),rep(noms[2],dim(dd2)[2]))
  
  #ajustem el model
  ## fit<-lmFit(mdd,design=model.matrix(~grupsf),na.rm=TRUE)
  fit<-lmFit(mdd,design=model.matrix(~as.factor(grupsf)),na.rm=TRUE)
  colnames(fit$coefficients)<-row.names(contrasts(as.factor(grupsf)))
  #moderated-t
  fit2<-contrasts.fit(fit,contrasts=contrasts(as.factor(grupsf)))
  qwe<-eBayes(fit2)
  
  
  #taula final
  results<-topTable(qwe,number=dim(dd)[1],adjust="BH")
  results<-data.frame(id=row.names(results),results,stringsAsFactors=FALSE)
  results_stats<-merge(probes2,results,by.y="id",by.x="probe")
  
  ########################################################################################################
  #plot
  if(plotting){
    
    #axis limits
    gmin <- unique(probes2$genestart[probes2[,paste0(geneNamesType)] == geneName])
    gmax <- unique(probes2$geneend[probes2[,paste0(geneNamesType)] == geneName])
    gstrand <- unique(probes2$genestrand[probes2[,paste0(geneNamesType)] == geneName])
    gchr <- unique(probes2$chr[probes2[,paste0(geneNamesType)] == geneName])
    
    
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
    
    par(mai = par()$mai + c(0.3,0,1,0))
    layout(matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(3.2,0.8))
    
    plot(probes2$cg_start, mddN, type="b",xlim = c(xmin, xmax), 
         pch =  1, cex = 0.7, axes = FALSE, col = "forestgreen", ylim = c(-0.2, 1), 
         ylab = "Mean Methylation (beta value)", xlab = "", las = 1, lwd=1.2) 
    lines(probes2$cg_start, mddT, col="darkred", lwd=1.2)
    points(probes2$cg_start, mddT, col="darkred", pch=1)
    title(paste0("Mean Methylation of ", geneName, "\n" ,tissue_label, "\n", gchr, ": ", xmin, " - ", xmax))
    axis(1, labels = probes2$probe, at = probes2$cg_start, col.axis = FALSE)
    axis(2, labels = seq(0, 1, 0.2), at = seq(0, 1, 0.2), las = 1)
    axis(3, labels = positions, at = positions, col.axis = FALSE) 
    mtext(side = 1, text = probes2$probe, at = probes2$cg_start, las = 3, col = colort, line = 1) 
    mtext(side = 3, text = format(positions, big.mark=','), at = positions, las = 1, line = 1, cex = 1) 
    if(geneLine){
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
  
