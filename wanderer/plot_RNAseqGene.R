plot_RNAseqGene <- function(dd, geneName, tissue_label){
  

  ddN <- as.numeric(dd$Normal[-1])
  ddT <- as.numeric(dd$Tumor[-1])
  ymin <- min(c(ddN, ddT), na.rm = TRUE)
  ymax <- max(c(ddN, ddT), na.rm = TRUE)
  
  if(length(ddN)!=0 & length(ddT)!=0){
    pval <- sprintf('%e',wilcox.test(ddN, ddT)$p.val)
  } else{
    pval <- NULL
  }
  
  par(mfrow = c(1,2))
  
  boxplot(ddN, at = 0, pch = 20, border = "dodgerblue", ylab = "log2(normalized rsem + 1)", xlim = c(-0.5, 1.5), ylim = c(ymin, ymax), las = 1)
  boxplot(ddT, at = 1, pch = 20, border = "darkred", add = TRUE, axes=FALSE)
  axis(1, at = c(0, 1), tick = FALSE, c(paste0("Normal\nn=",length(ddN)), paste0("Tumor\nn=",length(ddT))))
  title(paste0("Expression of ", geneName, " in \n", tissue_label), sub = ifelse(is.null(pval),"",paste0("wilcoxon p.val = ", pval)))
  box(lwd=1.5)
  
  stripchart(ddN, at = 0, pch = 1, ylab = "log2(normalized rsem + 1)", cex = 0.5, vertical = TRUE, method = "jitter", jitter = 0.1, col = "dodgerblue", xlim = c(-0.5, 1.5), ylim = c(ymin, ymax), las = 1)
  stripchart(ddT, at = 1, pch = 1, vertical = TRUE, method = "jitter", jitter = 0.1, col = "darkred", add = TRUE, axes=FALSE)
  axis(1, at = c(0, 1), tick = FALSE, c(paste0("Normal\nn=",length(ddN)), paste0("Tumor\nn=",length(ddT))))
  title(paste0("Expression of ", geneName, " in \n", tissue_label), sub = ifelse(is.null(pval),"",paste0("wilcoxon p.val = ", pval)))
  box(lwd=1.5)
  
}
