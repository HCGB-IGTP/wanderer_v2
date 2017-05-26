
## @author Anna Diez
## @fixup Izaskun Mallona
correl_meth_express <- function(geneName, probeID, ddmeth, ddGene, tissue_label, regressLine,
                                correlMethod, plotting, datareturn){

    save(file = '/tmp/wanderer_debug.RData',
         list = c('geneName', 'probeID', 'ddmeth', 'ddGene', 'tissue_label',
             'regressLine', 'correlMethod', 'plotting', 'datareturn'))
      
    ## [tn][em] stands for tumor/normal and expression/meth
    ## kept adiez's variable names too
    te <- ddN <- ddGene$Tumor
    tm <- ddT <-  ddmeth$ddT2
    ne <- ddNE <- ddGene$Normal
    nm <- ddTE <- ddmeth$ddN2

    if (is.null(ddN) | is.null(ddNE)) missatgeN <- "No data available"
    if (is.null(ddT) | is.null(ddTE)) missatgeT <- "No data available"

    for (item in c('te', 'tm', 'ne', 'nm')) {
        curr <- get(item)
        if (!is.null(curr)) {
            colnames(curr) <- strtrim(colnames(curr), 16)
            rownames(curr) <- curr[,1]
        }
        assign(x = item, value = curr)
    }

    common_tumors <- intersect(colnames(te), colnames(tm))
    if (length(common_tumors) > 0) {
        tumor_data <- data.frame(expr = as.numeric(te[,common_tumors]),
                                 meth = as.numeric(tm[probeID, common_tumors]))
        rownames(tumor_data) <- common_tumors
        tumor_data <- na.omit(tumor_data)
        correlT <- cor(tumor_data$meth, tumor_data$expr, method = correlMethod)
    } else {
        missatgeT <- "No common patients available"
    }

    common_normals <- intersect(colnames(ne), colnames(nm))
    if (length(common_normals) > 0) {
        normal_data <- data.frame(expr = as.numeric(ne[, common_normals]),
                                  meth = as.numeric(nm[probeID, common_normals]))
        rownames(normal_data) <- common_normals
        normal_data <- na.omit(normal_data)
        correlN <- cor(normal_data$meth, normal_data$expr, method = correlMethod)
    } else {
        missatgeN <- "No common patients available"
    }

    ## recovering adiez's naming
    if (exists('normal_data')) {
        ddNE <- ddNE2 <- normal_data$expr
        ddN <- ddN2 <- normal_data$meth
    }
    if (exists('tumor_data')) {
        ddTE <- ddTE2 <- tumor_data$expr
        ddT <- ddT2 <-tumor_data$meth        
    }

    if (plotting){
        if(!is.null(ddNE) & !is.null(ddTE)){
            ymax <- max(as.numeric(c(ddNE,ddTE)),na.rm=TRUE)
            ymin <- min(as.numeric(c(ddNE,ddTE)),na.rm=TRUE)
        }
        else if(is.null(ddNE)){
            ymax <- max(as.numeric(ddTE),na.rm=TRUE)
            ymin <- min(as.numeric(ddTE),na.rm=TRUE)
        }
        else if(is.null(ddTE)){
            ymax <- max(as.numeric(ddNE),na.rm=TRUE)
            ymin <- min(as.numeric(ddNE),na.rm=TRUE)
        }
        else if(is.null(ddNE) & is.null(ddTE)){
            ymax <- 1
            ymin <- 0
        }
        
        if (correlMethod=="spearman")
            titolcorrel <- "    Spearman rho="
        else if (correlMethod=="kendall")
            titolcorrel <- "    Kendall tau="
        else if (correlMethod=="pearson")
            titolcorrel <- "    Pearson r="

        ## plot start

        par(mfrow=c(1,2))
        
        if(!is.null(ddN) & !is.null(ddNE)){
            plot(as.numeric(ddN), as.numeric(ddNE), col="dodgerblue",
                 ylab=paste0("log2(normalized rsem + 1) of ", geneName),
                 xlab=paste0("methylation beta value of ", probeID),
                 xlim=c(0,1), ylim=c(ymin,ymax), las=1)
            if(length(ddN2)>1){
                title(paste0(geneName,"   Normal\n",tissue_label,"\nn=", length(ddN2),
                             titolcorrel, round(correlN,3)))
                if(regressLine) abline(lm(as.numeric(ddNE)~as.numeric(ddN)), col="dodgerblue", lwd=2)
            }
            if(length(ddN2)==1) title(paste0(geneName,"   Normal\n",tissue_label,"\nn=",dim(ddN)[2]))
            
        } else{
            plot(0.5, (ymax+ymin)/2, type = "n",
                 ylab=paste0("log2(normalized rsem + 1) of ", geneName), xlab=paste0("methylation beta value of ", probeID), xlim=c(0,1), ylim=c(ymin,ymax), las=1)
            text(0.5, (ymax+ymin)/2, paste0(missatgeN), col = "dodgerblue")
            title(paste0(geneName,"   Normal\n",tissue_label,"\nn=0"))
            
        }
        box(lwd=2)
        
        if(!is.null(ddT) & !is.null(ddTE)){
            plot(as.numeric(ddT),as.numeric(ddTE),col="darkred",
                 ylab=paste0("log2(normalized rsem + 1) of ", geneName),
                 xlab=paste0("methylation beta value of ", probeID),xlim=c(0,1),ylim=c(ymin,ymax),las=1)
            if(length(ddT2)>1){
                title(paste0(geneName,"   Tumor\n",tissue_label,"\nn=", length(ddT2), titolcorrel, round(correlT,3)))
                if(regressLine) abline(lm(as.numeric(ddTE)~as.numeric(ddT)), col="darkred", lwd=2)
            }
            if(length(ddT2)==1) title(paste0(geneName,"   Tumor\n",tissue_label,"\nn=",dim(ddT)[2]))
        } else{
            plot(0.5, (ymax+ymin)/2, type = "n",
                 ylab=paste0("log2(normalized rsem + 1) of ", geneName),
                 xlab=paste0("methylation beta value of ", probeID), xlim=c(0,1),
                 ylim=c(ymin,ymax), las=1)
            text(0.5, (ymax+ymin)/2, paste0(missatgeT), col = "dodgerblue")
            title(paste0(geneName,"   Tumor\n",tissue_label,"\nn=0"))
            
        }
        box(lwd=2)
    }

    if(!is.null(ddNE)){
        ddNE2 <- data.frame(PATIENT = rownames(normal_data),
                            gene = normal_data$expr)
        colnames(ddNE2)[2] <- geneName
    } else {
        ddNE2 <- NULL
    }

    if(!is.null(ddTE)){
        ddTE2 <- data.frame(PATIENT = rownames(tumor_data),
                            gene = tumor_data$expr)
        colnames(ddTE2)[2] <- geneName
        
    } else{
        ddTE2 <- NULL
    }

    tryCatch({
        save(file = '/tmp/wanderer_debug.RData',
             list = c('geneName', 'probeID', 'ddmeth', 'ddGene', 'tissue_label',
                 'regressLine', 'correlMethod', 'plotting', 'datareturn', 'ddNE2', 'ddTE2'))
    }, error = function(x) print('foo'))
    
    if(datareturn)  return(list(ddNE=ddNE2, ddTE=ddTE2))
}


