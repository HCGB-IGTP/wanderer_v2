#!/usr/bin/env R

load('/tmp/wanderer_debug.RData')
ls()

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
    colnames(curr) <- strtrim(colnames(curr), 16)
    rownames(curr) <- curr[,1]
    assign(x = item, value = curr)
}


common_tumors <- intersect(colnames(te), colnames(tm))
if (length(common_tumors) > 0) {
    tumor_data <- data.frame(expr = as.numeric(te[,common_tumors]),
                             meth = as.numeric(tm[probeID, common_tumors]))
    tumor_data <- na.omit(tumor_data)
    correlT <- cor(tumor_data$meth, tumor_data$expr, method = correlMethod)
} else {
    missatgeT <- "No common patients available"
}


common_normals <- intersect(colnames(ne), colnames(nm))
if (length(common_normals) > 0) {
    normal_data <- data.frame(expr = as.numeric(ne[, common_normals]),
                              meth = as.numeric(nm[probeID, common_normals]))
    normal_data <- na.omit(normal_data)
    correlN <- cor(normal_data$meth, normal_data$expr, method = correlMethod)
} else {
    missatgeN <- "No common patients available"
}

## recovering adiez's naming
ddTE <- ddTE2 <- tumor_data$expr
ddNE <- ddNE2 <- normal_data$expr
ddT <- ddT2 <-tumor_data$meth
ddN <- ddN2 <- normal_data$meth

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

    
    
## plot end

## let's try to replicate ddNE2
dim(ddNE2)
dim(ddTE2)



head(tm)

## anna start
  
ddN <- ddmeth$ddN2
ddT <- ddmeth$ddT2
ddNE <- ddGene$Normal
ddTE <- ddGene$Tumor

if(is.null(ddN) | is.null(ddNE)) missatgeN <- "No data available"
if(is.null(ddT) | is.null(ddTE)) missatgeT <- "No data available"

if(!is.null(ddN)){
    coln <- colnames(ddN)[-1]
    ddN <- ddN[ddN$PROBE==probeID,]
    row.names(ddN) <- ddN[,1]
    ddN <- as.data.frame(ddN[,-1])
    colnames(ddN) <- coln
}

if(!is.null(ddNE)){
    coln <- colnames(ddNE)[-1]
    row.names(ddNE) <- ddNE[,1]
    ddNE <- as.data.frame(ddNE[,-1])
    colnames(ddNE) <- coln
}

if(!is.null(ddT)){
    coln <- colnames(ddT)[-1]
    ddT <- ddT[ddT$PROBE==probeID,]
    row.names(ddT) <- ddT[,1]
    ddT <- as.data.frame(ddT[,-1])
    colnames(ddT) <- coln
}

if(!is.null(ddTE)){
    coln <- colnames(ddTE)[-1]
    row.names(ddTE) <- ddTE[,1]
    ddTE <- as.data.frame(ddTE[,-1])
    colnames(ddTE) <- coln
}

pTE <- sapply(strsplit(colnames(ddTE),"-"),"[",3)
pT <- sapply(strsplit(colnames(ddT),"-"),"[",3)

colT <- colnames(ddT)[pT%in%pTE]
colTE <- colnames(ddTE)[pTE%in%pT]
ddT <- as.data.frame(ddT[,pT%in%pTE])
colnames(ddT) <- colT
ddTE <- as.data.frame(ddTE[,pTE%in%pT])
colnames(ddTE) <- colTE
if(length(ddT)==0 | length(ddTE)==0){
    ddT <- ddTE <-NULL
    missatgeT <- "No common patients available"
} else{
    pTE <- sapply(strsplit(colnames(ddTE),"-"),"[",3)
    pT <- sapply(strsplit(colnames(ddT),"-"),"[",3)
    if(length(pT)!=length(pTE)){
        numpatients <- c(length(pT),length(pTE))
        aux <- which(numpatients==max(numpatients))
        if(aux==1){
            ddnew <- as.data.frame(matrix(0,ncol=numpatients[aux],nrow=dim(ddTE)[1]))
            for(i in 1:numpatients[aux]){
                aux2 <- which(pTE==pT[i])
                ddnew[,i] <- ddTE[,aux2]
                colnames(ddnew)[i] <- colnames(ddTE)[aux2]
            }
            ddTE <- ddnew
        }
        if(aux==2){
            ddnew <- as.data.frame(matrix(0,ncol=numpatients[aux],nrow=dim(ddT)[1]))
            for(i in 1:numpatients[aux]){
                aux2 <- which(pT==pTE[i])
                ddnew[,i] <- ddT[,aux2]
                colnames(ddnew)[i] <- colnames(ddT)[aux2]
            }
            ddT <- ddnew
        } 
    }
}
if(!is.null(ddT) & !is.null(ddTE)){
    ddT2<-as.numeric(ddT)
    ddTE2<-as.numeric(ddTE)
    ddTE2<-ddTE2[!is.na(ddT2)]
    ddT2<-ddT2[!is.na(ddT2)]
    ddT2<-ddT2[!is.na(ddTE2)]
    ddTE2<-ddTE2[!is.na(ddTE2)]
    if(length(ddT2)>1) correlT <- cor(ddT2, ddTE2, method=correlMethod)
}

head(ddT2)
head(ddTE2)

plot(ddT2, ddTE2)

dev.off()


if(plotting==TRUE){
    if(!is.null(ddNE) & !is.null(ddTE)){
        ymax<-max(as.numeric(c(ddNE,ddTE)),na.rm=TRUE)
        ymin<-min(as.numeric(c(ddNE,ddTE)),na.rm=TRUE)
    }
    if(is.null(ddNE)){
        ymax<-max(as.numeric(ddTE),na.rm=TRUE)
        ymin<-min(as.numeric(ddTE),na.rm=TRUE)
    }
    if(is.null(ddTE)){
        ymax<-max(as.numeric(ddNE),na.rm=TRUE)
        ymin<-min(as.numeric(ddNE),na.rm=TRUE)
    }
    if(is.null(ddNE) & is.null(ddTE)){
        ymax<-1
        ymin<-0
    }
    
    if(correlMethod=="spearman") titolcorrel <- "    Spearman rho="
    if(correlMethod=="kendall") titolcorrel <- "    Kendall tau="
    if(correlMethod=="pearson") titolcorrel <- "    Pearson r="
    

    
    
    par(mfrow=c(1,2))
    
    if(!is.null(ddN) & !is.null(ddNE)){
        plot(as.numeric(ddN), as.numeric(ddNE), col="dodgerblue", ylab=paste0("log2(normalized rsem + 1) of ", geneName), xlab=paste0("methylation beta value of ", probeID), xlim=c(0,1), ylim=c(ymin,ymax), las=1)
        if(length(ddN2)>1){
            title(paste0(geneName,"   Normal\n",tissue_label,"\nn=",dim(ddN)[2], titolcorrel, round(correlN,3)))
            if(regressLine) abline(lm(as.numeric(ddNE)~as.numeric(ddN)), col="dodgerblue", lwd=2)
        }
        if(length(ddN2)==1) title(paste0(geneName,"   Normal\n",tissue_label,"\nn=",dim(ddN)[2]))
        
    } else{
        plot(0.5, (ymax+ymin)/2, type = "n", ylab=paste0("log2(normalized rsem + 1) of ", geneName), xlab=paste0("methylation beta value of ", probeID), xlim=c(0,1), ylim=c(ymin,ymax), las=1)
        text(0.5, (ymax+ymin)/2, paste0(missatgeN), col = "dodgerblue")
        title(paste0(geneName,"   Normal\n",tissue_label,"\nn=0"))
        
    }
    box(lwd=2)
    
    if(!is.null(ddT) & !is.null(ddTE)){
        plot(as.numeric(ddT),as.numeric(ddTE),col="darkred",ylab=paste0("log2(normalized rsem + 1) of ", geneName), xlab=paste0("methylation beta value of ", probeID),xlim=c(0,1),ylim=c(ymin,ymax),las=1)
        if(length(ddT2)>1){
            title(paste0(geneName,"   Tumor\n",tissue_label,"\nn=",dim(ddT)[2], titolcorrel, round(correlT,3)))
            if(regressLine) abline(lm(as.numeric(ddTE)~as.numeric(ddT)), col="darkred", lwd=2)
        }
        if(length(ddT2)==1) title(paste0(geneName,"   Tumor\n",tissue_label,"\nn=",dim(ddT)[2]))
    } else{
        plot(0.5, (ymax+ymin)/2, type = "n", ylab=paste0("log2(normalized rsem + 1) of ", geneName), xlab=paste0("methylation beta value of ", probeID), xlim=c(0,1), ylim=c(ymin,ymax), las=1)
        text(0.5, (ymax+ymin)/2, paste0(missatgeT), col = "dodgerblue")
        title(paste0(geneName,"   Tumor\n",tissue_label,"\nn=0"))
        
    }
    box(lwd=2)
    
}
