

correl_meth_express <- function(geneName, probeID, ddmeth, ddGene, tissue_label, regressLine, correlMethod, plotting, datareturn){
  
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
  
  ##############################################
  #common patients in normal
  if(!is.null(ddN) & !is.null(ddNE)){
    pNE <- sapply(strsplit(colnames(ddNE),"-"),"[",3)
    pN <- sapply(strsplit(colnames(ddN),"-"),"[",3)
    
    colN <- colnames(ddN)[pN%in%pNE]
    colNE <- colnames(ddNE)[pNE%in%pN]
    ddN <- as.data.frame(ddN[,pN%in%pNE])
    colnames(ddN) <- colN
    ddNE <- as.data.frame(ddNE[,pNE%in%pN])
    colnames(ddNE) <- colNE
    if(length(ddN)==0 | length(ddNE)==0){
      ddN <- ddNE <- NULL
      missatgeN <- "No common patients available"
    } else{
      pNE <- sapply(strsplit(colnames(ddNE),"-"),"[",3)
      pN <- sapply(strsplit(colnames(ddN),"-"),"[",3)
      if(length(pN)!=length(pNE)){
        numpatients <- c(length(pN),length(pNE))
        aux <- which(numpatients==max(numpatients))
        if(aux==1){
          ddnew <- as.data.frame(matrix(0,ncol=numpatients[aux],nrow=dim(ddNE)[1]))
          for(i in 1:numpatients[aux]){
            aux2 <- which(pNE==pN[i])
            ddnew[,i] <- ddNE[,aux2]
            colnames(ddnew)[i] <- colnames(ddNE)[aux2]
          }
          ddNE <- ddnew
        }
        if(aux==2){
          ddnew <- as.data.frame(matrix(0,ncol=numpatients[aux],nrow=dim(ddN)[1]))
          for(i in 1:numpatients[aux]){
            aux2 <- which(pN==pNE[i])
            ddnew[,i] <- ddN[,aux2]
            colnames(ddnew)[i] <- colnames(ddN)[aux2]
          }
          ddN <- ddnew
        }
      }
    }
    if(!is.null(ddN) & !is.null(ddNE)){
      ddN2<-as.numeric(ddN)
      ddNE2<-as.numeric(ddNE)
      ddNE2<-ddNE2[!is.na(ddN2)]
      ddN2<-ddN2[!is.na(ddN2)]
      ddN2<-ddN2[!is.na(ddNE2)]
      ddNE2<-ddNE2[!is.na(ddNE2)]
      if(length(ddN2)>1)  correlN <- cor(ddN2, ddNE2, method=correlMethod)
    }
  }   
  
  
  
  #######################################
  #common patients in tumor
  if(!is.null(ddT) & !is.null(ddTE)){
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
  } 
  
  
  #plot
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
  
  if(!is.null(ddNE)){
    options(warn=-1)
    ddNE2 <- data.frame(PATIENT=rownames(t(ddNE)), t(ddNE))
    options(warn=0)
    colnames(ddNE2)[2] <- geneName
    ddNE2 <- ddNE2[!duplicated(ddNE2),]
  } else{
    ddNE2 <- NULL
  }
  
  if(!is.null(ddTE)){
    options(warn=-1)
    ddTE2 <- data.frame(PATIENT=rownames(t(ddTE)), t(ddTE))
    options(warn=0)
    colnames(ddTE2)[2] <- geneName
    ddTE2 <- ddTE2[!duplicated(ddTE2),]
  } else{
    ddTE2 <- NULL
  }
  
  if(datareturn)  return(list(ddNE=ddNE2, ddTE=ddTE2))
}
