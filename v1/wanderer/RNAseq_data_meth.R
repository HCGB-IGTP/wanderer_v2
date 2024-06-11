
RNAseq_data_meth <- function(ddGene, dd){
  
  rownames(dd) <- dd[,1]
  dd <- dd[,-1]
  
  rownames(ddGene) <- ddGene[,1]
  ddGene <- ddGene[,-1]
  
  #common patients
  pNE <- sapply(strsplit(colnames(ddGene),"-"),"[",3)
  pN <- sapply(strsplit(colnames(dd),"-"),"[",3)
  dd <- dd[,pN%in%pNE]
  ddGene<-ddGene[,pNE%in%pN]
  
  ddnew<-as.data.frame(matrix(0,ncol=1,nrow=dim(dd)[2]))
  for(i in 1:numpatients[aux]){
    aux2<-which(pNE==pN[i])
    ddnew[,i]<-ddGene[,aux2]
    colnames(ddnew)[i]<-colnames(ddGene)[aux2]
  }
  rownames(ddnew)<-rownames(ddGene)
  ddGene<-ddnew
}
if(aux==2){
  ddnew<-as.data.frame(matrix(0,ncol=numpatients[aux],nrow=dim(dd)[1]))
  for(i in 1:numpatients[aux]){
    aux2<-which(pN==pNE[i])
    ddnew[,i]<-dd[,aux2]
    colnames(ddnew)[i]<-colnames(dd)[aux2]
  }
  rownames(ddnew)<-rownames(dd)
  dd<-ddnew
}
}

aux1<-order(rownames(dd))
dd<-dd[aux1,]
ddT<-ddT[aux1,]

aux1<-order(rownames(ddGene))
ddGene<-ddGene[aux1,]
ddTE<-ddTE[aux1,]















pdd<-sapply(strsplit(colnames(dd)[-1],"-"),"[",3)
pGene<-sapply(strsplit(colnames(ddGene)[-1],"-"),"[",3)
aux<-pdd%in%pGene

if(length(aux)==0) dataGene<-NULL

if(length(aux)>0) {  
  namsG<-colnames(ddGene)[-1]
  namsdd<-colnames(dd)[-1]
  
  order_columns <- match(pGene[aux], pdd)
  nams<-nams[order_columns]
  dataGene <- as.data.frame(ddGene[,-1])
  dataGene <- data.frame(ddGene[,1],dataGene[,order_columns])
  colnames(dataGene)<-c("GENE",nams)
  for(i in 1:6) colnames(dataGene) <- sub("\\.","-",colnames(dataGene))
  dataGene <- t(dataGene)
  dataGene <- data.frame(rownames(dataGene), dataGene[,1])
  colnames(dataGene) <- c("PATIENT",paste0(dataGene[1,2]))
  dataGene <- dataGene[-1,]
}
return(dataGene)
}




