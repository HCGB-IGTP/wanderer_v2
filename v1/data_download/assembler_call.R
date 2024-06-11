# sessionInfo()
# 
library("HGNChelper")
library("RCurl")
library("httr")
library("stringr")
library("digest")
library("bitops")
# 
# sessionInfo()
setwd("./")
dataDir<-"./data/"
source("./Module_A.r")
source("./Module_B.r")


tissue_types <- read.table("./tissue_types.tab", sep = "\t", header = FALSE)

tiss <- 1:dim(tissue_types)[1]
for(i in tiss){
  print(i)
  
  ##########################################################################################################################
  #RNAseqGene Gene level
  
  dir.create(paste0(dataDir, tissue_types[i,2], "/RNAseqGene/"))
  
  ##################################
  #Tumors
  dir.create(paste0(dataDir, tissue_types[i,2], "/RNAseqGene/Tumor/"))
  Rtumor <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./db/", cancerType = tissue_types[i,2],
                               assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results", tissueType = "TP", outputFileName = paste0(tissue_types[i,2], "_RNAseq_Gene_Tumor"))
  
  
  if(length(Rtumor)>0){
    aux<-grep("hiseq",names(Rtumor))
    gene <- Rtumor[[aux]][,1]
    gene <- gene[-c(1, 2)]
    Rtumor <- as.data.frame(Rtumor[[aux]][,Rtumor[[aux]][2,]=="normalized_count"])
    noms <- as.character(t(Rtumor[1,]))
    Rtumor <- Rtumor[-c(1,2),]
    Rtumor <- data.frame(gene, Rtumor, stringsAsFactors = FALSE)
    colnames(Rtumor) <- c("gene", paste0(noms))
    
    write.table(Rtumor, paste0(dataDir, tissue_types[i,2], "/RNAseqGene/Tumor/", tissue_types[i,2], "_RNAseq_Gene_Tumor.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    create<-c(paste0("CREATE TABLE ", tissue_types[i,2], ".RNAseq.Gene.Tumor(\ngene VARCHAR(27),"), paste0(colnames(Rtumor)[-c(1,dim(Rtumor)[2])], " FLOAT(4),"), paste0(colnames(Rtumor)[dim(Rtumor)[2]], " FLOAT(4));"))
    write.table(create, paste0(dataDir, tissue_types[i,2], "/RNAseqGene/Tumor/", tissue_types[i,2], "_RNAseq_Gene_Tumor_createTable.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  ###########################
  #Normals
  dir.create(paste0(dataDir, tissue_types[i,2], "/RNAseqGene/Normal/"))
  Rnormal<-DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./db/", cancerType = tissue_types[i,2],
                              assayPlatform = "RNASeqV2", dataType = "rsem.genes.normalized_results", tissueType = "NT", outputFileName = paste0(tissue_types[i,2], "_RNAseq_Gene_Normal"))
  
  if(length(Rnormal)>0){
    aux<-grep("hiseq",names(Rnormal))
    gene <- Rnormal[[aux]][,1]
    gene <- gene[-c(1, 2)]
    Rnormal <- as.data.frame(Rnormal[[aux]][,Rnormal[[aux]][2,]=="normalized_count"])
    noms <- as.character(t(Rnormal[1,]))
    Rnormal <- Rnormal[-c(1,2),]
    Rnormal <- data.frame(gene, Rnormal, stringsAsFactors = FALSE)
    colnames(Rnormal) <- c("gene", paste0(noms))
    
    write.table(Rnormal, paste0(dataDir, tissue_types[i,2], "/RNAseqGene/Normal/", tissue_types[i,2], "_RNAseq_Gene_Normal.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    create<-c(paste0("CREATE TABLE ", tissue_types[i,2], ".RNAseq.Gene.Normal(\ngene VARCHAR(27),"), paste0(colnames(Rnormal)[-c(1,dim(Rnormal)[2])], " FLOAT(4),"), paste0(colnames(Rnormal)[dim(Rnormal)[2]], " FLOAT(4));"))
    write.table(create, paste0(dataDir, tissue_types[i,2], "/RNAseqGene/Normal/", tissue_types[i,2], "_RNAseq_Gene_Normal_createTable.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  #######################################################################################################################################
  
  
  
  #######################################################################################################################################
  #RNAseq at exon level
  
  dir.create(paste0(dataDir, tissue_types[i,2], "/RNAseq/"))
  
  #################################
  # RNAseq Tumors
  dir.create(paste0(dataDir, tissue_types[i,2], "/RNAseq/Tumor/"))
  Rtumor <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./db/", cancerType = tissue_types[i,2],
                               assayPlatform = "RNASeqV2", dataType = "exon_quantification", tissueType = "TP", outputFileName = paste0(tissue_types[i,2], "_RNAseq_Tumor"))
  if(length(Rtumor)>0){
    aux<-grep("hiseq",names(Rtumor))
    exon <- Rtumor[[aux]][,1]
    exon <- exon[-c(1, 2)]
    Rtumor <- Rtumor[[aux]][,Rtumor[[aux]][2,]=="RPKM"]
    colnames(Rtumor) <- Rtumor[1,]
    Rtumor <- Rtumor[-c(1,2),]
    Rtumor <- data.frame(exon, Rtumor, stringsAsFactors = FALSE)
    
    write.table(Rtumor, paste0(dataDir, tissue_types[i,2], "/RNAseq/Tumor/", tissue_types[i,2], "_RNAseq_Tumor.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    create<-c(paste0("CREATE TABLE ", tissue_types[i,2], ".RNAseq.Tumor(\nexon VARCHAR(27),"), paste0(colnames(Rtumor)[-c(1,dim(Rtumor)[2])], " FLOAT(4),"), paste0(colnames(Rtumor)[dim(Rtumor)[2]], " FLOAT(4));"))
    write.table(create, paste0(dataDir, tissue_types[i,2], "/RNAseq/Tumor/", tissue_types[i,2], "_RNAseq_Tumor_createTable.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  
  ###########################  ###########################  ###########################
  # RNAseq normals
  dir.create(paste0(dataDir, tissue_types[i,2], "/RNAseq/Normal/"))
  Rnormal<-DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./db/", cancerType = tissue_types[i,2],
                              assayPlatform = "RNASeqV2", dataType = "exon_quantification", tissueType = "NT", outputFileName = paste0(tissue_types[i,2], "_RNAseq_Normal"))
  
  
  if(length(Rnormal)>0){
    aux<-grep("hiseq",names(Rnormal))
    exon <- Rnormal[[aux]][,1]
    exon <- exon[-c(1, 2)]
    Rnormal <- Rnormal[[aux]][,Rnormal[[aux]][2,]=="RPKM"]
    Rnormal<-as.data.frame(Rnormal,stringsAsFactors=FALSE)
    coln <- Rnormal[1,]
    Rnormal <- as.data.frame(Rnormal[-c(1,2),],stringsAsFactors=FALSE)
    colnames(Rnormal) <- coln
    Rnormal <- data.frame(exon, Rnormal, stringsAsFactors = FALSE)
    
    write.table(Rnormal, paste0(dataDir, tissue_types[i,2], "2/RNAseq/Normal/", tissue_types[i,2], "_RNAseq_Normal.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    create<-c(paste0("CREATE TABLE ", tissue_types[i,2], ".RNAseq.Normal(\nexon VARCHAR(27),"), paste0(colnames(Rnormal)[-c(1,dim(Rnormal)[2])], " FLOAT(4),"), paste0(colnames(Rnormal)[dim(Rnormal)[2]], " FLOAT(4));"))
    write.table(create, paste0(dataDir, tissue_types[i,2], "2/RNAseq/Normal/", tissue_types[i,2], "_RNAseq_Normal_createTable.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  #######################################################################################################################################
  #450k methylation arrays
  dir.create(paste0(dataDir, tissue_types[i,2], "/450KMeth/"))
  
  ###########################
  #Meth Tumor
  dir.create(paste0(dataDir, tissue_types[i,2], "/450KMeth/Tumor/"))
  Mtumor<-DownloadMethylationData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./db/", cancerType = tissue_types[i,2],
                                  assayPlatform = "humanmethylation450", tissueType = "TP", outputFileName = paste0(tissue_types[i,2], "_450KMeth_Tumor"))
  if(!is.null(Mtumor)){
    Mtumor <- Mtumor[,c(1,5:dim(Mtumor)[2])]
    colnames(Mtumor) <- Mtumor[1,]
    colnames(Mtumor)[1] <- "probe"
    Mtumor <- Mtumor[-1,]
    
    write.table(Mtumor, paste0(dataDir, tissue_types[i,2], "/450KMeth/Tumor/", tissue_types[i,2], "_450KMeth_Tumor.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    create<-c(paste0("CREATE TABLE ", tissue_types[i,2], ".450KMeth.Tumor(\nprobe VARCHAR(27),"), paste0(colnames(Mtumor)[-c(1,dim(Mtumor)[2])], " FLOAT(4),"), paste0(colnames(Mtumor)[dim(Mtumor)[2]], " FLOAT(4));"))
    write.table(create, paste0(dataDir, tissue_types[i,2], "/450KMeth/Tumor/", tissue_types[i,2], "_450KMeth_Tumor_createTable.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  ###########################
  #Meth Normal
  dir.create(paste0(dataDir, tissue_types[i,2], "/450KMeth/Normal/"))
  Mnormal<-DownloadMethylationData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./db/", cancerType = tissue_types[i,2],
                                   assayPlatform = "humanmethylation450", tissueType = "NT", outputFileName = paste0(tissue_types[i,2], "_450KMeth_Normal"))
  
  if(!is.null(Mnormal)){
    Mnormal <- Mnormal[,c(1,5:dim(Mnormal)[2])]
    colnames(Mnormal) <- Mnormal[1,]
    colnames(Mnormal)[1] <- "probe"
    Mnormal <- Mnormal[-1,]
    
    write.table(Mnormal, paste0(dataDir, tissue_types[i,2], "/450KMeth/Normal/", tissue_types[i,2], "_450KMeth_Normal.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    create<-c(paste0("CREATE TABLE ", tissue_types[i,2], ".450KMeth.Normal(\nprobe VARCHAR(27),"), paste0(colnames(Mnormal)[-c(1,dim(Mnormal)[2])], " FLOAT(4),"), paste0(colnames(Mnormal)[dim(Mnormal)[2]], " FLOAT(4));"))
    write.table(create, paste0(dataDir, tissue_types[i,2], "/450KMeth/Normal/", tissue_types[i,2], "_450KMeth_Normal_createTable.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  ##########################################################################################################################
  #Clinical Data
  dir.create(paste0(dataDir, tissue_types[i,2], "/clinical/"))
  Clinical <- DownloadClinicalData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./db/Clinical/", cancerType = tissue_types[i,2],
                                   clinicalDataType = "patient", outputFileName = paste0(tissue_types[i,2], "_Clinical"))
  
  if(length(Clinical)>0){
    gene <- Clinical[[1]][,1]
    gene <- gene[-c(1, 2)]
    Clinical <- as.data.frame(Clinical[[1]][,Clinical[[1]][2,]=="normalized_count"])
    noms <- Clinical[1,]
    Clinical <- Clinical[-c(1,2),]
    Clinical <- data.frame(gene, Clinical, stringsAsFactors = FALSE)
    colnames(Clinical) <- c("gene", paste0(noms))
    
    write.table(Clinical, paste0(dataDir, tissue_types[i,2], "/RNAseqGene/Tumor/", tissue_types[i,2], "_RNAseq_Gene_Tumor.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    create<-c(paste0("CREATE TABLE ", tissue_types[i,2], ".RNAseq.Gene.Tumor(\ngene VARCHAR(27),"), paste0(colnames(Rtumor)[-c(1,dim(Rtumor)[2])], " FLOAT(4),"), paste0(colnames(Rtumor)[dim(Rtumor)[2]], " FLOAT(4));"))
    write.table(create, paste0(dataDir, tissue_types[i,2], "/RNAseqGene/Tumor/", tissue_types[i,2], "_RNAseq_Gene_Tumor_createTable.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  print(i)
}
