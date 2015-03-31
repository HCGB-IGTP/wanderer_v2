
library("HGNChelper")
library("RCurl")
library("httr")
library("stringr")
library("digest")
library("bitops")

setwd("/imppc/labs/maplab/adiez/region_profile/TCGA_assembler/TCGA-Assembler/")
dataDir<-"/imppc/labs/maplab/share/anna2izaskun/db_region_profile_data/"
source("./Module_A.r")
source("./Module_B.r")


tissue_types <- read.table("/imppc/labs/maplab/adiez/region_profile/web/tissue_types.tab", sep = "\t", header = FALSE)

for(i in 1:length(tissue_types[,2])){
  print(i)
  dir.create(paste0(, tissue_types[i,2], "/"))
  
  dir.create(paste0(dataDir, tissue_types[i,2], "/RNAseq/"))
  
  ##################################
  #RNAseq tumors
  dir.create(paste0(dataDir, tissue_types[i,2], "/RNAseq/Tumor/"))
  Rtumor <- DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./db/", cancerType = tissue_types[i,2],
                               assayPlatform = "RNASeqV2", dataType = "exon_quantification", tissueType = "TP", outputFileName = paste0(tissue_types[i,2], "_RNAseq_Tumor"))
  if(length(Rtumor)>0){
    exon <- Rtumor[[1]][,1]
    exon <- exon[-c(1, 2)]
    Rtumor <- Rtumor[[1]][,Rtumor[[1]][2,]=="RPKM"]
    colnames(Rtumor) <- Rtumor[1,]
    Rtumor <- Rtumor[-c(1,2),]
    Rtumor <- data.frame(exon, Rtumor, stringsAsFactors = FALSE)
    
    write.table(Rtumor, paste0(dataDir, tissue_types[i,2], "/RNAseq/Tumor/", tissue_types[i,2], "_RNAseq_Tumor.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    create<-c(paste0("CREATE TABLE ", tissue_types[i,2], ".RNAseq.Tumor(\nexon VARCHAR(27),"), paste0(colnames(Rtumor)[-c(1,dim(Rtumor)[2])], " FLOAT(4),"), paste0(colnames(Rtumor)[dim(Rtumor)[2]], " FLOAT(4));"))
    write.table(create, paste0(dataDir, tissue_types[i,2], "/RNAseq/Tumor/", tissue_types[i,2], "_RNAseq_Tumor_createTable.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  ###########################
  #RNAseq normals
  dir.create(paste0(dataDir, tissue_types[i,2], "/RNAseq/Normal/"))
  Rnormal<-DownloadRNASeqData(traverseResultFile = "./DirectoryTraverseResult_Jul-08-2014.rda", saveFolderName = "./db/", cancerType = tissue_types[i,2],
                              assayPlatform = "RNASeqV2", dataType = "exon_quantification", tissueType = "NT", outputFileName = paste0(tissue_types[i,2], "_RNAseq_Normal"))
  
  if(length(Rnormal)>0){
    exon <- Rnormal[[1]][,1]
    exon <- exon[-c(1, 2)]
    Rnormal <- Rnormal[[1]][,Rnormal[[1]][2,]=="RPKM"]
    colnames(Rnormal) <- Rnormal[1,]
    Rnormal <- Rnormal[-c(1,2),]
    Rnormal <- data.frame(exon, Rnormal, stringsAsFactors = FALSE)
    
    write.table(Rnormal, paste0(dataDir, tissue_types[i,2], "/RNAseq/Normal/", tissue_types[i,2], "_RNAseq_Normal.csv"), sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
    create<-c(paste0("CREATE TABLE ", tissue_types[i,2], ".RNAseq.Normal(\nexon VARCHAR(27),"), paste0(colnames(Rnormal)[-c(1,dim(Rnormal)[2])], " FLOAT(4),"), paste0(colnames(Rnormal)[dim(Rnormal)[2]], " FLOAT(4));"))
    write.table(create, paste0(dataDir, tissue_types[i,2], "/RNAseq/Normal/", tissue_types[i,2], "_RNAseq_Normal_createTable.txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  #######################################################################################################################################
  dir.create(paste0(dataDir, tissue_types[i,2], "/450KMeth/"))
  
  ###########################
  #Meth tumors
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
  print(i)
}
