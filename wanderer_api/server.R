#!/usr/bin/env R

library(shiny)
library(RPostgreSQL)
library(Cairo)

SRC <- '.'

DB_CONF <- file.path(SRC, 'db.txt')

source(file.path(SRC, 'GeneSize.R'))
source(file.path(SRC, 'expression_data.R'))
source(file.path(SRC, 'data_RNAseqGene.R'))
source(file.path(SRC, 'RNAseq_data_table.R'))
source(file.path(SRC, 'plot_RNAseqGene.R'))
source(file.path(SRC, 'correl_meth_express.R'))
source(file.path(SRC, 'methylation_data.R'))
source(file.path(SRC, 'wanderer_expression.R'))
source(file.path(SRC, 'wanderer_methylation.R'))
source(file.path(SRC, 'database.R'))
source(file.path(SRC, 'max_sample.R'))
source(file.path(SRC, 'data_meth_filtering.R'))
source(file.path(SRC, 'data_expr_filtering.R'))
source(file.path(SRC, 'stat_analysis_meth.R'))
source(file.path(SRC, 'stat_analysis_expr.R'))
source(file.path(SRC, 'permalink_modal.R'))
source(file.path(SRC, 'utils.R'))


sample_size <- read.table(file.path(SRC, "NumberOfSamples.csv"), sep = ",", stringsAsFactors = FALSE, header = TRUE)

error <- 'Malformed query'

shinyServer(function(input, output, session){
  
  #database connection
  con <- db_connect(DB_CONF)

  query<- list()

  output$urlText <- renderText({
    paste(sep = "",
      "protocol: ", session$clientData$url_protocol, "\n",
      "hostname: ", session$clientData$url_hostname, "\n",
      "pathname: ", session$clientData$url_pathname, "\n",
      "port: ",     session$clientData$url_port,     "\n",
      "search: ",   session$clientData$url_search,   "\n"
    )
  })

  output$queryText <- renderText({
    query <<- parseQueryString(session$clientData$url_search)

    # Return a string with key-value pairs
    stop(paste(names(query), query, sep = "=", collapse=", "))
  })

  query <- isolate(parseQueryString(session$clientData$url_search))

  ## print(paste(names(query), query, sep = "=", collapse=", "))

  ## test/?Gene=BRCA1&start=41198000&end=41278000&TissueType=brca&goButton=1&DataType=methylation&plotmean=FALSE&geneLine=TRUE&CpGi=FALSE&nN=30&nT=10"
  
  ## parameter processing
  ## query$Gene <- toupper(query$gene)
  query$Zoom <- as.numeric(c(query$start, query$end))
  query$start <- as.numeric(query$start)
  query$end <- as.numeric(query$end)
  query$nN <- as.numeric(query$nN)
  query$nT <- as.numeric(query$nT)
  query$plotmean <- as.logical(query$plotmean)
  query$geneLine <- as.logical(query$geneLine)
  query$CpGi <- as.logical(query$CpGi)
  ## query$region <- as.logical(query$region)
  query$region <- TRUE
  query$goButton <- 1
  query$distribute_uniformly <- as.logical(query$distribute_uniformly)
  query$pvalThres <- as.numeric(query$pvalThres)

  api_arguments_allowances <- sort(c('Gene', 'start','end', 'TissueType', 'DataType', 'plotmean',
                                     'geneLine', 'CpGi', 'nN', 'nT', 'region', 'goButton', 'Zoom',
                                     'distribute_uniformly', 'pvalThres'))
  
  
  #################################################
  #detect gene format
  geneFormat <- reactive({
      if(query$goButton == 0) {
          GeneFormat <- "genename"
      } else{
          if(substr(toupper(query$Gene), 1, 4) == "ENSG"){
              GeneFormat <- "emsemblgeneid"
          } else{
              GeneFormat <- "genename"
          }
      }
  })

  
  #################################################
  #Filtering Methylation data
  geneSize <- reactive({
      GeneSize(con = con, geneName =toupper(query$Gene), geneNamesType = geneFormat())
  })
   
  
  #################################################
  #Reading Methylation data
  datameth <- reactive({
    if(query$DataType == 'methylation'){
        if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1)
            methylation_data(con = con, geneName = toupper(query$Gene), geneNamesType = geneFormat(), tissue = query$TissueType)
    }
  })
  
  #################################################
  #Filtering Methylation data
  datamethfilt <- reactive({
    if(query$DataType == 'methylation'){
      if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
        if(!datameth()[['empty']])
            data_meth_filtering(results = datameth(), sampleSize = sample_size, tissue = query$TissueType, zoom = c(query$start, query$end))
      }
    }
  })
  
  
  #################################################
  #Reading Expression data
  dataexpr <- reactive({
    if(query$DataType == 'expression'){
        if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1)
            expression_data(con = con, geneName = toupper(query$Gene), geneNamesType = geneFormat(), tissue = query$TissueType)
      }
  })
  
  #################################################
  #Filtering Expression data
  dataexprfilt <- reactive({
    if(query$DataType == 'expression'){
      if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
        if(!dataexpr()[['empty']])
            data_expr_filtering(results = dataexpr(), sampleSize = sample_size, tissue = query$TissueType, zoom = c(query$start, query$end))    
      }
    }
  })
  
  
  #################################################
  #print wanderer plot
  output$plot1 <- renderPlot({
    if(!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(geneNameSaved()) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      if(input$region & ((input$end > input$Zoom[2]) | (input$end < input$Zoom[1]) | (input$start < input$Zoom[1]) | (input$start > input$Zoom[2]))) stop(print(paste0("The region must be between ", input$Zoom[1], " and ", input$Zoom[2])))
      
      if(input$DataType == 'methylation'){
        if(dim(datamethfilt()[['probes2']])[1]>0){
          if(is.null(datamethfilt()$ddN2) & is.null(datamethfilt()$ddT2)){
            stop("There are no samples in this tissue type")
          } else{
            wanderer_methylation(results_filt = datamethfilt(), geneName = geneNameSaved(),
                                 geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                 CpGislands = input$CpGi, plotmean = input$plotmean,
                                 plotting = TRUE, geneLine = input$geneLine, proportional = !(input$distribute_uniformly))
          }
        }
      }
      if(input$DataType == 'expression'){
        if(dim(dataexprfilt()[['exons2']])[1]>0){
          if(is.null(dataexprfilt()$ddN2) & is.null(dataexprfilt()$ddT2)){
            stop("There are no samples in this tissue type")
          } else{
            wanderer_expression(results_filt = dataexprfilt(), geneName = geneNameSaved(),
                                geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine, proportional = !(input$distribute_uniformly))
          }
        }
      }
    }
  }, height = 1000, width = 1000)
  
  
  ##################################################
  #print summary plot
  output$plotStat <- renderPlot({
    if(!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(geneNameSaved()) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      
      if(input$DataType == 'methylation'){
        if(dim(datamethfilt()[['probes2']])[1]>0){
          if(is.null(datamethfilt()$ddN2) | is.null(datamethfilt()$ddT2)){
            stop("There are not enough samples to perform the statistical analysis")
          } else{
            
            if(dim(datamethfilt()$ddN2)[2]<=2 | dim(datamethfilt()$ddT2)[2]<=2){
              stop("There are not enough samples to perform the statistical analysis")
            } else{
              stat_analysis_meth(results_filt = datamethfilt(), geneName = geneNameSaved(),
                                 geneNamesType = geneFormat(), CpGislands = input$CpGi, pvalThres = input$pvalThres,
                                 geneLine = input$geneLine, plotting = TRUE, proportional = !(input$distribute_uniformly))
            }
          }
        }
      }
      else if(input$DataType == 'expression'){
        if(dim(dataexprfilt()[['exons2']])[1]>0){
          if(is.null(dataexprfilt()$ddN2) | is.null(dataexprfilt()$ddT2)){
            stop("There are not enough samples to perform the statistical analysis")
          } else{
            
            if(dim(dataexprfilt()$ddN2)[2]<=2 | dim(dataexprfilt()$ddT2)[2]<=2){
              stop("There are not enough samples to perform the statistical analysis")
            } else{
              stat_analysis_expr(results_filt = dataexprfilt(), geneName = geneNameSaved(),
                                 geneNamesType = geneFormat(), pvalThres = input$pvalThres,
                                 geneLine = input$geneLine, plotting = TRUE, proportional = !(input$distribute_uniformly))
            }
          }
        }
      }
    }
  }, height = 500, width = 1000)
  
  #################################################
  #index html stat parameters to print (adjusted pval)
  output$pvalParam <- reactive({
    pvalparam <- 0
    if(!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(geneNameSaved()) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      if(input$DataType == 'methylation'){
        if(dim(datamethfilt()[['probes2']])[1]>0){
          if(is.null(datamethfilt()$ddN2) | is.null(datamethfilt()$ddT2)){
            pvalparam <- 0
          } else{
            if(dim(datamethfilt()$ddN2)[2]<=2 | dim(datamethfilt()$ddT2)[2]<=2){
              pvalparam <- 0
            } else{
              pvalparam <- 1
            }
          }
        }
      } else if(input$DataType == 'expression'){
        if(dim(dataexprfilt()[['exons2']])[1]>0){
          if(is.null(dataexprfilt()$ddN2) | is.null(dataexprfilt()$ddT2)){
            pvalparam <- 0
          } else{
            if(dim(dataexprfilt()$ddN2)[2]<=2 | dim(dataexprfilt()$ddT2)[2]<=2){
              pvalparam <- 0
            } else{
              pvalparam <- 1
            }
          }
        }
      }
    }
    return(pvalparam)
  })
  
  #################################################
  #index html stat parameters to print (adjusted pval)
  output$downloadParam <- reactive({
    downloadparam <- 0
    if(!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(geneNameSaved()) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      if(input$DataType == 'expression'){
        if(dim(dataexprfilt()[['exons2']])[1]>0){
          if(is.null(dataexprfilt()$ddN2) & is.null(dataexprfilt()$ddT2)){
            downloadparam <- 0
          } else{
            downloadparam <- 1
          }
        }
      } else downloadparam <- 1
    }
    return(downloadparam)
  })
  
  #################################################
  #print the number of probes or exons
  output$numberpoints <- renderText({
    if(!is.null(input$TissueType) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
      if(input$DataType == 'methylation'){
        if(!is.null(datamethfilt())){
          if(dim(datamethfilt()[['probes2']])[1]==0) printa <- "There are not probes in this region"
          if(dim(datamethfilt()[['probes2']])[1]>0)  printa <- paste0("There are ", dim(datamethfilt()[['probes2']])[1] ," probes in the selected region")
        }
      }
      if(input$DataType == 'expression'){
        if(!is.null(dataexprfilt())){
          if(dim(dataexprfilt()[['exons2']])[1]==0) printa <- "There are not exons in this region"
          if(dim(dataexprfilt()[['exons2']])[1]>0)  printa <- paste0("There are ", dim(dataexprfilt()[['exons2']])[1] ," exons in the selected region")
        }
      }
      return(printa)
    }
  })
  
  
  
  #################################################
  #DOWNLOAD RESULTS
  
  output$downloadResults <- downloadHandler(
    
    filename = function(){
      file.path(paste0("Wanderer_", geneNameSaved(), '_', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.zip'))
    },
    
    content = function(file) {
      
      ######################################
      #wanderer plot in png      
      f1 <- paste0("1_Wanderer_", geneNameSaved(), '_', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.png')
      CairoPNG(file.path(tempdir(),f1), width = 1000, height = 1000)
      if(input$DataType == 'methylation'){
        regplot <- wanderer_methylation(results_filt = datamethfilt(), geneName = geneNameSaved(),
                                        geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                        CpGislands = input$CpGi, plotmean = input$plotmean,
                                        plotting = TRUE, geneLine = input$geneLine, proportional = !(input$distribute_uniformly))
      } else if(input$DataType == 'expression'){
        regplot <- wanderer_expression(results_filt = dataexprfilt(), geneName = geneNameSaved(),
                                       geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                       plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine, proportional = !(input$distribute_uniformly))
      }
      dev.off()
      
      #####################################
      #wanderer plot in pdf
      f2 <- paste0("2_Wanderer_", geneNameSaved(), '_', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.pdf')
      pdf(file.path(tempdir(), f2), width = 14, height = 14)
      if(input$DataType == 'methylation'){
        regplot <- wanderer_methylation(results_filt = datamethfilt(), geneName = geneNameSaved(),
                                        geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                        CpGislands = input$CpGi, plotmean = input$plotmean,
                                        plotting = TRUE, geneLine = input$geneLine, proportional = !(input$distribute_uniformly))
      } else if(input$DataType == 'expression'){
        regplot <- wanderer_expression(results_filt = dataexprfilt(), geneName = geneNameSaved(),
                                       geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                       plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine, proportional = !(input$distribute_uniformly))
      }
      dev.off()
      
      
      #################################################
      #dowload Normal data
      fileN <- paste0("3_Wanderer_", geneNameSaved(), '_', input$DataType, '_', input$TissueType, '_Normal_', client_timestamp(input$clientTime), '.csv')
      if(input$DataType == 'methylation'){
        if(!is.null(datamethfilt()$ddN2)) write.table(datamethfilt()$ddN2, file = file.path(tempdir(), fileN), sep = ",", row.names = FALSE, quote = FALSE)
      } else if(input$DataType == 'expression'){
        if(!is.null(dataexprfilt()$ddN2)) write.table(dataexprfilt()$ddN2, file = file.path(tempdir(), fileN), sep = ",", row.names = FALSE, quote = FALSE)
      }
      
      
      #################################################
      #dowload RNAseq Gene level Normal data
      if(input$DataType == 'expression'){
        fileNG <- paste0("10_Wanderer_", geneNameSaved(), '_', input$DataType, '_', input$TissueType, '_Normal_RNAseqGENE_', client_timestamp(input$clientTime), '.csv')
        if(!is.null(dataRNAseqGene()$Normal)){
          ddN <- RNAseq_data_table(dataRNAseqGene()$Normal, dataexprfilt()$ddN2)
          write.table(ddN, file = file.path(tempdir(), fileNG), sep = ",", row.names = FALSE, quote = FALSE)        
        }
      }
      
      #################################################
      #dowload RNAseq Gene level Normal data common with methylation
      if(input$DataType == 'methylation'){
        fileNG <- paste0("10_Wanderer_", geneNameSaved(), '_', input$DataType, '_', input$TissueType, '_Normal_RNAseqGENE_CommonWithMethylation_', client_timestamp(input$clientTime), '.csv')
        if(!is.null(dataRNAseqGene()[['Normal']]) & !is.null(datamethfilt()[['ddN2']])){        
          RNAseqdata <- correl_meth_express(geneName = geneNameSaved(), probeID = input$ProbeSelection, ddmeth = datamethfilt(), ddGene = dataRNAseqGene(),  tissue_label = datamethfilt()[['tissue_label']], regressLine = input$regressionLine, correlMethod = input$correlationMethod, plotting = FALSE, datareturn = TRUE)
          if(!is.null(RNAseqdata$ddNE)) write.table(RNAseqdata$ddNE, file = file.path(tempdir(), fileNG), sep = ",", row.names = FALSE, quote = FALSE)        
        }
      }
      
      #################################################
      #dowload Tumor data
      fileT <- paste0("4_Wanderer_", geneNameSaved(), '_', input$DataType, '_', input$TissueType, '_Tumor_', client_timestamp(input$clientTime), '.csv')
      if(input$DataType == 'methylation'){
        if(!is.null(datamethfilt()$ddT2)) write.table(datamethfilt()$ddT2, file = file.path(tempdir(), fileT), sep = ",", row.names = FALSE, quote = FALSE)
      } else if(input$DataType == 'expression'){
        if(!is.null(dataexprfilt()$ddT2)) write.table(dataexprfilt()$ddT2, file = file.path(tempdir(), fileT), sep = ",", row.names = FALSE, quote = FALSE)
      }
      
      #################################################
      #dowload RNAseq Gene level Tumor data
      if(input$DataType == 'expression'){
        fileTG <- paste0("11_Wanderer_", geneNameSaved(), '_', input$DataType, '_', input$TissueType, '_Tumor_RNAseqGENE_', client_timestamp(input$clientTime), '.csv')
        if(!is.null(dataRNAseqGene()$Tumor)){
          ddT <- RNAseq_data_table(dataRNAseqGene()$Tumor, dataexprfilt()$ddT2)
          write.table(ddT, file = file.path(tempdir(), fileTG), sep = ",", row.names = FALSE, quote = FALSE)        
        }
      }
      
      #################################################
      #dowload RNAseq Gene level Tumor data common with methylation
      if(input$DataType == 'methylation'){
        fileTG <- paste0("11_Wanderer_", geneNameSaved(), '_', input$DataType, '_', input$TissueType, '_Tumor_RNAseqGENE_CommonWithMethylation_', client_timestamp(input$clientTime), '.csv')
        if(!is.null(dataRNAseqGene()[['Tumor']]) & !is.null(datamethfilt()[['ddT2']])){        
          RNAseqdata <- correl_meth_express(geneName = geneNameSaved(), probeID = input$ProbeSelection, ddmeth = datamethfilt(), ddGene = dataRNAseqGene(),  tissue_label = datamethfilt()[['tissue_label']], regressLine = input$regressionLine, correlMethod = input$correlationMethod, plotting = FALSE, datareturn = TRUE)
          if(!is.null(RNAseqdata$ddTE)) write.table(RNAseqdata$ddTE, file = file.path(tempdir(), fileTG), sep = ",", row.names = FALSE, quote = FALSE)        
        }
      }
      
      #################################################
      #dowload probe annotation and statistical analysis
      fileA <- paste0("7_Wanderer_", geneNameSaved(), '_', input$DataType, '_', input$TissueType, '_annotations_and_statistical_analysis_', client_timestamp(input$clientTime), '.csv')
      if(input$DataType == 'methylation'){
        results <-stat_analysis_meth(results_filt = datamethfilt(), geneName = geneNameSaved(),
                                     geneNamesType = geneFormat(), CpGislands = input$CpGi,
                                     geneLine = input$geneLine, plotting = FALSE, proportional = !(input$distribute_uniformly) )
      } else if(input$DataType == 'expression'){
        results <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = geneNameSaved(),
                                      geneNamesType = geneFormat(),
                                      geneLine = input$geneLine, plotting = FALSE, proportional = !(input$distribute_uniformly))
      }
      write.table(results, file = file.path(tempdir(), fileA), sep = ",", row.names = FALSE, quote = FALSE)        
      
      
      #################################################
      #dowload boxplot of RNAseq gene data as png
      if(input$DataType == 'expression'){
        fbox1 <- paste0("8_Wanderer_", geneNameSaved(), '_boxplot_', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.png')
        if(!is.null(dataRNAseqGene())){
          CairoPNG(file.path(tempdir(), fbox1), width = 1000, height = 500)
          RNAseqboxplot <- plot_RNAseqGene(dd = dataRNAseqGene(), geneName = geneNameSaved(), tissue_label = dataexprfilt()[['tissue_label']]) 
          print(RNAseqboxplot)
          dev.off()
        } 
      }
      
      #################################################
      #dowload boxplot of RNAseq gene data as pdf
      if(input$DataType == 'expression'){
        fbox2 <- paste0("9_Wanderer_", geneNameSaved(), '_boxplot_', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.pdf')
        if(!is.null(dataRNAseqGene())){
          pdf(file.path(tempdir(),fbox2), width = 12, height = 6)
          RNAseqboxplot <- plot_RNAseqGene(dd = dataRNAseqGene(), geneName = geneNameSaved(), tissue_label = dataexprfilt()[['tissue_label']]) 
          dev.off()
        } 
      }
      
      #################################################
      #dowload correl of meth vs RNAseq gene data as png
      if(input$DataType == 'methylation'){
        fbox1 <- paste0("8_Wanderer_", geneNameSaved(), '_correl_RNAseqGeneVS', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.png')
        if((!is.null(dataRNAseqGene()[['Normal']]) & !is.null(datamethfilt()[['ddN2']])) | (!is.null(dataRNAseqGene()[['Tumor']]) & !is.null(datamethfilt()[['ddT2']]))){
          CairoPNG(file.path(tempdir(), fbox1), width = 1000, height = 500)
          RNAseqcorrelplot <- correl_meth_express(geneName = geneNameSaved(), probeID = input$ProbeSelection, ddmeth = datamethfilt(), ddGene = dataRNAseqGene(),  tissue_label = datamethfilt()[['tissue_label']], regressLine = input$regressionLine, correlMethod = input$correlationMethod, plotting = TRUE, datareturn = FALSE) 
          print(RNAseqcorrelplot)
          dev.off()
        } 
      }
      
      #################################################
      #dowload correl of meth vs RNAseq gene data as pdf
      if(input$DataType == 'methylation'){
        fbox2 <- paste0("9_Wanderer_", geneNameSaved(), '_correl_RNAseqGeneVS', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.pdf')
        if((!is.null(dataRNAseqGene()[['Normal']]) & !is.null(datamethfilt()[['ddN2']])) | (!is.null(dataRNAseqGene()[['Tumor']]) & !is.null(datamethfilt()[['ddT2']]))){
          pdf(file.path(tempdir(),fbox2), width = 12, height = 6)
          RNAseqcorrelplot <- correl_meth_express(geneName = geneNameSaved(), probeID = input$ProbeSelection, ddmeth = datamethfilt(), ddGene = dataRNAseqGene(),  tissue_label = datamethfilt()[['tissue_label']], regressLine = input$regressionLine, correlMethod = input$correlationMethod, plotting = TRUE, datareturn = FALSE) 
          dev.off()
        } 
      }
      
      #################################################
      #dowload mean plot as png
      if(input$DataType == 'methylation'){
        if(!is.null(datamethfilt()$ddN2) & !is.null(datamethfilt()$ddT2)){
          if(dim(datamethfilt()$ddN2)[2]>2 & dim(datamethfilt()$ddT2)[2]>2){   
            fmean1 <- paste0("5_Wanderer_", geneNameSaved(), '_Mean_', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.png')
            CairoPNG(file.path(tempdir(), fmean1), width = 1000, height = 500)
            regplot <- stat_analysis_meth(results_filt = datamethfilt(), geneName = geneNameSaved(),
                                          geneNamesType = geneFormat(), CpGislands = input$CpGi, pvalThres = input$pvalThres,
                                          geneLine = input$geneLine, plotting = TRUE, proportional = !(input$distribute_uniformly))
            dev.off()
          } else fmean1 <- NULL
        } else fmean1 <- NULL
      } else if(input$DataType == 'expression'){
        if(!is.null(dataexprfilt()$ddN2) & !is.null(dataexprfilt()$ddT2)){
          if(dim(dataexprfilt()$ddN2)[2]>2 & dim(dataexprfilt()$ddT2)[2]>2){   
            fmean1 <- paste0("5_Wanderer_", geneNameSaved(), '_Mean_', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.png')
            CairoPNG(file.path(tempdir(), fmean1), width = 1000, height = 500)
            regplot <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = geneNameSaved(),
                                          geneNamesType = geneFormat(), pvalThres = input$pvalThres,
                                          geneLine = input$geneLine, plotting = TRUE, proportional = !(input$distribute_uniformly))
            dev.off()
          } else fmean1 <- NULL
        } else fmean1 <- NULL
      }
      
      
      
      #################################################
      #dowload mean plot as pdf
      if(input$DataType == 'methylation'){
        if(!is.null(datamethfilt()$ddN2) & !is.null(datamethfilt()$ddT2)){
          if(dim(datamethfilt()$ddN2)[2]>2 & dim(datamethfilt()$ddT2)[2]>2){   
            fmean2 <- paste0("6_Wanderer_", geneNameSaved(), '_Mean_', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.pdf')
            pdf(file.path(tempdir(), fmean2), width = 14, height = 6)
            regplot <- stat_analysis_meth(results_filt = datamethfilt(), geneName = geneNameSaved(),
                                          geneNamesType = geneFormat(), CpGislands = input$CpGi, pvalThres = input$pvalThres,
                                          geneLine = input$geneLine, plotting = TRUE, proportional = !(input$distribute_uniformly))
            dev.off()
          } else fmean2 <- NULL
        } else fmean2 <- NULL
      } else if(input$DataType == 'expression'){
        if(!is.null(dataexprfilt()$ddN2) & !is.null(dataexprfilt()$ddT2)){
          if(dim(dataexprfilt()$ddN2)[2]>2 & dim(dataexprfilt()$ddT2)[2]>2){   
            fmean2 <- paste0("6_Wanderer_", geneNameSaved(), '_Mean_', input$DataType, '_', input$TissueType, '_', client_timestamp(input$clientTime), '.pdf')
            pdf(file.path(tempdir(), fmean2), width = 14, height = 6)
            regplot <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = geneNameSaved(),
                                          geneNamesType = geneFormat(), pvalThres = input$pvalThres,
                                          geneLine = input$geneLine, plotting = TRUE, proportional = !(input$distribute_uniformly))
            dev.off()
          } else fmean2 <- NULL
        } else fmean2<- NULL
      }
      
      
      #################################################
      #dowload documentation
      fileDoc <- file.path(SRC,"Wanderer_Documentation.txt")
      fileDoc2 <- file.path(SRC, "Wanderer_Documentation.pdf")
      
      #################################################
      #dowload documentation

      fileClinic <- file.path(SRC, 'Clinical', paste0(toupper(input$TissueType),"_Clinical__nationwidechildrens.org_clinical_patient_",input$TissueType,".txt"))

      
      #################################################
      
      zip(zipfile =  file, files = c(file.path(tempdir(), c(f1, f2, fileN, fileT, fileA, fbox1, fbox2, fmean1, fmean2, fileNG, fileTG)), fileClinic, fileDoc, fileDoc2), flags = "-j")
      
      # stop(file)
      if (file.exists(paste0( file, ".zip")))
        file.rename(paste0(file, ".zip"), file)
      
      ## cleaning the files after compressing and packing them	
      to_delete <-  c(f1, f2, fileN, fileT, fileA, fbox1, fbox2, fmean1, fmean2, fileNG, fileTG)
      for (fn in to_delete) 
        unlink(file.path(tempdir(),fn))
      
    },
    contentType = "application/zip"
  )

  
  cancel.onSessionEnded <- session$onSessionEnded(function() {
      dbDisconnect(con)
  })
  
})
