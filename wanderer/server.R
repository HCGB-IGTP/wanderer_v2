#!/usr/bin/env R

library(shiny)
library(RPostgreSQL)
library(Cairo)

# the file containing the db parameters
#SRC <- '/imppc/labs/maplab/imallona/src/regional_profiler/wanderer'
## SRC <- '/imppc/labs/maplab/adiez/region_profile/Wanderer_190515/'
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

shinyServer(function(input, output, session){
  
  #database connection
  con <- db_connect(DB_CONF)
  
  
  #################################################
  #Gene Name
  geneNameSaved <- reactive({
    if(input$goButton == 0) {
      GeneNameSaved <- 'BRCA1'
    }else{
      GeneNameSaved <- isolate(toupper(input$Gene))
    }
  })
  
  
  #################################################
  ##detect gene format
  geneFormat <- reactive({
    if(input$goButton == 0) {
      GeneFormat <- "genename"
    } else{
      if(substr(geneNameSaved(), 1, 4) == "ENSG"){
        GeneFormat <- "emsemblgeneid"
      } else{
        GeneFormat <- "genename"
      }
    }
  })
  
  #################################################
  #Gene Size
  geneSize <- reactive({
    if(input$goButton == 0) {
      GeneSize(con = con, geneName = 'BRCA1', geneNamesType = 'genename')  
    }else{
      ## if (!is.null(input$Gene))
      GeneSize(con = con, geneName = geneNameSaved(), geneNamesType = geneFormat())  
    }
  })

  ## Sync gene symbol identifier for ENSG genenames; used when generating URL to external services
  gene_symbol <- reactive({
       if(input$goButton == 0) {
           curr_gene_symbol <- 'BRCA1'
       } else if (geneFormat() == 'genename') {
           curr_gene_symbol <- geneNameSaved()
       } else if (geneFormat() == 'emsemblgeneid') {
           curr_gene_symbol <- ensembl_to_gene_symbol(con, geneNameSaved())
       }
  })
  
  #################################################
  #zoom
  output$ZoomControl <- renderUI({
    
    if (geneFormat() == "genename") geneNamesType_label <- "Gene Name"
    if (geneFormat() == "emsemblgeneid") geneNamesType_label <- "Ensembl Gene ID"
    
    if(geneSize()[[1]]==0) stop(paste0("The gene ", geneNameSaved(), " is not in the ", geneNamesType_label, " annotation."))
    
    if(geneSize()[[1]]==1) stop(paste0("The gene ", geneNameSaved(), " appears more than once in the genome. Please introduce an Ensembl (ENSG) identifier instead."))
    
    if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
      
      sGene <- geneSize()[[1]]
      eGene <- geneSize()[[2]]
      sGene <- ((sGene%/%1000)-1)*1000
      eGene <- ((eGene%/%1000)+1)*1000
      
      
      if(input$DataType == 'methylation'){
        minGene <- sGene - 100000
        maxGene <- eGene + 100000
        tcks <- c(seq(minGene, sGene,(sGene-minGene)/3), seq(eGene, maxGene, (maxGene-eGene)/3))
        tcks <- (tcks%/%1000)*1000
        tcks <- format(tcks, big.mark = ',')        
      }
      if(input$DataType == 'expression'){
        minGene <- sGene
        maxGene <- eGene
        tcks <- seq(sGene, eGene, (eGene-sGene)/5)
        tcks <- (tcks%/%1000)*1000
        tcks <- format(tcks, big.mark = ',')        
      }
      
      sliderInput("Zoom", label = h5("Zoom"), min = minGene, max = maxGene, value = c(sGene, eGene), ticks = tcks, step = 1000, width = "800px")
    }
  })
  
  #################################################
  #region specification
  output$regionlimit <- renderUI({
    if(!is.null(input$Zoom)) conditionalPanel("input.region == true", helpText(paste0("Define a start and an end within the slider's values (min = ", input$Zoom[1],"; max = ", input$Zoom[2],")")))
  })
  
  output$start <- renderUI({
    if(!is.null(input$Zoom)) conditionalPanel("input.region == true", numericInput("start", "Start", value = input$Zoom[1], min = as.numeric(input$Zoom[1]), max = as.numeric(input$Zoom[2])))
  })
  
  output$end <- renderUI({
    if(!is.null(input$Zoom)) conditionalPanel("input.region == true", numericInput("end", "End", value = input$Zoom[2], min = as.numeric(input$Zoom[1]), max = as.numeric(input$Zoom[2])))
  })
  
  
  #################################################
  #Reading Methylation data
  datameth <- reactive({
    if(input$goButton == 0) {
      if(input$DataType == 'methylation' & !is.null(input$TissueType)){
        methylation_data(con = con, geneName = 'BRCA1', geneNamesType = 'genename', tissue = input$TissueType)  
      }
    }
    else if(input$goButton > 0 & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      methylation_data(con = con, geneName = geneNameSaved(), geneNamesType = geneFormat(), tissue = input$TissueType)
      
    }
  })
  
  #################################################
  #Filtering Methylation data
  datamethfilt <- reactive({
    if(input$DataType == 'methylation'){
      if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1 & !is.null(datameth()[['empty']])){
        if(!datameth()[['empty']])  data_meth_filtering(results = datameth(), sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end))
      }
    }
  })
  
  
  #################################################
  #Reading Expression data
  dataexpr <- reactive({
    if(input$goButton == 0) {      
      if(input$DataType == 'expression' & !is.null(input$TissueType)){
        expression_data(con = con, geneName = 'BRCA1', geneNamesType = 'genename', tissue = input$TissueType)  
      }
    } else if(input$goButton > 0 & geneSize()[[1]]!=0 & geneSize()[[1]]!=1 ) {
      expression_data(con = con, geneName = geneNameSaved(), geneNamesType = geneFormat(), tissue = input$TissueType)
    }
  })
  
  #################################################
  #Filtering Expression data
  dataexprfilt <- reactive({
    if(input$DataType == 'expression'){
      if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
        if(!dataexpr()[['empty']])  data_expr_filtering(results = dataexpr(), sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end))    
      }
    }
  })
  
  #################################################
  #RNAseq data at gene level
  dataRNAseqGene <- reactive({
    if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
      data_RNAseqGene(con = con, geneName = geneNameSaved(), geneNamesType = geneFormat(), tissue = input$TissueType)
    }
  })
  
  #################################################
  #RNAseq plot at gene level for expression data
  output$plotRNAseqGene <- renderPlot({
    if(!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(geneNameSaved()) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      if(input$DataType == 'expression'){
        if(!is.null(dataRNAseqGene()[['Normal']]) | !is.null(dataRNAseqGene()[['Tumor']])){
          plot_RNAseqGene(dd = dataRNAseqGene(), geneName = geneNameSaved(), tissue_label = dataexprfilt()[['tissue_label']]) 
        } 
      }
    }
  }, height = 500, width = 1000)
  
  
  #################################################
  #Probe selection
  output$ProbeSelected <- renderUI({
    if (!is.null(input$DataType) & !is.null(input$TissueType) & !is.null(geneNameSaved())) {
      if(input$DataType == 'methylation'){
        probes2Select <- datamethfilt()[['probes2']]
        probes2Select <- probes2Select$probe
        selectInput("ProbeSelection", label = h5("Select a probe:"), choices = probes2Select, selected = probes2Select[1])
      }
    }
  })
  
  #################################################
  #correl plot of meth and RNAseq at gene level for methylation data
  output$correlmethexpress <- renderPlot({
    if(!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(geneNameSaved()) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      if(input$DataType == 'methylation'){
        #if(!is.null(input$ProbeSelection) & !is.null(dataRNAseqGene()) & !is.null(datamethfilt()) & dim(datamethfilt()[['ddN2']])[1]>0){
        if((is.null(dataRNAseqGene()[['Normal']]) | is.null(datamethfilt()[['ddN2']])) & (is.null(dataRNAseqGene()[['Tumor']]) | is.null(datamethfilt()[['ddT2']]))){
          stop("No common RNAseq data available")
        }
        
        if((!is.null(dataRNAseqGene()[['Normal']]) & !is.null(datamethfilt()[['ddN2']])) | (!is.null(dataRNAseqGene()[['Tumor']]) & !is.null(datamethfilt()[['ddT2']]))){
          correl_meth_express(geneName = geneNameSaved(), probeID = input$ProbeSelection, ddmeth = datamethfilt(), ddGene = dataRNAseqGene(),  tissue_label = datamethfilt()[['tissue_label']], regressLine = input$regressionLine, correlMethod = input$correlationMethod, plotting = TRUE, datareturn = FALSE) 
        } 
      }
    }
  }, height = 500, width = 1000)
  
  #################################################
  #index html correl parameters to print
  output$correlParam <- reactive({
    correlparam <- 0
    if(input$DataType == 'methylation'){
      if((is.null(dataRNAseqGene()[['Normal']]) | is.null(datamethfilt()[['ddN2']])) & (is.null(dataRNAseqGene()[['Tumor']]) | is.null(datamethfilt()[['ddT2']]))){
        correlparam <- 0
      }
      if((!is.null(dataRNAseqGene()[['Normal']]) & !is.null(datamethfilt()[['ddN2']])) | (!is.null(dataRNAseqGene()[['Tumor']]) & !is.null(datamethfilt()[['ddT2']]))){
        correlparam <- 1
      }
    }
    return(correlparam)
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
  #Tissue selection
  output$Tissues <- renderUI({
    tissues <- paste0("'",sample_size[,1], " (", sample_size[,2], ")", "'='", sample_size[,3], "'")
    tissues <- paste0(tissues,collapse=",")
    tissues <- paste0("c(",tissues,")")
    tissues <- eval(parse(text=tissues))
    selectInput("TissueType", label = h5("Project:"), choices = tissues, selected = "brca")
  })
  
  #################################################
  ##number of samples to plot
  output$nNmax <- renderUI({
    if (!is.null(input$DataType) & !is.null(input$TissueType) & !is.null(geneNameSaved())) {
      if(is.null(datamethfilt()$ddN2) | is.null(dataexprfilt()$ddN2)) minn <- 0
      if(!is.null(datamethfilt()$ddN2) & !is.null(dataexprfilt()$ddN2)) minn <- 1
      
      maxn <- max_sample(sample_size, input$DataType, input$TissueType)[1]
      valor <- min(maxn, 30)
      conditionalPanel("input.plotmean == false", numericInput("nN", h5(paste0("Number of normal samples to plot (max = ", maxn, ")")), value = valor, min = minn, max = maxn))
    }
  })
  
  output$nTmax <- renderUI({
    if (!is.null(input$DataType) & !is.null(input$TissueType) & !is.null(geneNameSaved())) {
      maxt <- max_sample(sample_size, input$DataType, input$TissueType)[2]
      valor <- min(maxt, 30)
      mint <- 1
      conditionalPanel("input.plotmean == false", numericInput("nT", h5(paste0("Number of tumoral samples to plot (max = ", maxt, ")")), value = valor, min = mint, max = maxt))
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
      pdf(file.path(tempdir(), f2), width = 14, height = 14, useDingbats = FALSE)
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
          pdf(file.path(tempdir(),fbox2), width = 12, height = 6, useDingbats = FALSE)
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
          pdf(file.path(tempdir(),fbox2), width = 12, height = 6, useDingbats = FALSE)
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
            pdf(file.path(tempdir(), fmean2), width = 14, height = 6, useDingbats = FALSE)
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
            pdf(file.path(tempdir(), fmean2), width = 14, height = 6, useDingbats = FALSE)
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
  
  
  
  output$permalink_modal <- renderText({
      ## location <- 'http://gattaca.imppc.org:3939/betawanderer_api'
      ## location <- 'http://192.168.3.67:3939/new_wanderer_api/'
      ## location <- 'http://gattaca.imppc.org:3838/wanderer_api'
      location <- 'http://www.maplab.cat/wanderer_api'
      ## if (input$goButton > 0) {
      
      if (input$region & !is.null(input$Gene) & !is.null(input$start) & !is.null(input$end) & !is.null(input$TissueType) & !is.null(input$DataType)
          & !is.null(input$plotmean) & !is.null(input$geneLine) & !is.null(input$CpGi) & !is.null(input$nN) & !is.null(input$nT)) {  
          url <- sprintf('%s?Gene=%s&start=%s&end=%s&TissueType=%s&DataType=%s&plotmean=%s&geneLine=%s&CpGi=%s&nN=%s&nT=%s&distribute_uniformly=%s&pvalThres=%s',
                         location,
                         input$Gene,
                         input$start,
                         input$end,
                         input$TissueType,
                         input$DataType,
                         input$plotmean,
                         input$geneLine,
                         input$CpGi,
                         input$nN,
                         input$nT,
                         input$distribute_uniformly,
                         input$pvalThres)
      }
      else if (!input$region & !is.null(input$Gene) & !is.null(input$Zoom[0]) & !is.null(input$Zoom[1]) & !is.null(input$TissueType) & !is.null(input$DataType)
               & !is.null(input$plotmean) & !is.null(input$geneLine) & !is.null(input$CpGi) & !is.null(input$nN) & !is.null(input$nT)) {          
          
          url = sprintf('%s?Gene=%s&start=%s&end=%s&TissueType=%s&DataType=%s&plotmean=%s&geneLine=%s&CpGi=%s&nN=%s&nT=%s&distribute_uniformly=%s&pvalThres=%s',
              location,
              input$Gene,
              input$Zoom[1],
              input$Zoom[2],
              input$TissueType,
              input$DataType,
              input$plotmean,
              input$geneLine,
              input$CpGi,
              input$nN,
              input$nT,
              input$distribute_uniformly,
              input$pvalThres)
      } else {
          url <- sprintf('%s/?Gene=BRCA1&start=41195000&end=41278000&TissueType=brca&DataType=methylation&plotmean=FALSE&geneLine=TRUE&CpGi=TRUE&nN=30&nT=30&distribute_uniformly=FALSE&pvalThres=0.05', location)
          
      }
      
      generate_modal(url)
  })

    ## external link to genome browser
  output$genome_browser <- renderText({
    generate_genome_browser_link(chromosome = geneSize()[[3]],
                                 start = input$Zoom[1],
                                 end = input$Zoom[2])    
  })

  ## external link to cbioportal
  output$cbioportal <- renderText({
      if (!is.null(input$Gene) & !is.null(input$TissueType) & !is.null(input$DataType)) {
          if (tolower(input$TissueType) %in% names(cbioportal_datasets())) {
              ## the equivalence to the ENSG identifier, if any
              gene_symbol <- gene_symbol()
              
              if (is.character(gene_symbol) & nchar(gene_symbol) > 0) {
                  generate_cbioportal_link(dataset = input$TissueType,
                                           gene = gene_symbol)
              }
          }
      }
  })

  ## external link to regulome explorer
  output$regulome_explorer <- renderText({
      if (!is.null(input$Gene) & !is.null(input$TissueType) & !is.null(input$DataType)) {
          if (tolower(input$TissueType) %in% names(regulome_explorer_datasets())) {
               ## the equivalence to the ENSG identifier, if any
              gene_symbol <- gene_symbol()
              
              if (is.character(gene_symbol) & nchar(gene_symbol) > 0) {
              
                  generate_regulome_explorer_link(dataset = input$TissueType,
                                                  gene = gene_symbol)
              }
          }
      }          
  })
  
  
  cancel.onSessionEnded <- session$onSessionEnded(function() {
    dbDisconnect(con)
      if (length(dev.list()) > 0){          
          for (i in 1:length(dev.list())) {
              dev.off()
          }
      }
  })
  
})
