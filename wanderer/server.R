#!/usr/bin/env R

library(shiny)
library(RPostgreSQL)
library(Cairo)

# the file containing the db parameters
#SRC <- '/imppc/labs/maplab/imallona/src/regional_profiler/wanderer'
#SRC <- '/imppc/labs/maplab/adiez/region_profile/Wanderer_060215/'
SRC <- '.'

DB_CONF <- file.path(SRC, 'db.txt')

source(file.path(SRC, 'GeneSize.R'))
source(file.path(SRC, 'expression_data.R'))
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


sample_size <- read.table(file.path(SRC, "samplesN_filtered2.csv"), sep = ",", stringsAsFactors = FALSE, header = TRUE)

shinyServer(function(input, output, session){
  
  #database connection
  con <- db_connect(DB_CONF)
  ## print(dbGetInfo(con))
  #################################################
  ##detect gene format
  geneFormat <- reactive({
    if(input$goButton == 0) {
      GeneFormat <- "genename"
    } else{
      if(substr(isolate(toupper(input$Gene)), 1, 4) == "ENSG"){
        GeneFormat <- "emsemblgeneid"
      } else{
        GeneFormat <- "genename"
      }
    }
  })
  
  #################################################
  #Filtering Methylation data
  geneSize <- reactive({
    if(input$goButton == 0) {
      GeneSize(con = con, geneName = 'BRCA1', geneNamesType = 'genename')  
    }else{
      ## if (!is.null(input$Gene))
      GeneSize(con = con, geneName = isolate(toupper(input$Gene)), geneNamesType = geneFormat())  
    }
  })
  
  #################################################
  #zoom
  output$ZoomControl <- renderUI({
    
    if (geneFormat() == "genename") geneNamesType_label <- "Gene Name"
    if (geneFormat() == "emsemblgeneid") geneNamesType_label <- "Ensembl Gene ID"
    
    if(geneSize()[[1]]==0) stop(paste0("The gene ", isolate(toupper(input$Gene)), " is not in the ", geneNamesType_label, " annotation."))
    
    if(geneSize()[[1]]==1) stop(paste0("The gene ", isolate(toupper(input$Gene)), " appears more than once in the genome. Please introduce an Ensembl (ENSG) identifier instead."))
    
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
      methylation_data(con = con, geneName = isolate(toupper(input$Gene)), geneNamesType = geneFormat(), tissue = input$TissueType)
      
    }
  })
  
  #################################################
  #Filtering Methylation data
  datamethfilt <- reactive({
    if(input$DataType == 'methylation'){
      if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
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
      expression_data(con = con, geneName = isolate(toupper(input$Gene)), geneNamesType = geneFormat(), tissue = input$TissueType)
      ## }
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
    if (!is.null(input$DataType) & !is.null(input$TissueType) & !is.null(isolate(toupper(input$Gene)))) {
      if(is.null(datamethfilt()$ddN2) | is.null(dataexprfilt()$ddN2)) minn <- 0
      if(!is.null(datamethfilt()$ddN2) & !is.null(dataexprfilt()$ddN2)) minn <- 1
      
      maxn <- max_sample(sample_size, input$DataType, input$TissueType)[1]
      valor <- min(maxn, 30)
      conditionalPanel("input.plotmean == false", numericInput("nN", h5(paste0("Number of normal samples to plot (max = ", maxn, ")")), value = valor, min = minn, max = maxn))
    }
  })
  
  output$nTmax <- renderUI({
    if (!is.null(input$DataType) & !is.null(input$TissueType) & !is.null(isolate(toupper(input$Gene)))) {
      maxt <- max_sample(sample_size, input$DataType, input$TissueType)[2]
      valor <- min(maxt, 30)
      mint <- 1
      conditionalPanel("input.plotmean == false", numericInput("nT", h5(paste0("Number of tumoral samples to plot (max = ", maxt, ")")), value = valor, min = mint, max = maxt))
    }
    
  })
  
  
  #################################################
  #print wanderer plot
  output$plot1 <- renderPlot({
    if(!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(isolate(toupper(input$Gene))) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      if(input$region & ((input$end > input$Zoom[2]) | (input$end < input$Zoom[1]) | (input$start < input$Zoom[1]) | (input$start > input$Zoom[2]))) stop(print(paste0("The region must be between ", input$Zoom[1], " and ", input$Zoom[2])))
      
      if(input$DataType == 'methylation'){
        if(dim(datamethfilt()[['probes2']])[1]>0){
          if(is.null(datamethfilt()$ddN2) & is.null(datamethfilt()$ddT2)){
            stop("There are no samples in this tissue type")
          } else{
            wanderer_methylation(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                 geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                 CpGislands = input$CpGi, plotmean = input$plotmean,
                                 plotting = TRUE, geneLine = input$geneLine)
          }
        }
      }
      if(input$DataType == 'expression'){
        if(dim(dataexprfilt()[['exons2']])[1]>0){
          if(is.null(dataexprfilt()$ddN2) & is.null(dataexprfilt()$ddT2)){
            stop("There are no samples in this tissue type")
          } else{
            wanderer_expression(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                                geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
          }
        }
      }
    }
  }, height = 1000, width = 1000)
  
  
  
  ##################################################
  #print summary plot
  output$plotStat <- renderPlot({
    if(!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(isolate(toupper(input$Gene))) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      
      if(input$DataType == 'methylation'){
        if(dim(datamethfilt()[['probes2']])[1]>0){
          if(is.null(datamethfilt()$ddN2) | is.null(datamethfilt()$ddT2)){
            stop("There are not enough samples to perform the statistical analysis")
          } else{
            
            if(dim(datamethfilt()$ddN2)[2]==2 | dim(datamethfilt()$ddT2)[2]==2){
              stop("There are not enough samples to perform the statistical analysis")
            } else{
              stat_analysis_meth(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                 geneNamesType = geneFormat(), CpGislands = input$CpGi,
                                 geneLine = input$geneLine, plotting = TRUE)
            }
          }
        }
      }
      else if(input$DataType == 'expression'){
        if(dim(dataexprfilt()[['exons2']])[1]>0){
          if(is.null(dataexprfilt()$ddN2) | is.null(dataexprfilt()$ddT2)){
            stop("There are not enough samples to perform the statistical analysis")
          } else{
            
            if(dim(dataexprfilt()$ddN2)[2]==2 | dim(dataexprfilt()$ddT2)[2]==2){
              stop("There are not enough samples to perform the statistical analysis")
            } else{
              stat_analysis_expr(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                                 geneNamesType = geneFormat(),
                                 geneLine = input$geneLine, plotting = TRUE)
            }
          }
        }
      }
    }
  }, height = 500, width = 1000)
  
  
  #################################################
  #dowload plot as png & pdf
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0("Wanderer_", isolate(toupper(input$Gene)), '_', input$DataType, '_', input$TissueType, '_', Sys.Date(), '.png') },
    content = function(file) {
      CairoPNG(file, width = 1000, height = 1000)
      if(input$DataType == 'methylation'){
        regplot <- wanderer_methylation(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                        geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                        CpGislands = input$CpGi, plotmean = input$plotmean,
                                        plotting = TRUE, geneLine = input$geneLine)
      }
      else if(input$DataType == 'expression'){
        regplot <- wanderer_expression(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                                       geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                       plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
      }
      print(regplot)
      dev.off()
    }
  )
  output$downloadPlotPDF <- downloadHandler(
    filename = function() { paste0("Wanderer_", isolate(toupper(input$Gene)), '_', input$DataType, '_', input$TissueType, '_', Sys.Date(), '.pdf') },
    content = function(file) {
      pdf(file, width = 14, height = 14, useDingbats = FALSE)
      if(input$DataType == 'methylation'){
        regplot <- wanderer_methylation(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                        geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                        CpGislands = input$CpGi, plotmean = input$plotmean,
                                        plotting = TRUE, geneLine = input$geneLine)
      }
      else if(input$DataType == 'expression'){
        regplot <- wanderer_expression(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                                       geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                       plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
      }
      print(regplot)
      dev.off()
    }
  )
  
  
  #################################################
  #dowload Normal data
  output$downloadNData <- downloadHandler(
    
    filename = function() { paste0("Wanderer_", isolate(toupper(input$Gene)), '_', input$DataType, '_', input$TissueType, '_Normal_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation')  write.table(datamethfilt()$ddN2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      if(input$DataType == 'expression')  write.table(dataexprfilt()$ddN2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
    }
  )
  
  #################################################
  #dowload Tumor data
  output$downloadTData <- downloadHandler(
    
    filename = function() { paste0("Wanderer_", isolate(toupper(input$Gene)), '_', input$DataType, '_', input$TissueType, '_Tumor_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation')  write.table(datamethfilt()$ddT2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      if(input$DataType == 'expression')  write.table(dataexprfilt()$ddT2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)      
    }
  )
  
  #################################################
  #dowload probe annotation and statistical analysis
  output$downloadPData <- downloadHandler(
    
    filename = function() { paste0("Wanderer_", isolate(toupper(input$Gene)), '_', input$DataType, '_', input$TissueType, '_annotations_and_statistical_analysis_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation'){
        results <-stat_analysis_meth(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                     geneNamesType = geneFormat(), CpGislands = input$CpGi,
                                     geneLine = input$geneLine, plotting = FALSE) 
        write.table(results, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      if(input$DataType == 'expression'){
        results <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                                      geneNamesType = geneFormat(),
                                      geneLine = input$geneLine, plotting = FALSE)
        write.table(results, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
  
  #################################################
  #dowload plot as png & pdf
  output$downloadMeanPlot <- downloadHandler(
    filename = function() { paste0("Wanderer_", isolate(toupper(input$Gene)), '_Mean_', input$DataType, '_', input$TissueType, '_', Sys.Date(), '.png') },
    content = function(file) {
      CairoPNG(file, width = 1000, height = 500)
      if(input$DataType == 'methylation'){
        regplot <- stat_analysis_meth(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                      geneNamesType = geneFormat(), CpGislands = input$CpGi,
                                      geneLine = input$geneLine, plotting = TRUE)
      }
      else if(input$DataType == 'expression'){
        regplot <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                                      geneNamesType = geneFormat(),
                                      geneLine = input$geneLine, plotting = TRUE)
      }
      ## print(regplot)
      dev.off()
    }
  )
  output$downloadMeanPlotPDF <- downloadHandler(
    filename = function() { paste0("Wanderer_", isolate(toupper(input$Gene)), '_Mean_', input$DataType, '_', input$TissueType, '_', Sys.Date(), '.pdf') },
    content = function(file) {
      pdf(file, width = 14, height = 6, useDingbats = FALSE)
      if(input$DataType == 'methylation'){
        regplot <- stat_analysis_meth(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                      geneNamesType = geneFormat(), CpGislands = input$CpGi,
                                      geneLine = input$geneLine, plotting = TRUE)
      }
      else if(input$DataType == 'expression'){
        regplot <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                                      geneNamesType = geneFormat(),
                                      geneLine = input$geneLine, plotting = TRUE)
      }
      print(regplot)
      dev.off()
    }
  )
  
  output$permalink_modal <- renderText({
    ## location <- 'http://gattaca.imppc.org:3939/betawanderer_api'
    location <- 'http://192.168.3.67:3939/new_wanderer_api/'
    
    ## if (input$goButton > 0) {
    
    if (input$region & !is.null(input$Gene) & !is.null(input$start) & !is.null(input$end) & !is.null(input$TissueType) & !is.null(input$DataType)
        & !is.null(input$plotmean) & !is.null(input$geneLine) & !is.null(input$CpGi) & !is.null(input$nN) & !is.null(input$nT)) {  
      url <- sprintf('%s?Gene=%s&start=%s&end=%s&TissueType=%s&DataType=%s&plotmean=%s&geneLine=%s&CpGi=%s&nN=%s&nT=%s',
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
                     input$nT)
    }
    else if (!input$region & !is.null(input$Gene) & !is.null(input$Zoom[0]) & !is.null(input$Zoom[1]) & !is.null(input$TissueType) & !is.null(input$DataType)
             & !is.null(input$plotmean) & !is.null(input$geneLine) & !is.null(input$CpGi) & !is.null(input$nN) & !is.null(input$nT)) {          
      
      url = sprintf('%s?Gene=%s&start=%s&end=%s&TissueType=%s&DataType=%s&plotmean=%s&geneLine=%s&CpGi=%s&nN=%s&nT=%s',
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
                    input$nT)
    } else {
      url <- sprintf('%s/?Gene=BRCA1&start=41195000&end=41278000&TissueType=brca&DataType=methylation&plotmean=FALSE&geneLine=TRUE&CpGi=TRUE&nN=30&nT=30', location)
      
    }
    
    generate_modal(url)
  })
  
  output$genome_browser <- renderText({
    generate_genome_browser_link(chromosome = geneSize()[[3]],
                                 start = input$Zoom[1],
                                 end = input$Zoom[2])
    
  })

  
  output$proportionality <- renderText({
     ## generate_modal('foobar')
      pop_modal_plot('foobar')
    
  })

  
  cancel.onSessionEnded <- session$onSessionEnded(function() {
    dbDisconnect(con)
  })

  output$plot_test <- renderPlot({
      hist(runif(20))
  })
  
  #on.exit(dbDisconnect(con), add = TRUE)
  
})
