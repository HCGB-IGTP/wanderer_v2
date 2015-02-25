#!/usr/bin/env R

library(shiny)
library(RPostgreSQL)
library(Cairo)

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


sample_size <- read.table(file.path(SRC, "samplesN_filtered2.csv"), sep = ",", stringsAsFactors = FALSE, header = TRUE)

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

  api_arguments_allowances <- sort(c('Gene', 'start','end', 'TissueType', 'DataType', 'plotmean',
                                'geneLine', 'CpGi', 'nN', 'nT', 'region', 'goButton', 'Zoom'))
  
  
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
  
  
    
  ## #################################################
  ## #number of samples to plot
  ## output$nNmax <- renderUI({
  ##   if (!is.null(query$DataType) & !is.null(query$TissueType) & !is.null(toupper(query$Gene))) {
  ##     maxn <- max_sample(sample_size, query$DataType, query$TissueType)[1]
  ##     valor <- min(maxn, 30)
  ##     minn <- 1
  ##     conditionalPanel("query.plotmean == false", numericInput("nN", h5(paste0("Number of normal samples to plot (max = ", maxn, ")")), value = valor, min = minn, max = maxn))
  ##   }
  ## })
  
  ## output$nTmax <- renderUI({
  ##   if (!is.null(query$DataType) & !is.null(query$TissueType) & !is.null(toupper(query$Gene))) {
  ##     maxt <- max_sample(sample_size, query$DataType, query$TissueType)[2]
  ##     valor <- min(maxt, 30)
  ##     mint <- 1
  ##     conditionalPanel("query.plotmean == false", numericInput("nT", h5(paste0("Number of tumoral samples to plot (max = ", maxt, ")")), value = valor, min = mint, max = maxt))
  ##   } 
  ## })
  
  
  #################################################
  #print wanderer plot
  output$plot1 <- renderPlot({
      if (!all(api_arguments_allowances == sort(names(query))) | !(toupper(query$TissueType) %in% sample_size[,2]))
          stop(error)

      if(!is.null(query$TissueType) & !is.null(query$nN) & !is.null(query$nT) & !is.null(toupper(query$Gene)) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) { 
          if(query$region & ((query$end > query$Zoom[2]) | (query$end < query$Zoom[1]) | (query$start < query$Zoom[1]) | (query$start > query$Zoom[2])))
              stop(error)
          
          if(query$DataType == 'methylation'){
              if(dim(datamethfilt()[['probes2']])[1]<=0){
                  stop(error)
              }else{
                  wanderer_methylation(results_filt = datamethfilt(), geneName = toupper(query$Gene),
                                       geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                       CpGislands = query$CpGi, plotmean = query$plotmean,
                                       plotting = TRUE, geneLine = query$geneLine)
              }
          }
          else if(query$DataType == 'expression'){
              if(dim(dataexprfilt()[['exons2']])[1]<=0){
                  stop(error)
              }else{
                  wanderer_expression(results_filt = dataexprfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                      plotmean = query$plotmean, plotting = TRUE, geneLine = query$geneLine)
              }
          }

      }
  }, height = 1000, width = 1000)
  
  
  
  
  ##################################################
  #print summary plot
  output$plotStat <- renderPlot({

      if (!all(api_arguments_allowances == sort(names(query))) | !(toupper(query$TissueType) %in% sample_size[,2]))
          stop(error)

    if(!is.null(query$TissueType) & !is.null(query$nN) & !is.null(query$nT) & !is.null(toupper(query$Gene)) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      if(query$DataType == 'methylation'){
        if(dim(datamethfilt()[['probes2']])[1]>0){
          stat_analysis_meth(results_filt = datamethfilt(), geneName = toupper(query$Gene),
                             geneNamesType = geneFormat(), CpGislands = query$CpGi,
                             geneLine = query$geneLine, plotting = TRUE)
        }
      }
      if(query$DataType == 'expression'){
        if(dim(dataexprfilt()[['exons2']])[1]>0){
          stat_analysis_expr(results_filt = dataexprfilt(), geneName =toupper(query$Gene),
                             geneNamesType = geneFormat(),
                             geneLine = query$geneLine, plotting = TRUE)
        }
      }
    }
  }, height = 500, width = 1000)
  
  
  #################################################
  #dowload plot as png & pdf
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0("Wanderer_", toupper(query$Gene), '_', query$DataType, '_', query$TissueType, '_', Sys.Date(), '.png') },
    content = function(file) {
      CairoPNG(file, width = 1000, height = 1000)
      if(query$DataType == 'methylation'){
        regplot <- wanderer_methylation(results_filt = datamethfilt(), geneName = toupper(query$Gene),
                                        geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                        CpGislands = query$CpGi, plotmean = query$plotmean,
                                        plotting = TRUE, geneLine = query$geneLine)
      }
      else if(query$DataType == 'expression'){
        regplot <- wanderer_expression(results_filt = dataexprfilt(), geneName = toupper(query$Gene),
                                       geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                       plotmean = query$plotmean, plotting = TRUE, geneLine = query$geneLine)
      }

      dev.off()
    }
  )
  output$downloadPlotPDF <- downloadHandler(
    filename = function() { paste0("Wanderer_", toupper(query$Gene), '_', query$DataType, '_', query$TissueType, '_', Sys.Date(), '.pdf') },
    content = function(file) {
      pdf(file, width = 14, height = 14, useDingbats = FALSE)
      if(query$DataType == 'methylation'){
        regplot <- wanderer_methylation(results_filt = datamethfilt(), geneName = (toupper(query$Gene)),
                                        geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                        CpGislands = query$CpGi, plotmean = query$plotmean,
                                        plotting = TRUE, geneLine = query$geneLine)
      }
      else if(query$DataType == 'expression'){
        regplot <- wanderer_expression(results_filt = dataexprfilt(), geneName = (toupper(query$Gene)),
                                       geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                       plotmean = query$plotmean, plotting = TRUE, geneLine = query$geneLine)
      }

      dev.off()
    }
  )
  
  
  #################################################
  #dowload Normal data
  output$downloadNData <- downloadHandler(
    
    filename = function() { paste0("Wanderer_", isolate(toupper(query$Gene)), '_', query$DataType, '_', query$TissueType, '_Normal_', Sys.Date(), '.txt') },
    content = function(file) {
      if(query$DataType == 'methylation')
          write.table(datamethfilt()$ddN2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      else if(query$DataType == 'expression')
          write.table(dataexprfilt()$ddN2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
    }
  )
  
  #################################################
  #dowload Tumor data
  output$downloadTData <- downloadHandler(
    
    filename = function() { paste0("Wanderer_", (toupper(query$Gene)), '_', query$DataType, '_', query$TissueType, '_Tumor_', Sys.Date(), '.txt') },
    content = function(file) {
      if(query$DataType == 'methylation')
          write.table(datamethfilt()$ddT2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      else if(query$DataType == 'expression')
          write.table(dataexprfilt()$ddT2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)      
    }
  )
  
  #################################################
  #dowload probe annotation and statistical analysis
  output$downloadPData <- downloadHandler(
    
    filename = function() { paste0("Wanderer_", (toupper(query$Gene)), '_', query$DataType, '_', query$TissueType, '_annotations_and_statistical_analysis_', Sys.Date(), '.txt') },
    content = function(file) {
      if(query$DataType == 'methylation'){
        results <-stat_analysis_meth(results_filt = datamethfilt(), geneName = (toupper(query$Gene)),
                                     geneNamesType = geneFormat(), CpGislands = query$CpGi,
                                     geneLine = query$geneLine, plotting = FALSE) 
        write.table(results, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      else if(query$DataType == 'expression'){
        results <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(),
                                      geneLine = query$geneLine, plotting = FALSE)
        write.table(results, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
  
  #################################################
  #dowload plot as png & pdf
  output$downloadMeanPlot <- downloadHandler(
    filename = function() { paste0("Wanderer_", (toupper(query$Gene)), '_Mean_', query$DataType, '_', query$TissueType, '_', Sys.Date(), '.png') },
    content = function(file) {
      CairoPNG(file, width = 1000, height = 500)
      if(query$DataType == 'methylation'){
        regplot <- stat_analysis_meth(results_filt = datamethfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(), CpGislands = query$CpGi,
                                      geneLine = query$geneLine, plotting = TRUE)
      }
      else if(query$DataType == 'expression'){
        regplot <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(),
                                      geneLine = query$geneLine, plotting = TRUE)
      }
      dev.off()
    }
  )
  output$downloadMeanPlotPDF <- downloadHandler(
    filename = function() { paste0("Wanderer_", (toupper(query$Gene)), '_Mean_', query$DataType, '_', query$TissueType, '_', Sys.Date(), '.pdf') },
    content = function(file) {
      pdf(file, width = 14, height = 6, useDingbats = FALSE)
      if(query$DataType == 'methylation'){
        regplot <- stat_analysis_meth(results_filt = datamethfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(), CpGislands = query$CpGi,
                                      geneLine = query$geneLine, plotting = TRUE)
      }
      else if(query$DataType == 'expression'){
        regplot <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(),
                                      geneLine = query$geneLine, plotting = TRUE)
      }
      
      dev.off()
    }
  )

  cancel.onSessionEnded <- session$onSessionEnded(function() {
      dbDisconnect(con)
  })  
})
