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
  

  ## ## ###############################################
  ## #Gene Name
  ## geneNameSaved <- function(){
  ##     return(toupper(query$Gene))
  ## }
  
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
  
  ## main profile plot
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
                                       plotting = TRUE, geneLine = query$geneLine,
                                       proportional = !(query$distribute_uniformly))
              }
          }
          else if(query$DataType == 'expression'){
              if(dim(dataexprfilt()[['exons2']])[1]<=0){
                  stop(error)
              }else{
                  wanderer_expression(results_filt = dataexprfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                      plotmean = query$plotmean, plotting = TRUE, geneLine = query$geneLine,
                                      proportional = !(query$distribute_uniformly))
              }
          }

      }
  }, height = 1000, width = 1000)

  
  ## profile NT comparison
  output$plotStat <- renderPlot({

      if (!all(api_arguments_allowances == sort(names(query))) | !(toupper(query$TissueType) %in% sample_size[,2]))
          stop(error)

      if(!is.null(query$TissueType) & !is.null(query$nN) & !is.null(query$nT) & !is.null(toupper(query$Gene)) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
          if(query$DataType == 'methylation'){
              if(dim(datamethfilt()[['probes2']])[1]>0){
                  stat_analysis_meth(results_filt = datamethfilt(), geneName = toupper(query$Gene),
                                     geneNamesType = geneFormat(), CpGislands = query$CpGi,
                                     geneLine = query$geneLine, plotting = TRUE,
                                     proportional = !(query$distribute_uniformly),
                                     pvalThres = query$pvalThres)
              }
          }
          if(query$DataType == 'expression'){
              if(dim(dataexprfilt()[['exons2']])[1]>0){
                  stat_analysis_expr(results_filt = dataexprfilt(), geneName =toupper(query$Gene),
                                     geneNamesType = geneFormat(),
                                     geneLine = query$geneLine, plotting = TRUE,
                                     proportional = !(query$distribute_uniformly),
                                     pvalThres = query$pvalThres)
              }
          }
      }
  }, height = 500, width = 1000)

    
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
  
  
  cancel.onSessionEnded <- session$onSessionEnded(function() {
      dbDisconnect(con)
  })
  
})
