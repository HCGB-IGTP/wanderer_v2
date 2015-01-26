
library(shiny)
library(RPostgreSQL)

# the file containing the db parameters
SRC <- '/imppc/labs/maplab/imallona/src/regional_profiler/wanderer'
#SRC <- '/data/shiny/apps/wanderer'

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


sample_size <- read.table(file.path(SRC, "samplesN_filtered.csv"), sep = ",", stringsAsFactors = FALSE, header = TRUE)


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

  query<-isolate( parseQueryString(session$clientData$url_search))
  
  ## query <- parseQueryString(session, reactive(session$clientData$url_search))

  print(paste(names(query), query, sep = "=", collapse=", "))

  ## 127.0.0.1:6681?Gene=BRCA1&sGene=41195000&eGene=41278000&start=41195000&end=41278000&TissueType=brca&goButton=1&DataType=methylation&goButton=1&plotmean=FALSE&geneLine=TRUE&CpGi=FALSE&nN=30&nT=30&region=TRUE

  ## http://127.0.0.1:3013/?Gene=BRCA1&start=41198000&end=41278000&TissueType=brca&goButton=1&DataType=methylation&goButton=1&plotmean=FALSE&geneLine=TRUE&CpGi=FALSE&nN=30&nT=10&region=TRUE

  query$Zoom <- as.numeric(c(query$start, query$end))
  query$start <- as.numeric(query$start)
  query$end <- as.numeric(query$end)
  ## query$eGene <- as.numeric(query$eGene)
  ## query$sGene <- as.numeric(query$sGene)
  query$nN <- as.numeric(query$nN)
  query$nT <- as.numeric(query$nT)
  query$plotmean <- as.logical(query$plotmean)
  query$geneLine <- as.logical(query$geneLine)
  query$CpGi <- as.logical(query$CpGi)
  query$region <- as.logical(query$region)
  
  print(query)
  
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
    if(query$goButton == 0) {
      GeneSize(con = con, geneName = 'TP53', geneNamesType = 'genename')  
    }else{
      GeneSize(con = con, geneName =toupper(query$Gene), geneNamesType = geneFormat())  
    }
  })
  
  #################################################
  #zoom
  output$ZoomControl <- renderUI({
    
    if (geneFormat() == "genename") geneNamesType_label <- "Gene Name"
    if (geneFormat() == "emsemblgeneid") geneNamesType_label <- "Ensembl Gene ID"
    
    if(geneSize()[[1]]==0) stop(paste0("The gene ", toupper(query$Gene), " is not in the ", geneNamesType_label, " annotation."))
    
    if(geneSize()[[1]]==1) stop(paste0("The gene ", toupper(query$Gene), " appears more than once in the genome. Please introduce an Ensembl (ENSG) identifier instead."))
    
    if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
      
      sGene <- geneSize()[[1]]
      eGene <- geneSize()[[2]]
      sGene <- ((sGene%/%1000)-1)*1000
      eGene <- ((eGene%/%1000)+1)*1000
      
      
      if(query$DataType == 'methylation'){
        minGene <- sGene - 100000
        maxGene <- eGene + 100000
        tcks <- c(seq(minGene, sGene,(sGene-minGene)/3), seq(eGene, maxGene, (maxGene-eGene)/3))
        tcks <- (tcks%/%1000)*1000
        tcks <- format(tcks, big.mark = ',')        
      }
      if(query$DataType == 'expression'){
        minGene <- sGene
        maxGene <- eGene
        tcks <- seq(sGene, eGene, (eGene-sGene)/5)
        tcks <- (tcks%/%1000)*1000
        tcks <- format(tcks, big.mark = ',')        
      }
      
      sliderInput("Zoom", label = h5("Zoom"), min = minGene, max = maxGene, value = c(sGene, eGene), ticks = tcks, step = 1000, width = "800px")
    }
  })
  
  
  ## region specification
  output$regionlimit <- renderUI({
    if(!is.null(query$Zoom)) conditionalPanel("query.region == true", helpText(paste0("Define a start and an end within the slider's values (min = ", query$Zoom[1],"; max = ", query$Zoom[2],")")))
  })
  
  output$start <- renderUI({
    if(!is.null(query$Zoom)) conditionalPanel("query.region == true", numericInput("start", "Start", value = query$Zoom[1], min = as.numeric(query$Zoom[1]), max = as.numeric(query$Zoom[2])))
  })
  
  output$end <- renderUI({
    if(!is.null(query$Zoom)) conditionalPanel("query.region == true", numericInput("end", "End", value = query$Zoom[2], min = as.numeric(query$Zoom[1]), max = as.numeric(query$Zoom[2])))
  })
  
  
  #################################################
  #Reading Methylation data
  datameth <- reactive({
    if(query$DataType == 'methylation'){
      if(query$goButton == 0) {
        methylation_data(con = con, geneName = 'TP53', geneNamesType = 'genename', tissue = 'brca')  
      }else{
        if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1)  methylation_data(con = con, geneName = toupper(query$Gene), geneNamesType = geneFormat(), tissue = query$TissueType)
      }
    }
  })
  
  #################################################
  #Filtering Methylation data
  datamethfilt <- reactive({
    if(query$DataType == 'methylation'){
      if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
        if(!datameth()[['empty']])  data_meth_filtering(results = datameth(), sampleSize = sample_size, tissue = query$TissueType, zoom = c(query$start, query$end))
      }
    }
  })
  
  
  #################################################
  #Reading Expression data
  dataexpr <- reactive({
    if(query$DataType == 'expression'){
      if(query$goButton == 0) {
        expression_data(con = con, geneName = 'TP53', geneNamesType = 'genename', tissue = 'brca')  
      }else{
        if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1)  expression_data(con = con, geneName = toupper(query$Gene), geneNamesType = geneFormat(), tissue = query$TissueType)
      }
    }
  })
  
  #################################################
  #Filtering Expression data
  dataexprfilt <- reactive({
    if(query$DataType == 'expression'){
      if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
        if(!dataexpr()[['empty']])  data_expr_filtering(results = dataexpr(), sampleSize = sample_size, tissue = query$TissueType, zoom = c(query$start, query$end))    
      }
    }
  })
  
  #################################################
  #print the number of probes or exons
  output$numberpoints <- renderText({
    if(geneSize()[[1]]!=0 & geneSize()[[1]]!=1){
      if(query$DataType == 'methylation'){
        if(!is.null(datamethfilt())){
          if(dim(datamethfilt()[['probes2']])[1]==0) printa <- "There are not probes in this region"
          if(dim(datamethfilt()[['probes2']])[1]>0)  printa <- paste0("There are ", dim(datamethfilt()[['probes2']])[1] ," probes in the selected region")
        }
      }
      if(query$DataType == 'expression'){
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
  #number of samples to plot
  output$nNmax <- renderUI({
    if (!is.null(query$DataType) & !is.null(query$TissueType) & !is.null(toupper(query$Gene))) {
      maxn <- max_sample(sample_size, query$DataType, query$TissueType)[1]
      valor <- min(maxn, 30)
      minn <- 1
      conditionalPanel("query.plotmean == false", numericInput("nN", h5(paste0("Number of normal samples to plot (max = ", maxn, ")")), value = valor, min = minn, max = maxn))
    }
  })
  
  output$nTmax <- renderUI({
    if (!is.null(query$DataType) & !is.null(query$TissueType) & !is.null(toupper(query$Gene))) {
      maxt <- max_sample(sample_size, query$DataType, query$TissueType)[2]
      valor <- min(maxt, 30)
      mint <- 1
      conditionalPanel("query.plotmean == false", numericInput("nT", h5(paste0("Number of tumoral samples to plot (max = ", maxt, ")")), value = valor, min = mint, max = maxt))
    } 
  })
  
  
  #################################################
  #print wanderer plot
  output$plot1 <- renderPlot({
    if(!is.null(query$TissueType) & !is.null(query$nN) & !is.null(query$nT) & !is.null(toupper(query$Gene)) & geneSize()[[1]]!=0 & geneSize()[[1]]!=1) {
      if(query$region & ((query$end > query$Zoom[2]) | (query$end < query$Zoom[1]) | (query$start < query$Zoom[1]) | (query$start > query$Zoom[2]))) stop(print(paste0("The region must be between ", query$Zoom[1], " and ", query$Zoom[2])))
      
      if(query$DataType == 'methylation'){
        if(dim(datamethfilt()[['probes2']])[1]>0){
          #           stop(print(paste0("There are not probes in this region")))
          #         }else{
          wanderer_methylation(results_filt = datamethfilt(), geneName = toupper(query$Gene),
                               geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                               CpGislands = query$CpGi, plotmean = query$plotmean,
                               plotting = TRUE, geneLine = query$geneLine)
        }
      }
      if(query$DataType == 'expression'){
        if(dim(dataexprfilt()[['exons2']])[1]>0){
          #           stop(print(paste0("There are not exons in this region")))
          #         }else{
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
      png(file, width = 1000, height = 1000)
      if(query$DataType == 'methylation'){
        regplot <- wanderer_methylation(results_filt = datamethfilt(), geneName = toupper(query$Gene),
                                        geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                        CpGislands = query$CpGi, plotmean = query$plotmean,
                                        plotting = TRUE, geneLine = query$geneLine)
      }
      if(query$DataType == 'expression'){
        regplot <- wanderer_expression(results_filt = dataexprfilt(), geneName = toupper(query$Gene),
                                       geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                       plotmean = query$plotmean, plotting = TRUE, geneLine = query$geneLine)
      }
      print(regplot)
      dev.off()
    }
  )
  output$downloadPlotPDF <- downloadHandler(
    filename = function() { paste0("Wanderer_", toupper(query$Gene), '_', query$DataType, '_', query$TissueType, '_', Sys.Date(), '.pdf') },
    content = function(file) {
      pdf(file, width = 10, height = 13)
      if(query$DataType == 'methylation'){
        regplot <- wanderer_methylation(results_filt = datamethfilt(), geneName = (toupper(query$Gene)),
                                        geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                        CpGislands = query$CpGi, plotmean = query$plotmean,
                                        plotting = TRUE, geneLine = query$geneLine)
      }
      if(query$DataType == 'expression'){
        regplot <- wanderer_expression(results_filt = dataexprfilt(), geneName = (toupper(query$Gene)),
                                       geneNamesType = geneFormat(), npointsN = query$nN, npointsT = query$nT,
                                       plotmean = query$plotmean, plotting = TRUE, geneLine = query$geneLine)
      }
      print(regplot)
      dev.off()
    }
  )
  
  
  #################################################
  #dowload Normal data
  output$downloadNData <- downloadHandler(
    
    filename = function() { paste0("Wanderer_", isolate(toupper(query$Gene)), '_', query$DataType, '_', query$TissueType, '_Normal_', Sys.Date(), '.txt') },
    content = function(file) {
      if(query$DataType == 'methylation')  write.table(datamethfilt()$ddN2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      if(query$DataType == 'expression')  write.table(dataexprfilt()$ddN2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
    }
  )
  
  #################################################
  #dowload Tumor data
  output$downloadTData <- downloadHandler(
    
    filename = function() { paste0("Wanderer_", (toupper(query$Gene)), '_', query$DataType, '_', query$TissueType, '_Tumor_', Sys.Date(), '.txt') },
    content = function(file) {
      if(query$DataType == 'methylation')  write.table(datamethfilt()$ddT2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      if(query$DataType == 'expression')  write.table(dataexprfilt()$ddT2, file = file, sep = "\t", row.names = FALSE, quote = FALSE)      
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
      if(query$DataType == 'expression'){
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
      png(file, width = 1000, height = 500)
      if(query$DataType == 'methylation'){
        regplot <- stat_analysis_meth(results_filt = datamethfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(), CpGislands = query$CpGi,
                                      geneLine = query$geneLine, plotting = TRUE)
      }
      if(query$DataType == 'expression'){
        regplot <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(),
                                      geneLine = query$geneLine, plotting = TRUE)
      }
      print(regplot)
      dev.off()
    }
  )
  output$downloadMeanPlotPDF <- downloadHandler(
    filename = function() { paste0("Wanderer_", (toupper(query$Gene)), '_Mean_', query$DataType, '_', query$TissueType, '_', Sys.Date(), '.pdf') },
    content = function(file) {
      pdf(file, width = 12, height = 5)
      if(query$DataType == 'methylation'){
        regplot <- stat_analysis_meth(results_filt = datamethfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(), CpGislands = query$CpGi,
                                      geneLine = query$geneLine, plotting = TRUE)
      }
      if(query$DataType == 'expression'){
        regplot <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = (toupper(query$Gene)),
                                      geneNamesType = geneFormat(),
                                      geneLine = query$geneLine, plotting = TRUE)
      }
      print(regplot)
      dev.off()
    }
  )

  cancel.onSessionEnded <- session$onSessionEnded(function() {
      dbDisconnect(con)
  })

  
  #on.exit(dbDisconnect(con), add = TRUE)
  
})
