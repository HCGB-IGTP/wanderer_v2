#!/usr/bin/env R
#
## @package wanderer
## @author Anna Diez
## @author Izaskun Mallona

library(shiny)
library(RPostgreSQL)

# the file containing the db parameters
SRC <- '/imppc/labs/maplab/imallona/src/regional_profiler/wanderer'
#SRC <- '/imppc/labs/maplab/share/izaskun2anna/wanderer/current'
DB_CONF <- file.path(SRC, 'db.txt')

source(file.path(SRC, 'region_profile_methylation.R'))
source(file.path(SRC, 'GeneSize_methylation.R'))
source(file.path(SRC, 'expression_data.R'))
source(file.path(SRC, 'methylation_data.R'))
source(file.path(SRC, 'wanderer_expression.R'))
source(file.path(SRC, 'wanderer_methylation.R'))
source(file.path(SRC, 'database.R'))
source(file.path(SRC, 'max_sample.R'))
source(file.path(SRC, 'help_messages.R'))


sample_size <- read.table(file.path(SRC, "samplesN_filtered.csv"), sep = ",", stringsAsFactors = FALSE, header = TRUE)

shinyServer(function(input, output, session){
  ## db_conf <- get_db_parameters(DB_CONF)
  ## con <- dbConnect(drv, host=dbhost, port=dbport, dbname=dbname, user=dbuser, password=dbpass)
    
  con <- db_connect(DB_CONF)


  ## by default values start
  user_values <- reactiveValues(Gene = 'BRCA1',
                                GeneFormat = 'genename',
                                Tissue = 'brca')

  user_values$Gene <-'BRCA1'
  user_values$GeneFormat <- 'genename'
  user_values$Tissue <- 'brca'

 
  meth_data_curr_query <<- methylation_data(con = con,
                                            geneName = 'BRCA1',
                                            geneNamesType = 'genename',
                                            tissue = 'brca')
  
  expr_data_curr_query <<- expression_data(con = con,
                                           geneName ='BRCA1',
                                           geneNamesType = 'genename',
                                           tissue = 'brca')

  ## by default values end

  ## current query data fetch from database start

  
  observe({        
      if(input$goButton > 0) {
            user_values$Gene <- isolate(toupper(input$Gene))
            user_values$GeneFormat <- isolate(input$GeneFormat)

            meth_data_curr_query <<- methylation_data(con = con,
                                                     geneName = user_values$Gene,
                                                     geneNamesType = user_values$GeneFormat,
                                                     tissue = user_values$Tissue)
            ## stop(meth_data_curr_query)
            
            expr_data_curr_query <<- expression_data(con = con,
                                                    geneName = user_values$Gene,
                                                    geneNamesType = user_values$GeneFormat,
                                                    tissue = user_values$Tissue)
        }
    })
  
  ## current query data fetch from database end

  output$Tissues <- renderUI({
    tissues <- paste0("'",sample_size[,1], " (", sample_size[,2], ")", "'='", sample_size[,3], "'")
    tissues <- paste0(tissues,collapse=",")
    tissues <- paste0("c(",tissues,")")
    tissues <- eval(parse(text=tissues))
    selectInput("TissueType", label = h5("Project:"), choices = tissues, selected = "brca")
  })
   
  output$nNmax <- renderUI({

    if (!is.null(input$DataType) & !is.null(input$TissueType)) {
        maxn <- max_sample(sample_size, input$DataType, input$TissueType)[1]
        valor <- min(maxn, 30)
        conditionalPanel("input.plotmean == false", numericInput("nN", h5(paste0("Number of normal samples to plot (max = ", maxn, ")"), help_popup('Number of normal samples to plot')), value = valor, min = 1, max = maxn))
    }
  })
    
  output$nTmax <- renderUI({
      if (!is.null(input$DataType) & !is.null(input$TissueType)) {
          maxt <- max_sample(sample_size, input$DataType, input$TissueType)[2]
          valor <- min(maxt, 30)
          conditionalPanel("input.plotmean == false", numericInput("nT", h5(paste0("Number of tumoral samples to plot (max = ", maxt, ")"), help_popup('Number of tumoral samples to plot')), value = valor, min = 1, max = maxt))
      }
  })
  
  output$ZoomControl <- renderUI({
    ## lengthGene <- GeneSize_methylation(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat)
    lengthGene <- GeneSize_methylation(con = con,
                                       geneName = user_values$Gene,
                                       geneNamesType = user_values$GeneFormat)
    
    minval <- (lengthGene/2)-5
    sliderInput("Zoom", label = h5("Zoom"), value = 0, min = 0, max = minval, step = 100)
  })
  
  output$WalkControl <- renderUI({
    ## lengthGene <- GeneSize_methylation(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat)
      lengthGene <- GeneSize_methylation(con = con,
                                         geneName = user_values$Gene,
                                         geneNamesType = user_values$GeneFormat)

   val <- (lengthGene/2)
    if(input$DataType == 'methylation'){
      maxval <- (lengthGene/2) + 20000
      minval <- -(lengthGene/2) - 20000
    }
    if(input$DataType == 'expression'){
      maxval <- lengthGene/2
      minval <- -lengthGene/2
    }
    sliderInput("Walk", label = h5("Slider"), value = 0, min = minval, max = maxval, step = 100)
  })


  output$plot1 <- renderPlot({
     if(input$DataType == 'methylation'){
      if (!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT)) {
          wanderer_methylation(results = meth_data_curr_query,geneName = as.character(user_values$Gene),
                                     geneNamesType = as.character(user_values$GeneFormat),
                                     sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom,
                                     walk = input$Walk ,npointsN = input$nN, npointsT = input$nT,
                                     CpGislands = input$CpGi, plotmean = input$plotmean,
                                     plotting = TRUE, geneLine = input$geneLine)
          
      }}
    if(input$DataType == 'expression'){
      wanderer_expression(results = expr_data_curr_query,
                                ## geneName = input$Gene, geneNamesType = input$GeneFormat,
                                geneName = user_values$Gene,
                                geneNamesType = user_values$GeneFormat,
                                sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom,
                                walk = input$Walk ,npointsN = input$nN, npointsT = input$nT,
                                plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
    }

  }, height = 1000, width = 1000)

  ## test end
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0("Wanderer", input$Gene, '_', input$DataType, '_', input$TissueType, '_', Sys.Date(), '.png') },
    content = function(file) {
      png(file, width = 1000, height = 1000)
      if(input$DataType == 'methylation'){
        regplot <- wanderer_methylation(results = meth_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
      }
      if(input$DataType == 'expression'){
        regplot <- wanderer_expression(results = expr_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
      }
      print(regplot)
      dev.off()
    }
  )
  output$downloadNData <- downloadHandler(
    filename = function() { paste0("Wanderer", input$Gene, '_', input$DataType, '_', input$TissueType, '_Normal_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation'){
        results <- wanderer_methylation(results = meth_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[1]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      if(input$DataType == 'expression'){
        results <- wanderer_expression(results = expr_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[1]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
  output$downloadTData <- downloadHandler(
    filename = function() { paste0("Wanderer", input$Gene, '_', input$DataType, '_', input$TissueType, '_Tumor_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation'){
        results <- wanderer_methylation(results = meth_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[2]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      if(input$DataType == 'expression'){
        results <- wanderer_expression(results = expr_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[2]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
  output$downloadPData <- downloadHandler(
    filename = function() { paste0("Wanderer", input$Gene, '_', input$DataType, '_', input$TissueType, '_annotations_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation'){
        results <- wanderer_methylation(results = meth_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[3]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      if(input$DataType == 'expression'){
        results <- wanderer_expression(results = expr_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[3]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
    #}
  
  #on.exit(dbDisconnect(con), add = TRUE)
  
})
