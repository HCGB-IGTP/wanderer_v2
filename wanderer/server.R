#!/usr/bin/env R
#
## @package wanderer
## @author Anna Diez
## @author Izaskun Mallona

library(shiny)
library(RPostgreSQL)

# the file containing the db parameters
#SRC <- '/imppc/labs/maplab/imallona/src/regional_profiler/wanderer'
#SRC <- '/imppc/labs/maplab/share/izaskun2anna/wanderer/current'
SRC <- '/imppc/labs/maplab/adiez/region_profile/wanderer/current'
DB_CONF <- file.path(SRC, 'db.txt')

source(file.path(SRC, 'GeneSize.R'))
source(file.path(SRC, 'expression_data.R'))
source(file.path(SRC, 'methylation_data.R'))
source(file.path(SRC, 'wanderer_expression.R'))
source(file.path(SRC, 'wanderer_methylation.R'))
source(file.path(SRC, 'database.R'))
source(file.path(SRC, 'max_sample.R'))


sample_size <- read.table(file.path(SRC, "samplesN_filtered.csv"), sep = ",", stringsAsFactors = FALSE, header = TRUE)

shinyServer(function(input, output, session){
  
  con <- db_connect(DB_CONF)
  
  
  ## by default values start
  user_values <- reactiveValues(Gene = 'BRCA1',
                                GeneFormat = 'genename',
                                Tissue = 'brca')
  
  
  meth_data_curr_query <<- methylation_data(con = con,
                                            geneName = 'BRCA1',
                                            geneNamesType = 'genename',
                                            tissue = 'brca')
  
  expr_data_curr_query <<- expression_data(con = con,
                                           geneName ='BRCA1',
                                           geneNamesType = 'genename',
                                           tissue = 'brca')
  lengthGene <<- GeneSize(con = con, geneName = 'BRCA1', geneNamesType = 'genename')
  
  
  ## by default values end
  
  ## current query data fetch from database start
  
  
  observe({        
    if(input$goButton > 0) {
      user_values$Gene <- isolate(toupper(input$Gene))
      user_values$GeneFormat <- isolate(input$GeneFormat)
      
      meth_data_curr_query <<- methylation_data(con = con,
                                                geneName = user_values$Gene,
                                                geneNamesType = user_values$GeneFormat,
                                                tissue = input$TissueType)
      
      expr_data_curr_query <<- expression_data(con = con,
                                               geneName = user_values$Gene,
                                               geneNamesType = user_values$GeneFormat,
                                               tissue = input$TissueType)
      
      lengthGene <<- GeneSize(con = con, geneName = user_values$Gene, geneNamesType = user_values$GeneFormat)
      
      
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
    if (lengthGene[[1]]!=0 & lengthGene[[1]]!=1 & !is.null(input$DataType) & !is.null(input$TissueType) & !is.null(user_values$Gene) & !is.null(user_values$GeneFormat)) {
      maxn <- max_sample(sample_size, input$DataType, input$TissueType)[1]
      valor <- min(maxn, 30)
      minn <- 1
      conditionalPanel("input.plotmean == false", numericInput("nN", h5(paste0("Number of normal samples to plot (max = ", maxn, ")")), value = valor, min = minn, max = maxn))
    }
  })
  
  output$nTmax <- renderUI({
    if (lengthGene[[1]]!=0 & lengthGene[[1]]!=1 & !is.null(input$DataType) & !is.null(input$TissueType) & !is.null(user_values$Gene) & !is.null(user_values$GeneFormat)) {
      maxt <- max_sample(sample_size, input$DataType, input$TissueType)[2]
      valor <- min(maxt, 30)
      mint <- 1
      conditionalPanel("input.plotmean == false", numericInput("nT", h5(paste0("Number of tumoral samples to plot (max = ", maxt, ")")), value = valor, min = mint, max = maxt))
    }
  })
  
  output$ZoomControl <- renderUI({
    
    if(lengthGene[[1]]==0) stop(paste0("The gene ", user_values$Gene, " is not in the exon and probe annotation."))
    
    if(lengthGene[[1]]==1) stop(paste0("The gene ", user_values$Gene, " appears more than once in the genome. Please introduce an Ensembl (ENSG) identifier instead."))
    
    if(lengthGene[[1]]!=0 & lengthGene[[1]]!=1 & !is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(user_values$Gene) & !is.null(user_values$GeneFormat)) {
      
      sGene <- lengthGene[[1]]
      eGene <- lengthGene[[2]]
      
      if(input$DataType == 'methylation'){
        minGene <- sGene - 100000
        maxGene <- eGene + 100000
        tcks <- round(c(seq(minGene, sGene,(sGene-minGene)/3), seq(eGene, maxGene, (maxGene-eGene)/3)),0)
      }
      if(input$DataType == 'expression'){
        minGene <- sGene
        maxGene <- eGene
        tcks <- round(seq(sGene, eGene, (eGene-sGene)/5),0)
      }
      sliderInput("Zoom", label = h5("Zoom"), min = minGene, max = maxGene, value = c(sGene, eGene), step = 1000, ticks = tcks, width = "800px")
      
    }
  })
  
  output$regionlimit <- renderUI({
    if(!is.null(input$Zoom))  conditionalPanel("input.region == true", helpText(paste0("Define a start and an end within the slider's values (min = ", input$Zoom[1],"; max = ", input$Zoom[2])))
  })
  
  output$start <- renderUI({
    if(!is.null(input$Zoom))  conditionalPanel("input.region == true", numericInput("start", "Start", value = input$Zoom[1], min = as.numeric(input$Zoom[1]), max = as.numeric(input$Zoom[2])))
  })
  
  output$end <- renderUI({
    if(!is.null(input$Zoom))  conditionalPanel("input.region == true", numericInput("end", "End", value = input$Zoom[2], min = as.numeric(input$Zoom[1]), max = as.numeric(input$Zoom[2])))
  })
  
  
  output$plot1 <- renderPlot({
    if(lengthGene[[1]]!=0 & lengthGene[[1]]!=1 & !is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(input$Zoom) & !is.null(user_values$Gene) & !is.null(user_values$GeneFormat)) {
      
      if(input$region & input$end > input$Zoom[2]) stop(print(paste0("Region end must be less than ",input$Zoom[2])))
      if(input$region & input$start < input$Zoom[1]) stop(print(paste0("Region start must be greater than ",input$Zoom[1])))
      
      
      if(input$DataType == 'methylation'){
        if(meth_data_curr_query$empty){
          stop(print(paste0("The gene ", user_values$Gene, " does not correspond to any ", meth_data_curr_query[[2]], " or is not in the probe annotation")))
        }else{
          wanderer_methylation(results = meth_data_curr_query,geneName = as.character(user_values$Gene),
                               geneNamesType = as.character(user_values$GeneFormat),
                               sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end),
                               npointsN = input$nN, npointsT = input$nT,
                               CpGislands = input$CpGi, plotmean = input$plotmean,
                               plotting = TRUE, geneLine = input$geneLine)
        }
      }
    }
    if(input$DataType == 'expression'){
      if(expr_data_curr_query$empty){
        stop(print(paste0("Your gene ", user_values$Gene, " does not correspond to any ", expr_data_curr_query[[2]], " or is not in the exon annotation")))
      }else{
        wanderer_expression(results = expr_data_curr_query,
                            geneName = user_values$Gene,
                            geneNamesType = user_values$GeneFormat,
                            sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end),
                            npointsN = input$nN, npointsT = input$nT,
                            plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
      }
      
    }
    
  }, height = 1000, width = 1000)
  
  ## test end
  output$downloadPlot <- downloadHandler(

    filename = function() { paste0("Wanderer_", input$Gene, '_', input$DataType, '_', input$TissueType, '_', Sys.Date(), '.png') },
    content = function(file) {
      png(file, width = 1000, height = 1000)
      if(input$DataType == 'methylation'){
        regplot <- wanderer_methylation(results = meth_data_curr_query, geneName = input$Gene, 
                                        geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, 
                                        zoom = c(input$start, input$end), npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, 
                                        plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
      }
      if(input$DataType == 'expression'){
        regplot <- wanderer_expression(results = expr_data_curr_query, geneName = input$Gene, 
                                       geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, 
                                       zoom = c(input$start, input$end), npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, 
                                       plotting = TRUE, geneLine = input$geneLine)
      }
      print(regplot)
      dev.off()
    }
  )
  output$downloadPlotPDF <- downloadHandler(

    filename = function() { paste0("Wanderer_", input$Gene, '_', input$DataType, '_', input$TissueType, '_', Sys.Date(), '.pdf') },
    content = function(file) {
      pdf(file, width = 10, height = 13)
      if(input$DataType == 'methylation'){
        regplot <- wanderer_methylation(results = meth_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat, 
                                        sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end), npointsN = input$nN,
                                        npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
      }
      if(input$DataType == 'expression'){
        regplot <- wanderer_expression(results = expr_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat,
                                       sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end), npointsN = input$nN, 
                                       npointsT = input$nT, plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
      }
      print(regplot)
      dev.off()
    }
  )
  output$downloadNData <- downloadHandler(

    filename = function() { paste0("Wanderer_", input$Gene, '_', input$DataType, '_', input$TissueType, '_Normal_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation'){
        results <- wanderer_methylation(results = meth_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat,
                                        sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end), npointsN = input$nN, npointsT = input$nT,
                                        CpGislands = input$CpGi, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[1]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      if(input$DataType == 'expression'){
        results <- wanderer_expression(results = expr_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat,
                                       sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end), npointsN = input$nN, npointsT = input$nT,
                                       plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[1]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
  output$downloadTData <- downloadHandler(

    filename = function() { paste0("Wanderer_", input$Gene, '_', input$DataType, '_', input$TissueType, '_Tumor_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation'){
        results <- wanderer_methylation(results = meth_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat,
                                        sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end), npointsN = input$nN, npointsT = input$nT,
                                        CpGislands = input$CpGi, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[2]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      if(input$DataType == 'expression'){
        results <- wanderer_expression(results = expr_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat,
                                       sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end), npointsN = input$nN, npointsT = input$nT,
                                       plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[2]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
  output$downloadPData <- downloadHandler(

    filename = function() { paste0("Wanderer_", input$Gene, '_', input$DataType, '_', input$TissueType, '_annotations_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation'){
        results <- wanderer_methylation(results = meth_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat,
                                        sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end), npointsN = input$nN, npointsT = input$nT,
                                        CpGislands = input$CpGi, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[3]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      if(input$DataType == 'expression'){
        results <- wanderer_expression(results = expr_data_curr_query, geneName = input$Gene, geneNamesType = input$GeneFormat,
                                       sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end), npointsN = input$nN, npointsT = input$nT,
                                       plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[3]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
  #}
  
  #on.exit(dbDisconnect(con), add = TRUE)
  
})
