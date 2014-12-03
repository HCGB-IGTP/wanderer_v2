#!/usr/bin/env R
#
## @package wanderer
# @author Anna Diez

library(shiny)
library(RPostgreSQL)

# the file containing the db parameters
## SRC <- '/imppc/labs/maplab/adiez/region_profile/web/'
SRC <- '/imppc/labs/maplab/imallona/src/regional_profiler/wanderer'
DB_CONF <- file.path(SRC, 'db.txt')

source(file.path(SRC, 'region_profile_methylation.R'))
source(file.path(SRC, 'GeneSize_methylation.R'))
source(file.path(SRC, 'region_profile_expression.R'))
source(file.path(SRC, 'database.R'))
source(file.path(SRC, 'max_sample.R'))


sample_size <- read.table(file.path(SRC, "samplesN_filtered.csv"), sep = ",", stringsAsFactors = FALSE, header = TRUE)


shinyServer(function(input, output){
  ## db_conf <- get_db_parameters(DB_CONF)
  ## con <- dbConnect(drv, host=dbhost, port=dbport, dbname=dbname, user=dbuser, password=dbpass)

  con <- db_connect(DB_CONF)

  
  user_values <- reactiveValues(Gene = NA,
                                GeneFormat = NA)
  ## user_values <- list(gene_name = 'ENSG0000001204',
  ##                     gene_format = 'emsemblgeneid')
  
  ## isolating tests start

  user_values$gene_name <-'ENSG0000001204'
  user_values$gene_format <- 'ensemblgeneid'
  
  observe({        
      if(input$goButton > 0) {
            user_values$Gene <- isolate(input$Gene)
            user_values$GeneFormat <- isolate(input$GeneFormat)
            ## stop(print(user_values))
        }
    })
  ## isolating tests end
  
    #if(exists(input$TissueType)){
      #geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType
  
  output$Tissues <- renderUI({
    tissues <- paste0("'",sample_size[,1], " (", sample_size[,2], ")", "'='", sample_size[,3], "'")
    tissues <- paste0(tissues,collapse=",")
    tissues <- paste0("c(",tissues,")")
    tissues <- eval(parse(text=tissues))
    selectInput("TissueType", label = h5("Select Tissue Type:"), choices = tissues, selected = "brca")
  })
  
  output$nNmax <- renderUI({
    maxn <- max_sample(sample_size, input$DataType, input$TissueType)[1]
    valor <- min(maxn, 30)
    conditionalPanel("input.plotmean == false", numericInput("nN", h5(paste0("Number of normal samples to plot (max = ", maxn, ")"), help_popup('Number of normal samples to plot')), value = valor, min = 1, max = maxn))
  })
    
  output$nTmax <- renderUI({
    maxt <- max_sample(sample_size, input$DataType, input$TissueType)[2]
    valor <- min(maxt, 30)
    conditionalPanel("input.plotmean == false", numericInput("nT", h5(paste0("Number of tumoral samples to plot (max = ", maxt, ")"), help_popup('Number of tumoral samples to plot')), value = valor, min = 1, max = maxt))
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

  ## ##nonisolated
  ## output$plot1_bak <- renderPlot({
  ##   ## if(input$goButton == 0) return()
  ##   ## input$goButton
  ##   ## isolate({
  ##   ## if(input$DataType == 'methylation'){
  ##   ##   print(input$nN)
  ##   ##   region_profile_methylation(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
  ##   ## }
  ##   ## if(input$DataType == 'expression'){
  ##   ##   region_profile_expression(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
  ##   ## }
  ##   ## })
  ##   if(input$DataType == 'methylation'){
  ##     print(input$nN)
  ##     print(user_values)
  ##     region_profile_methylation(con = con,
  ##                                ## geneName = input$Gene, geneNamesType = input$GeneFormat,
  ##                                geneName = user_values$geneName,
  ##                                geneNamesType = user_values$geneNamesType,
  ##                                sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom,
  ##                                walk = input$Walk ,npointsN = input$nN, npointsT = input$nT,
  ##                                CpGislands = input$CpGi, plotmean = input$plotmean,
  ##                                plotting = TRUE, geneLine = input$geneLine)
  ##   }
  ##   if(input$DataType == 'expression'){
  ##     region_profile_expression(con = con,
  ##                               ## geneName = input$Gene, geneNamesType = input$GeneFormat,
  ##                               geneName = user_values$geneName,
  ##                               geneNamesType = user_values$geneNamesType,
  ##                               sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom,
  ##                               walk = input$Walk ,npointsN = input$nN, npointsT = input$nT,
  ##                               plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
  ##   }

  ## }, height = 1000, width = 1000)


  ## test start

  output$plot1 <- renderPlot({
    ## if(input$goButton == 0) return()
    ## input$goButton
    ## isolate({
    ## if(input$DataType == 'methylation'){
    ##   print(input$nN)
    ##   region_profile_methylation(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
    ## }
    ## if(input$DataType == 'expression'){
    ##   region_profile_expression(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
    ## }
    ## })

    if(input$DataType == 'methylation'){
      print(input$nN)
      ## stop(user_values$gene_name)
      region_profile_methylation(con = con,
                                 ## geneName = input$Gene, geneNamesType = input$GeneFormat,
                                 geneName = as.character(user_values$Gene),
                                 geneNamesType = as.character(user_values$GeneFormat),
                                 sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom,
                                 walk = input$Walk ,npointsN = input$nN, npointsT = input$nT,
                                 CpGislands = input$CpGi, plotmean = input$plotmean,
                                 plotting = TRUE, geneLine = input$geneLine)
      
    }
    if(input$DataType == 'expression'){
        region_profile_expression(con = con,
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
    filename = function() { paste0(input$Gene, '_', input$DataType, '_', input$TissueType, '_', Sys.Date(), '.png') },
    content = function(file) {
      png(file, width = 1000, height = 1000)
      if(input$DataType == 'methylation'){
        regplot <- region_profile_methylation(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
      }
      if(input$DataType == 'expression'){
        regplot <- region_profile_expression(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
      }
      print(regplot)
      dev.off()
    }
  )
  output$downloadNData <- downloadHandler(
    filename = function() { paste0(input$Gene, '_', input$DataType, '_', input$TissueType, '_Normal_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation'){
        results <- region_profile_methylation(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[1]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      if(input$DataType == 'expression'){
        results <- region_profile_expression(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[1]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
  output$downloadTData <- downloadHandler(
    filename = function() { paste0(input$Gene, '_', input$DataType, '_', input$TissueType, '_Tumor_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation'){
        results <- region_profile_methylation(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[2]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      if(input$DataType == 'expression'){
        results <- region_profile_expression(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[2]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
  output$downloadPData <- downloadHandler(
    filename = function() { paste0(input$Gene, '_', input$DataType, '_', input$TissueType, '_annotations_', Sys.Date(), '.txt') },
    content = function(file) {
      if(input$DataType == 'methylation'){
        results <- region_profile_methylation(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, CpGislands = input$CpGi, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[3]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
      if(input$DataType == 'expression'){
        results <- region_profile_expression(con = con, geneName = input$Gene, geneNamesType = input$GeneFormat, sampleSize = sample_size, tissue = input$TissueType, zoom = input$Zoom, walk = input$Walk ,npointsN = input$nN, npointsT = input$nT, plotmean = input$plotmean, plotting = FALSE, geneLine = input$geneLine)
        write.table(results[[3]], file = file, sep = "\t", row.names = FALSE, quote = FALSE)        
      }
    }
  )
    #}
  
  #on.exit(dbDisconnect(con), add = TRUE)
  
})
