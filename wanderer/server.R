
library(shiny)
library(RPostgreSQL)
library(limma)

# the file containing the db parameters
SRC <- '/imppc/labs/maplab/imallona/src/regional_profiler/wanderer'
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
  
  #################################################
  #detect gene format
  geneFormat <- reactive({
    if(input$goButton == 0) {
      GeneFormat <- "genename"
    } else{
      if(substr(isolate(toupper(input$Gene)), 1, 4) == "ENSG"){
        GeneFormat <- "emsemblgeneid"
      } else{
        GeneFormat <- "genename"
      }
      ## updateSliderInput(session, 'Zoom',  value = c(sGene, eGene))
    }
  })

  ## test start

  rv <- reactiveValues(geneSize = NULL, geneNamesType_label = NULL,
                       eGene = NULL, sGene = NULL, minGene = NULL, maxGene = NULL)
  
  observe({

      if(input$goButton == 0) {
          geneSize <- GeneSize(con = con, geneName = 'BRCA1', geneNamesType = 'genename')   
      } else{
          geneSize <- GeneSize(con = con, geneName = isolate(toupper(input$Gene)), geneNamesType = geneFormat())  
      }
      
      ## geneSize <- GeneSize(con = con, geneName = isolate(toupper(input$Gene)), geneNamesType = geneFormat())
      if (geneFormat() == "genename")
          geneNamesType_label <- "Gene Name"
      if (geneFormat() == "emsemblgeneid")
          geneNamesType_label <- "Ensembl Gene ID"
      
      if(geneSize[[1]]==0)
          stop(paste0("The gene ", isolate(toupper(input$Gene)), " is not in the ", geneNamesType_label, " annotation."))
      
      if(geneSize[[1]]==1)
          stop(paste0("The gene ", isolate(toupper(input$Gene)), " appears more than once in the genome. Please introduce an Ensembl (ENSG) identifier instead."))
      
      if(geneSize[[1]]!=0 & geneSize[[1]]!=1){
          
          sGene <- geneSize[[1]]
          eGene <- geneSize[[2]]
          sGene <- ((sGene%/%1000)-1)*1000
          eGene <- ((eGene%/%1000)+1)*1000
          
          
          if(input$DataType == 'methylation'){
              minGene <- sGene - 100000
              maxGene <- eGene + 100000
              tcks <- c(seq(minGene, sGene,(sGene-minGene)/3), seq(eGene, maxGene, (maxGene-eGene)/3))
              tcks <- (tcks%/%1000)*1000
              tcks <- format(tcks, big.mark = ',')        
          }
          else if (input$DataType == 'expression'){
              minGene <- sGene
              maxGene <- eGene
              tcks <- seq(sGene, eGene, (eGene-sGene)/5)
              tcks <- (tcks%/%1000)*1000
              tcks <- format(tcks, big.mark = ',')        
          }
      }

      ## rv$geneSize <- geneSize
      ## rv$geneNamesType_label <- geneNamesType_label
      ## rv$eGene <- eGene
      ## rv$sGene <- sGene
      ## rv$minGene <- minGene
      ## rv$maxGene <- maxGene
      ## rv$tcks <- tcks
      
      ## rv <<- rv

      geneSize <<- geneSize
      geneNamesType_label <<- geneNamesType_label
      eGene <<- eGene
      sGene <<- sGene
      minGene <<- minGene
      maxGene <<- maxGene
      tcks <<- tcks
      
      ## rv <<- rv

      
  })

  
  ## observe({

  ##     if(input$goButton > 0) {
  ##         updateSliderInput(session, 'Zoom',  value = c(rv$sGene, rv$eGene))
          
  ##     }
  ## })

  ## test end
  
  
  #################################################
  #zoom
  output$ZoomControl <- renderUI({
    
  ##   if(input$goButton == 0) {
  ##     geneSize <- GeneSize(con = con, geneName = 'BRCA1', geneNamesType = 'genename')   
  ## } else{
  ##     geneSize <- GeneSize(con = con, geneName = isolate(toupper(input$Gene)), geneNamesType = geneFormat())  
  ##   }
    
    
    
  ##   if (geneFormat() == "genename")
  ##       geneNamesType_label <- "Gene Name"
  ##   if (geneFormat() == "emsemblgeneid")
  ##       geneNamesType_label <- "Ensembl Gene ID"
    
  ##   if(geneSize[[1]]==0)
  ##       stop(paste0("The gene ", isolate(toupper(input$Gene)), " is not in the ", geneNamesType_label, " annotation."))
    
  ##   if(geneSize[[1]]==1)
  ##       stop(paste0("The gene ", isolate(toupper(input$Gene)), " appears more than once in the genome. Please introduce an Ensembl (ENSG) identifier instead."))
    
  ##   if(geneSize[[1]]!=0 & geneSize[[1]]!=1){
      
  ##     sGene <- geneSize[[1]]
  ##     eGene <- geneSize[[2]]
  ##     sGene <- ((sGene%/%1000)-1)*1000
  ##     eGene <- ((eGene%/%1000)+1)*1000
      
      
  ##     if(input$DataType == 'methylation'){
  ##       minGene <- sGene - 100000
  ##       maxGene <- eGene + 100000
  ##       tcks <- c(seq(minGene, sGene,(sGene-minGene)/3), seq(eGene, maxGene, (maxGene-eGene)/3))
  ##       tcks <- (tcks%/%1000)*1000
  ##       tcks <- format(tcks, big.mark = ',')        
  ##     }
  ##     else if (input$DataType == 'expression'){
  ##       minGene <- sGene
  ##       maxGene <- eGene
  ##       tcks <- seq(sGene, eGene, (eGene-sGene)/5)
  ##       tcks <- (tcks%/%1000)*1000
  ##       tcks <- format(tcks, big.mark = ',')        
  ##     }
      if(input$goButton > -1)
      sliderInput("Zoom", label = h5("Zoom"), min = minGene, max = maxGene, value = c(sGene, eGene), ticks = tcks, step = 1000, width = "800px")

      ## sliderInput("Zoom", label = h5("Zoom2"), min = minGene, max = maxGene, value = c(sGene, eGene), ticks = tcks, step = 1000, width = "800px")

      ## updateSliderInput(session, 'Zoom',  value = c(sGene, eGene))
      ## sliderInput("Zoom", label = h5("Zoom"), min = rv$minGene, max = rv$maxGene, value = c(rv$sGene, rv$eGene), ticks = rv$tcks, step = 1000, width = "800px")

  ## }
  })

  ##  #zoom
  ## output$ZoomControl <- renderUI({
    
  ##   if(input$goButton == 0) {
  ##     geneSize <- GeneSize(con = con, geneName = 'BRCA1', geneNamesType = 'genename')  
  ##   }else{
  ##     geneSize <- GeneSize(con = con, geneName = isolate(toupper(input$Gene)), geneNamesType = geneFormat())  
  ##   }
    
    
  ##   if (geneFormat() == "genename")
  ##       geneNamesType_label <- "Gene Name"
  ##   if (geneFormat() == "emsemblgeneid")
  ##       geneNamesType_label <- "Ensembl Gene ID"
    
  ##   if(geneSize[[1]]==0)
  ##       stop(paste0("The gene ", isolate(toupper(input$Gene)), " is not in the ", geneNamesType_label, " annotation."))
    
  ##   if(geneSize[[1]]==1)
  ##       stop(paste0("The gene ", isolate(toupper(input$Gene)), " appears more than once in the genome. Please introduce an Ensembl (ENSG) identifier instead."))
    
  ##   if(geneSize[[1]]!=0 & geneSize[[1]]!=1){
      
  ##     sGene <- geneSize[[1]]
  ##     eGene <- geneSize[[2]]
  ##     sGene <- ((sGene%/%1000)-1)*1000
  ##     eGene <- ((eGene%/%1000)+1)*1000
      
      
  ##     if(input$DataType == 'methylation'){
  ##       minGene <- sGene - 100000
  ##       maxGene <- eGene + 100000
  ##       tcks <- c(seq(minGene, sGene,(sGene-minGene)/3), seq(eGene, maxGene, (maxGene-eGene)/3))
  ##       tcks <- (tcks%/%1000)*1000
  ##       tcks <- format(tcks, big.mark = ',')        
  ##     }
  ##     else if (input$DataType == 'expression'){
  ##       minGene <- sGene
  ##       maxGene <- eGene
  ##       tcks <- seq(sGene, eGene, (eGene-sGene)/5)
  ##       tcks <- (tcks%/%1000)*1000
  ##       tcks <- format(tcks, big.mark = ',')        
  ##     }
      
  ##     sliderInput("Zoom", label = h5("Zoom"), min = minGene, max = maxGene, value = c(sGene, eGene), ticks = tcks, step = 1000, width = "800px")

  ##     ## sliderInput("Zoom", label = h5("Zoom2"), min = minGene, max = maxGene, value = c(sGene, eGene), ticks = tcks, step = 1000, width = "800px")

  ##     ## updateSliderInput(session, 'Zoom',  value = c(sGene, eGene))
  ##   }
  ## })

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
    if(input$DataType == 'methylation'){
      if(input$goButton == 0) {
        methylation_data(con = con, geneName = 'BRCA1', geneNamesType = 'genename', tissue = 'brca')  
      }else{
        methylation_data(con = con, geneName = isolate(toupper(input$Gene)), geneNamesType = geneFormat(), tissue = input$TissueType)
      }
    }
  })
  
  #################################################
  #Filtering Methylation data
  datamethfilt <- reactive({
    if(input$region & ((input$end > input$Zoom[2]) | (input$end < input$Zoom[1]) | (input$start < input$Zoom[1]) | (input$start > input$Zoom[2])))
        stop(print(paste0("The region must be between ", input$Zoom[1], " and ", input$Zoom[2])))
    
    if(input$DataType == 'methylation'){
      if(datameth()[['empty']]){
        stop(print(paste0("The gene ", isolate(toupper(input$Gene)), " does not correspond to any ", datameth()[['geneNamesType_label']], " or is not in the probe annotation")))
      }else{
        data_meth_filtering(results = datameth(), sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end))
      }
    }
  })
  
  
  #################################################
  #Reading Expression data
  dataexpr <- reactive({
    if(input$DataType == 'expression'){
      if(input$goButton == 0) {
        expression_data(con = con, geneName = 'BRCA1', geneNamesType = 'genename', tissue = 'brca')  
      }else{
        expression_data(con = con, geneName = isolate(toupper(input$Gene)), geneNamesType = geneFormat(), tissue = input$TissueType)
      }
    }
  })
  
  #################################################
  #Filtering Expression data
  dataexprfilt <- reactive({
    if(input$region & ((input$end > input$Zoom[2]) | (input$end < input$Zoom[1]) | (input$start < input$Zoom[1]) | (input$start > input$Zoom[2]))) stop(print(paste0("The region must be between ", input$Zoom[1], " and ", input$Zoom[2])))
    
    if(input$DataType == 'expression'){
      if(dataexpr()[['empty']]){
        stop(print(paste0("The gene ", isolate(toupper(input$Gene)), " does not correspond to any ", datameth()[['geneNamesType_label']], " or is not in the exon annotation")))
      }else{
        data_expr_filtering(results = dataexpr(), sampleSize = sample_size, tissue = input$TissueType, zoom = c(input$start, input$end))
      }
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
    if (!is.null(input$DataType) & !is.null(input$TissueType) & !is.null(isolate(toupper(input$Gene)))) {
      maxn <- max_sample(sample_size, input$DataType, input$TissueType)[1]
      valor <- min(maxn, 30)
      minn <- 1
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
    if(!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(isolate(toupper(input$Gene)))) {
      
      if(input$DataType == 'methylation'){
        if(datamethfilt()[['empty']]){
          stop(print(paste0("There are not probes in this region")))
          ## here
        }else{
          wanderer_methylation(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                               geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                               CpGislands = input$CpGi, plotmean = input$plotmean,
                               plotting = TRUE, geneLine = input$geneLine)
        }
      } else if(input$DataType == 'expression'){
        if(dataexprfilt()[['empty']]){
          stop(print(paste0("There are not exons in this region")))
        }else{
          wanderer_expression(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                              geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                              plotmean = input$plotmean, plotting = TRUE, geneLine = input$geneLine)
        }
      }
    }
  }, height = 1000, width = 1000)
  
  
  
  ##################################################
  #print summary plot
  output$plotStat <- renderPlot({
    if(!is.null(input$TissueType) & !is.null(input$nN) & !is.null(input$nT) & !is.null(isolate(toupper(input$Gene)))) {
      
      
      if(input$DataType == 'methylation'){
        if(!datamethfilt()[['empty']]){
          stat_analysis_meth(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                             geneNamesType = geneFormat(), CpGislands = input$CpGi,
                             geneLine = input$geneLine, plotting = TRUE)
        }
      } else if(input$DataType == 'expression'){
        if(!dataexprfilt()[['empty']]){
          stat_analysis_expr(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                             geneNamesType = geneFormat(),
                             geneLine = input$geneLine, plotting = TRUE)
        }
      }
    }
  }, height = 500, width = 1000)
  
  
  #################################################
  #dowload plot as png & pdf
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0("Wanderer_", isolate(toupper(input$Gene)), '_', input$DataType, '_', input$TissueType, '_', Sys.Date(), '.png') },
    content = function(file) {
      png(file, width = 1000, height = 1000)
      if(input$DataType == 'methylation'){
        regplot <- wanderer_methylation(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                        geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                        CpGislands = input$CpGi, plotmean = input$plotmean,
                                        plotting = TRUE, geneLine = input$geneLine)
      } else if(input$DataType == 'expression'){
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
      pdf(file, width = 10, height = 13)
      if(input$DataType == 'methylation'){
        regplot <- wanderer_methylation(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                        geneNamesType = geneFormat(), npointsN = input$nN, npointsT = input$nT,
                                        CpGislands = input$CpGi, plotmean = input$plotmean,
                                        plotting = TRUE, geneLine = input$geneLine)
      } else if(input$DataType == 'expression'){
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
    png(file, width = 1000, height = 500)
    if(input$DataType == 'methylation'){
      regplot <- stat_analysis_meth(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                    geneNamesType = geneFormat(), CpGislands = input$CpGi,
                                    geneLine = input$geneLine, plotting = TRUE)
    }
    if(input$DataType == 'expression'){
      regplot <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                                    geneNamesType = geneFormat(),
                                    geneLine = input$geneLine, plotting = TRUE)
    }
    print(regplot)
    dev.off()
  }
)
output$downloadMeanPlotPDF <- downloadHandler(
  filename = function() { paste0("Wanderer_", isolate(toupper(input$Gene)), '_Mean_', input$DataType, '_', input$TissueType, '_', Sys.Date(), '.pdf') },
  content = function(file) {
    pdf(file, width = 12, height = 5)
    if(input$DataType == 'methylation'){
      regplot <- stat_analysis_meth(results_filt = datamethfilt(), geneName = isolate(toupper(input$Gene)),
                                    geneNamesType = geneFormat(), CpGislands = input$CpGi,
                                    geneLine = input$geneLine, plotting = TRUE)
    }
    if(input$DataType == 'expression'){
      regplot <- stat_analysis_expr(results_filt = dataexprfilt(), geneName = isolate(toupper(input$Gene)),
                                    geneNamesType = geneFormat(), CpGislands = input$CpGi,
                                    geneLine = input$geneLine, plotting = TRUE)
    }
    print(regplot)
    dev.off()
  }
)


#on.exit(dbDisconnect(con), add = TRUE)

})
