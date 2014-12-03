#!/usr/bin/env R

library(shiny)
source('help_messages.R')


WLEFT = 3    # left column width
OLEFT = 1    # left column offset
WRIGHT = 7   # rigth column width
ORIGHT = 0   # rigth column offset

shinyUI(
  #corner_element = HTML(paste0('<a href=',shQuote(paste0("http://maplab.imppc.org/",sep='')), '>', 'maplab.cat', '</a>')),
  
  navbarPage(id = "Region Profile (beta)", title="",
             header = div(img(src="regional_profile.png",style="float:left; padding: 100px 20px 100px 150px;"), h2("TCGA Wanderer (beta)", style="color:#dd4814; float:left; padding: 120px 20px 20px 20px;")),
             #padding: up right bottom left
             theme = 'united.css',
             #tabPanel(a(href="http://maplab.imppc.org/", "maplab.cat")),
             #tabPanel(HTML("</a></li><li><a href=\"http://maplab.imppc.org/\">maplab.cat")),
             #tabPanel("maplab.cat",includeHTML("http://maplab.imppc.org/")),
            tabPanel("TCGA Wanderer",

                      fluidRow(
                        column(width = WLEFT, offset = OLEFT,
                              actionButton("goButton", "Go!"),
                              hr(),
                              uiOutput("Tissues"),
                               hr(),
                               selectInput("DataType", label = h5("Select Data Type:"), 
                                           choices = c("450k Methylation Array" = 'methylation', "Illumina HiSeq RNAseq" = 'expression'), selected = "methylation"),
                               hr(),
                               selectInput("GeneFormat", label = h5("Select Gene Format:"), 
                                           choices = c("Ensembl Gene ID" = 'emsemblgeneid', "Gene Name" = 'genename')),
                               hr(),
                               textInput("Gene", h5("Gene"), value='ENSG00000012048'),
                               helpText("Examples: BRCA1 or ENSG00000012048"),
                               hr(),
                               checkboxInput("plotmean", label = h5("Boxplot Summary"), FALSE),
                               hr(),
                               checkboxInput("geneLine", label = h5("Show Gene", help_popup('Show Gene')), TRUE),
                               hr(),
                               conditionalPanel("input.DataType == 'methylation'", checkboxInput("CpGi", label = h5("Show CpG islands", help_popup('Show CpG islands')), TRUE)),
                               hr(),
                               uiOutput("nNmax"),
                               hr(),
                               uiOutput("nTmax"),
                               hr(),
                               actionButton("goButton", "Go!")
                        ),
                        column( width = WRIGHT,
                                inputPanel(downloadButton('downloadPlot', 'Download Plot'), br(), downloadButton('downloadNData', 'Download Normal Data'), br(), downloadButton('downloadTData', 'Download Tumor Data'), br(), downloadButton('downloadPData', 'Download Annotation Data')),
                                hr(),
                                inputPanel(uiOutput("ZoomControl"), br(), br(), br(), uiOutput("WalkControl")),
                                hr(),
                                plotOutput(outputId = 'plot1', height = "100%", width = "100%")
                        )
                      )),
             tabPanel("Download")
  ))