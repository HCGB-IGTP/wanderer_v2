#!/usr/bin/env R

PLACEMENT = 'right'
TRIGGER = 'hover'
  

# the help dictionary containing pairs of (key, values)
# the key is the same as the div associated to the action that must be documented
hd <- list()

## help messages start
hd[['Number of bp downstream']] <- 'Number of base pairs downstream the gene to be plot. Only allowed for Illumina 450k methylation data.'
hd[['Number of bp upstream']] <- 'Number of base pairs downstream the gene to be plot. Only allowed for Illumina 450k methylation data.'
hd[['Show all region']] <- 'If checked, the plot shows the region between the start of the gene - downstream bp and the end of the gene + upstream bp. If not checked, the plot only shows the parts of the region that includes a CpG site. Only allowed for Illumina 450k methylation data.'
hd[['Show CpG islands']] <- 'CpG islands are displayed in green color in the x-axis. Only allowed for Illumina 450k methylation data.'
hd[['Show Gene']] <- 'Gene is displayed as an orange arrow in the bottom of the plot.'
hd[['Number of normal samples to plot']] <- 'Choose the number of Normal samples to plot. The maximum number of Normal samples to plot for the choosen tissue is indicated in brackets.'
hd[['Number of tumoral samples to plot']] <- 'Choose the number of Tumoral samples to plot. The maximum number of Tumoral samples to plot for the choosen tissue is indicated in brackets.'
## help messages end



## popups end

# gist from jcheng5
# https://gist.github.com/jcheng5/5913297
help_popup_core <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual')) {
  tagList(
          singleton(
                    tags$head(
                              tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })")
                              )
                    ),
          tags$a(
                 href = "#", class = "btn btn-mini", `data-toggle` = "popover",
                 title = "", `data-content` = content, `data-animation` = TRUE,
                 `data-placement` = match.arg(placement, several.ok=TRUE)[1],
                 `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
                 #tags$i(class="icon-question-sign")
                 tags$b("?")
                 )
          )
}


# builds the popup querying the help dictionary
# @param key the key to fetch the hd dictionary values
help_popup <- function(key) {
   return(help_popup_core(key,
                          hd[[key]],
                          placement = PLACEMENT,
                          trigger= TRIGGER))
   
}

help_paragraph <- function(key) {
  return(hd[[key]])
}
