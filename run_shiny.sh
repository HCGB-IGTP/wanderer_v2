#!/bin/bash

R3=/soft/general/R-3.1.0/bin/R
#$R3 -e "shiny::runApp('/imppc/labs/maplab/adiez/region_profile/web/', launch.browser = TRUE, port = 5559)"

$R3 -e "shiny::runApp('/imppc/labs/maplab/imallona/src/regional_profiler/wanderer', launch.browser = TRUE, port = 5559)"
