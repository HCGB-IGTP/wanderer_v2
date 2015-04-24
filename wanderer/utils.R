#!/usr/bin/env R


## CEST date formatter
cest_timestamp <- function(){
    return(format(Sys.time(), "%a_%b_%d_%T_%Z_%Y"))
}
