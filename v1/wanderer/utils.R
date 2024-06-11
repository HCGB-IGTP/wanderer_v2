#!/usr/bin/env R


## CEST date formatter
cest_timestamp <- function(){
    return(format(Sys.time(), "%a_%b_%d_%T_%Z_%Y"))
}

## receives an epoch time in miliseconds
client_timestamp <- function(z){
    return(format( as.POSIXct(as.numeric(z)/1000, origin = "1970-01-01"), "%b_%d_%Y_at_%H%M%S_%Z"))

}
