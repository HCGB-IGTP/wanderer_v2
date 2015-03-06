#!/usr/bin/env R
#
## @package wanderer
## @author Izaskun Mallona

library(RPostgreSQL)

# the file containing the db parameters

#SRC <- '/imppc/labs/maplab/imallona/src/regional_profiler/wanderer'
#SRC <- '/data/shiny/apps/correlational'
## DB_CONF <- file.path(SRC, 'db.txt')

# returns a unique connection
db_connect <- function(db_conf_fn) {
    drv <- dbDriver("PostgreSQL", max.con = 100)

    db_conf <- get_db_parameters(db_conf_fn)
    
    con <- dbConnect(drv,
                     user = db_conf[['user']],
                     password = db_conf[['password']],
                     dbname = db_conf[['dbname']],
                     host = db_conf[['host']],
                     port = db_conf[['port']])

    ## on.exit(dbDisconnect(con), add = TRUE)
    return(con)

}

## sends a stmt and retrieves all the result
get_query <- function(con, stmt) {
    query <- dbSendQuery(con, stmt)
    return(fetch(query, n = -1))
}

get_db_parameters <- function(conf) {
  params <- read.table(conf, sep = ",", stringsAsFactors = FALSE)
  return(list(user = params$V2[1],
              ## password = params$V2[2],
              dbname = params$V2[3],
              host = params$V2[4],
              port = params$V2[5]))
}

sql_quote <- function(x, quote = "'") {
  y <- gsub(quote, paste0(quote, quote), x, fixed = TRUE)
  y <- paste0(quote, y, quote)
  y[is.na(x)] <- "NULL"
  names(y) <- names(x)  
  return(y)
}
