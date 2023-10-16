# DTedit Shiny Demo
library(shiny)
library(RSQLite)
library(DTedit)

##### Load books data.frame as a SQLite database
enssentials = c("chr","gene","start","end","type")
mat <- as.data.frame(matrix(data = rep(0,10), nrow = 2))
colnames(mat) <- enssentials

con <- dbConnect(RSQLite::SQLite(), "gene_tab.sqlite")

if(!'gene_tab' %in% dbListTables(con)) {
  genetable <- .load.gff("/Users/qianli/Desktop/workfiles/gene_circuits/depthfile/control.cir.gff")
  gene_tab <- genetable[, match(enssentials, colnames(genetable))]
  gene_tab$id <- 1:length(gene_tab$chr)
  dbWriteTable(con, "gene_tab", gene_tab, overwrite = TRUE)
}


#fetching information from SQLite objects
getgenetab <- function() {
  res <- dbSendQuery(con, "SELECT * FROM gene_tab")
  gene_tab <- dbFetch(res)
  dbClearResult(res)
  return(gene_tab)
}

## Callback functions.
#insert a new line
genetab.insert.callback <- function(data, row) {
  query <- paste0("INSERT INTO gene_tab (chr, gene, start, end, type, id) VALUES (",
                  "", as.character(data[row,]$chr), ", ",
                  "'", as.character(data[row,]$gene), "', ",
                  "'", as.numeric(data[row,]$start), "', ",
                  "'", as.numeric(data[row,]$end), "', ",
                  "'", as.character(data[row,]$type), "', ",
                  "'",  max(getgenetab()$id) + 1, "' ",
                  ")")
  print(query) # For debugging
  dbSendQuery(con, query)
  return(getgenetab())
}

#updata cells
genetab.update.callback <- function(data, olddata, row) {
  query <- paste0("UPDATE gene_tab SET ",
                  "chr = '", as.character(data[row,]$gene), "', ",
                  "gene = '", as.character(data[row,]$gene), "', ",
                  "start = '", as.numeric(data[row,]$start), "', ",
                  "end = '", as.numeric(data[row,]$end), "' ",
                  "type = '", as.character(data[row,]$type), "' ",
                  "WHERE id = ", data[row,]$id)
  print(query) # For debugging
  dbSendQuery(con, query)
  return(getgenetab())
}

#delete lines
genetab.delete.callback <- function(data, row) {
  query <- paste0('DELETE FROM gene_tab WHERE id = ', data[row,]$id)
  dbSendQuery(con, query)
  return(getgenetab())
}

##### Create the Shiny server
server <- function(input, output, session) {
  gene_tab <- getgenetab()
  
  dtedit(input, output,
         name = 'my_datatable',
         thedata = gene_tab,
         edit.cols = c('chr', 'gene', 'start', 'end', 'type'),
         edit.label.cols = c('chr', 'gene', 'start', 'end', 'type'),
         input.types = c(chr='textAreaInput', 
                         gene='textAreaInput', 
                         start='numericInput',
                         end='numericInput',
                         type = 'textAreaInput'),
         # view.cols = names(genetab)[1:6],
         callback.update = genetab.update.callback,
         callback.insert = genetab.insert.callback,
         callback.delete = genetab.delete.callback)
}

##### Create the shiny UI
ui <- fluidPage(
  h3('Books'),
  uiOutput('my_datatable'),

)

shinyApp(ui = ui, server = server)
