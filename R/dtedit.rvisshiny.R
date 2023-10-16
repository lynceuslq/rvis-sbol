#!/bin/Rscript
#author: Qian Li
#contact: 402146079@qq.com
#An interactive app to load and plot genetic circuits based on SBOL3 standards

library(shiny)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(ggpubr)
library(ggtree)
library(aplot)
library(ggfortify)
library(ggplotify)
library(easyGgplot2)
library(ComplexHeatmap)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
#library(reticulate)
library(gggenes)
library(ggfittext)
library(ggrepel)
library(plasmapR)
#library(caret)
library(dplyr)
library(DT)
library(RSQLite)
library(DTedit)

#loading and preprocess gff files
.load.gff <-  function(file) {
  gfff <- read.table(file, header=F, sep="\t")
  geneticparts <- c("gene", "promoter", "terminator")
  ribo <- c("RBS", "ribozyme")
  unit <- c("transcript", "promoter_unit")
  
  for(i in 1:length(gfff$V1)){
    # gfff$start[i] <- min(gfff$V4[i], gfff$V5[i])
    # gfff$end[i] <- max(gfff$V4[i], gfff$V5[i])
    gfff$start[i] <- gfff$V4[i]
    gfff$end[i] <- gfff$V5[i]
    gfff$gene[i] <- unlist(strsplit(gfff$V9[i], "="))[2]
    
  }
  
  gfff$strand <- c("s")
  
  if(length(gfff[gfff$V7 == "+", ]$strand) > 0 ) {
    gfff[gfff$V7 == "+", ]$strand <- c("forward")
  }
  if(length(gfff[gfff$V7 == "-", ]$strand) > 0) {
    gfff[gfff$V7 == "-", ]$strand <- c("reverse")
    ls <- gfff[gfff$V7 == "-", ]$start
    le <- gfff[gfff$V7 == "-", ]$end
    gfff[gfff$V7 == "-", ]$start <- le
    gfff[gfff$V7 == "-", ]$end <- ls
  }
  
  
  gfff$molecule <- gfff$V1
  
  if(length(gfff[gfff$V3 %in% unit, ]$molecule) > 0) {
    gfff[gfff$V3 %in% unit, ]$molecule <- c("promoter unit")
  }
  
  if(length(gfff[gfff$V3 %in% geneticparts, ]$molecule) > 0) {
    gfff[gfff$V3 %in% geneticparts, ]$molecule <- c("genetic elements")
  }
  
  if(length(gfff[gfff$V3 %in% ribo, ]$molecule) > 0) {
    gfff[gfff$V3 %in% ribo, ]$molecule <- c("ribozyme binding sites")
  }
  
  gfff$direction <- c("s")
  if(length(gfff[gfff$V7 == "+", ]$strand) > 0 ) {
    gfff[gfff$V7 == "+", ]$direction <- c("LEFT")
  }
  
  if(length(gfff[gfff$V7 == "-", ]$strand) > 0) {
    gfff[gfff$V7 == "-", ]$direction <- c("RIGHT")
  }
  
  gene_tab <- data.frame(molecule=gfff$molecule, 
                         chr = gfff$V1,
                         gene=gfff$gene, 
                         start=gfff$start, 
                         end=gfff$end, 
                         strand=gfff$strand, 
                         type=gfff$V3, 
                         name=gfff$gene, 
                         direction= gfff$direction,
                         index = c(1:length(gfff$V1)))
  
  
  return(gene_tab)
  
}

#circuit plot functions
.circuitplot <- function(gene_tab_plot, list) {
  with(list, {
  startloc=1000
  endloc=5000
  gene_tab_plot <- subset(gene_tab_plot, start >= startloc & end <= endloc)
  
  plot <- ggplot() 
  if(length(gene_tab_plot[gene_tab_plot$type=="gene",]$gene) > 0) {
    plot <- plot + 
      geom_cds(data = gene_tab_plot[gene_tab_plot$type=="gene",], 
               aes(xmin=start, xmax=end, y=0, group = gene, label = gene , fill =gene), 
               fontsize=5, 
               h =0.01, 
               show.legend = F) +
      geom_cds(data = gene_tab_plot[gene_tab_plot$type=="gene",], 
               aes(xmin=start, xmax=end, y=-1, group = gene, label = gene , fill =gene), 
               fontsize=0, 
               h =0.01, 
               show.legend = F) +
      geom_text_repel(data = gene_tab_plot[gene_tab_plot$type == "gene",], 
                      aes(x=(start+end)/2, y=-1, label = gene), 
                      check_overlap = F, 
                      show.legend = F, 
                      size=3, 
                      nudge_y = -0.1, 
                      angle = 0, 
                      segment.color = 'transparent')
    
  }
  if(length(gene_tab_plot[gene_tab_plot$type=="terminator",]$gene) > 0) {
    plot <- plot + 
      geom_terminator(data = gene_tab_plot[gene_tab_plot$type == "terminator",], 
                      aes(x=start, y=0, col = gene), 
                      h=0.03, 
                      show.legend = F) + 
      geom_terminator(data = gene_tab_plot[gene_tab_plot$type == "terminator",], 
                      aes(x=start, y=-2, col = gene, lwd = 5), 
                      h=0.03, 
                      show.legend = F)+
      geom_text_repel(data = gene_tab_plot[gene_tab_plot$type == "terminator",], 
                      aes(x=start, y=-2, label = gene), 
                      check_overlap = F, 
                      show.legend = F, 
                      size=2, 
                      nudge_y = -0.1, 
                      angle = 0, 
                      segment.color = 'transparent')
  }
  if(length(gene_tab_plot[gene_tab_plot$type=="promoter",]$gene) > 0) {
    plot <- plot + 
      geom_promoter(data = gene_tab_plot[gene_tab_plot$type == "promoter",], 
                    aes(x=start, y=0, label = gene, col = gene), 
                    fontsize=0, 
                    r = 0.01, 
                    show.legend = F) +
      geom_promoter(data = gene_tab_plot[gene_tab_plot$type == "promoter",], 
                    aes(x=start, y=-0.5, label = gene, col = gene), 
                    fontsize=0, 
                    r = 0.01, 
                    show.legend = F) +
      geom_text_repel(data = gene_tab_plot[gene_tab_plot$type == "promoter",], 
                      aes(x=start, y=-0.5, label = gene), 
                      check_overlap = F, 
                      show.legend = F, 
                      size=2, 
                      nudge_y = -0.1, 
                      angle = 0, 
                      segment.color = 'transparent')
    
  }
  if(length(gene_tab_plot[gene_tab_plot$type=="RBS",]$gene) > 0) {
    plot <- plot + 
      geom_rbs(data = gene_tab_plot[gene_tab_plot$type == "RBS",], 
               aes(x=start, y=0, label = gene, fill = gene),
               fontsize=0, 
               r = 0.005, 
               show.legend = F)   +
      geom_rbs(data = gene_tab_plot[gene_tab_plot$type == "RBS",], 
               aes(x=start, y=-1.5, label = gene, fill = gene),
               fontsize=1, r = 0.01, show.legend = F) +
      geom_text_repel(data = gene_tab_plot[gene_tab_plot$type == "RBS",], 
                      aes(x=start, y=-1.5, label = gene), 
                      check_overlap = F, 
                      show.legend = F, 
                      size=3, 
                      nudge_y = -0.1, 
                      angle = 0, 
                      segment.color = 'transparent')
  }
  
  plot + 
    facet_wrap(~chr) +
    theme_gc()
  
  })
}

.dtedit <- function(input, output, name, thedata, id, 
                       pltfunc = NULL, 
                       input.args = NULL,
                       view.cols = names(thedata),
                       edit.cols = names(thedata),
                       edit.label.cols = edit.cols,
                       input.types,
                       input.choices = NULL,
                       selectize = TRUE,
                       modal.size = 'm',
                       text.width = '100%',
                       textarea.width = '570px',
                       textarea.height = '200px',
                       date.width = '100px',
                       numeric.width = '100px',
                       select.width = '100%',
                       defaultPageLength = 10,
                       title.delete = 'Delete',
                       title.edit = 'Edit',
                       title.add = 'New',
                       label.delete = 'Delete',
                       label.edit = 'Edit',
                       label.add = 'New',
                       label.copy = 'Copy',
                       label.download = 'Download',
                       show.delete = TRUE,
                       show.update = TRUE,
                       show.insert = TRUE,
                       show.copy = TRUE,
                       show.download = TRUE,
                       callback.delete = function(data, row) { },
                       callback.update = function(data, olddata, row) { },
                       callback.insert = function(data, row) { },
                       click.time.threshold = 2, # in seconds
                       datatable.options = list(pageLength=defaultPageLength)
) {
  # Some basic parameter checking
  if(!is.data.frame(thedata) | ncol(thedata) < 1) {
    stop('Must provide a data frame with at least one column.')
  } else if(length(edit.cols) != length(edit.label.cols)) {
    stop('edit.cols and edit.label.cols must be the same length.')
  } else if(!all(view.cols %in% names(thedata))) {
    stop('Not all view.cols are in the data.')
  } else if(!all(edit.cols %in% names(thedata))) {
    stop('Not all edit.cols are in the data.')
  }
  
  if(missing(id)) {
    id <- ''
  } else {
    id <- paste0(id, '-')
  }
  
  
  DataTableName <- paste0(name, 'dt')
  PlotName <- paste0(name, 'plt')
  
  result <- shiny::reactiveValues()
  result$thedata <- thedata
  result$view.cols <- view.cols
  result$edit.cols <- edit.cols
  
  # observeEvent(input[[paste0(name, '_type')]], {
  #   print(input[[paste0(name, '_type')]])
  # })
  
  
  dt.proxy <- DT::dataTableProxy(DataTableName)
  
  selectInputMultiple <- function(...) {
    shiny::selectInput(multiple = TRUE, selectize = selectize, ...)
  }
  
  print(input.args)
  
  valid.input.types <- c('dateInput', 'selectInput', 'numericInput',
                         'textInput', 'textAreaInput', 'passwordInput',
                         'selectInputMultiple')
  inputTypes <- sapply(thedata[,edit.cols], FUN=function(x) {
    switch(class(x),
           list = 'selectInputMultiple',
           character = 'textInput',
           Date = 'dateInput',
           factor = 'selectInput',
           integer = 'numericInput',
           numeric = 'numericInput')
  })
  if(!missing(input.types)) {
    if(!all(names(input.types) %in% edit.cols)) {
      stop('input.types column not a valid editting column: ',
           paste0(names(input.types)[!names(input.types) %in% edit.cols]))
    }
    if(!all(input.types %in% valid.input.types)) {
      stop(paste0('input.types must only contain values of: ',
                  paste0(valid.input.types, collapse = ', ')))
    }
    inputTypes[names(input.types)] <- input.types
  }
  
  # Convert any list columns to characters before displaying
  for(i in 1:ncol(thedata)) {
    if(nrow(thedata) == 0) {
      thedata[,i] <- character()
    } else if(is.list(thedata[,i])) {
      thedata[,i] <- sapply(thedata[,i], FUN = function(x) { paste0(x, collapse = ', ') })
    }
  }
  
  
  output[[DataTableName]] <- DT::renderDataTable({
    print(isolate(thedata))
    thedata[,view.cols]
  }, options = datatable.options, server=TRUE, selection='single', rownames=FALSE)
  
  output[[PlotName]] <- renderPlot({
    inputargslist <- lapply(input.args, function(x)  input[[paste0(name, "_", x[[2]])]])
    
    names(inputargslist) <- foreach::foreach(a=input.args,.combine = "c") %do% a[[2]] 
    
    print(inputargslist)
    print(isolate(result$thedata))
    
    pltfunc(isolate(result$thedata),inputargslist)
    
    
  })
  
  output[[paste0(name, '_download')]] <- downloadHandler(
    filename = function() {
      paste("download", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(data.frame(isolate(thedata)), file, row.names = F)
    }
  )
  
  getFields <- function(typeName, values) {
    fields <- list()
    for(i in seq_along(edit.cols)) {
      if(inputTypes[i] == 'dateInput') {
        value <- ifelse(missing(values),
                        as.character(Sys.Date()),
                        as.character(values[,edit.cols[i]]))
        fields[[i]] <- dateInput(
          inputId = paste0(id, name, typeName, edit.cols[i]),
          label = edit.label.cols[i],
          value = value,
          width = date.width)
      } else if(inputTypes[i] == 'selectInputMultiple') {
        value <- ifelse(missing(values), '', values[,edit.cols[i]])
        if(is.list(value)) {
          value <- value[[1]]
        }
        choices <- ''
        if(!missing(values)) {
          choices <- unique(unlist(values[,edit.cols[i]]))
        }
        if(!is.null(input.choices)) {
          if(edit.cols[i] %in% names(input.choices)) {
            choices <- input.choices[[edit.cols[i]]]
          }
        }
        if(length(choices) == 1) {
          if(choices == '') {
            warning(paste0('No choices available for ', edit.cols[i],
                           '. Specify them using the input.choices parameter'))
          }
        }
        fields[[i]] <- selectInputMultiple(
          inputId = paste0(id, name, typeName, edit.cols[i]),
          label = edit.label.cols[i],
          choices = choices,
          selected = value,
          width = select.width)
      } else if(inputTypes[i] == 'selectInput') {
        value <- ifelse(missing(values), '', as.character(values[,edit.cols[i]]))
        fields[[i]] <- shiny::selectInput(
          inputId = paste0(id, name, typeName, edit.cols[i]),
          label = edit.label.cols[i],
          choices = levels(result$thedata[,edit.cols[i]]),
          selected = value,
          width = select.width)
      } else if(inputTypes[i] == 'numericInput') {
        value <- ifelse(missing(values), 0, values[,edit.cols[i]])
        fields[[i]] <- shiny::numericInput(
          inputId = paste0(id, name, typeName, edit.cols[i]),
          label = edit.label.cols[i],
          value = value,
          width = numeric.width)
      } else if(inputTypes[i] == 'textAreaInput') {
        value <- ifelse(missing(values), '', values[,edit.cols[i]])
        fields[[i]] <- shiny::textAreaInput(
          inputId = paste0(id, name, typeName, edit.cols[i]),
          label = edit.label.cols[i],
          value = value,
          width = textarea.width, height=textarea.height)
      } else if(inputTypes[i] == 'textInput') {
        value <- ifelse(missing(values), '', values[,edit.cols[i]])
        fields[[i]] <- shiny::textInput(
          inputId = paste0(id, name, typeName, edit.cols[i]),
          label = edit.label.cols[i],
          value = value,
          width = text.width)
      } else if(inputTypes[i] == 'passwordInput') {
        value <- ifelse(missing(values), '', values[,edit.cols[i]])
        fields[[i]] <- shiny::passwordInput(
          inputId = paste0(id, name, typeName, edit.cols[i]),
          label = edit.label.cols[i],
          value = value,
          width = text.width)
      } else {
        stop('Invalid input type!')
      }
    }
    return(fields)
  }
  
  output[[paste0(name, '_message')]] <- shiny::renderText('')
  
  updateData <- function(proxy, data, ...) {
    # Convert any list columns to characters before displaying
    for(i in 1:ncol(data)) {
      if(is.list(data[,i])) {
        data[,i] <- sapply(data[,i], FUN = function(x) { paste0(x, collapse = ', ') })
      }
    }
    DT::replaceData(proxy, data, ...)
  }
  
  ##### Insert functions #####################################################
  
  observeEvent(input[[paste0(name, '_add')]], {
    if(!is.null(row)) {
      shiny::showModal(addModal())
    }
  })
  
  insert.click <- NA
  
  observeEvent(input[[paste0(name, '_insert')]], {
    if(!is.na(insert.click)) {
      lastclick <- as.numeric(Sys.time() - insert.click, units = 'secs')
      if(lastclick < click.time.threshold) {
        warning(paste0('Double click detected. Ignoring insert call for ', name, '.'))
        return()
      }
    }
    insert.click <<- Sys.time()
    
    newdata <- result$thedata
    row <- nrow(newdata) + 1
    newdata[row,] <- NA
    
    for(i in edit.cols) {
      if(inputTypes[i] %in% c('selectInputMultiple')) {
        newdata[[i]][row] <- list(input[[paste0(name, '_add_', i)]])
      } else {
        newdata[row,i] <- input[[paste0(name, '_add_', i)]]
      }
    }
    tryCatch({
      callback.data <- callback.insert(data = newdata, row = row)
      if(!is.null(callback.data) & is.data.frame(callback.data)) {
        result$thedata <- callback.data
      } else {
        result$thedata <- newdata
      }
      updateData(dt.proxy,
                 result$thedata[,view.cols],
                 rownames = FALSE)
      
      output[[PlotName]] <- renderPlot({
        inputargslist <- lapply(input.args, function(x)  input[[paste0(name, "_", x[[2]])]])
        
        names(inputargslist) <- foreach::foreach(a=input.args,.combine = "c") %do% a[[2]] 
        
        print(inputargslist)
        print(isolate(result$thedata))
        
        pltfunc(isolate(result$thedata),inputargslist)
        
        
      })
      
      
      shiny::removeModal()
      return(TRUE)
    }, error = function(e) {
      output[[paste0(name, '_message')]] <<- shiny::renderText(geterrmessage())
      return(FALSE)
    })
  })
  
  addModal <- function(row, values) {
    output[[paste0(name, '_message')]] <- shiny::renderText('')
    fields <- getFields('_add_', values)
    shiny::modalDialog(title = title.add,
                       shiny::div(shiny::textOutput(paste0(name, '_message')), style='color:red'),
                       fields,
                       footer = shiny::column(shiny::modalButton('Cancel'),
                                              shiny::actionButton(paste0(id, name, '_insert'), 'Save'),
                                              width=12),
                       size = modal.size
    )
  }
  
  ##### Copy functions #######################################################
  
  observeEvent(input[[paste0(name, '_copy')]], {
    row <- input[[paste0(name, 'dt_rows_selected')]]
    if(!is.null(row)) {
      if(row > 0) {
        shiny::showModal(addModal(values=result$thedata[row,]))
      }
    }
  })
  
  ##### Update functions #####################################################
  
  observeEvent(input[[paste0(name, '_edit')]], {
    row <- input[[paste0(name, 'dt_rows_selected')]]
    if(!is.null(row)) {
      if(row > 0) {
        shiny::showModal(editModal(row))
      }
    }
    print(row)
  })
  
  update.click <- NA
  
  observeEvent(input[[paste0(name, '_update')]], {
    if(!is.na(update.click)) {
      lastclick <- as.numeric(Sys.time() - update.click, units = 'secs')
      if(lastclick < click.time.threshold) {
        warning(paste0('Double click detected. Ignoring update call for ', name, '.'))
        return()
      }
    }
    update.click <- Sys.time()
    
    row <- input[[paste0(name, 'dt_rows_selected')]]
    if(!is.null(row)) {
      if(row > 0) {
        newdata <- result$thedata
        for(i in edit.cols) {
          if(inputTypes[i] %in% c('selectInputMultiple')) {
            newdata[[i]][row] <- list(input[[paste0(name, '_edit_', i)]])
          } else {
            newdata[row,i] <- input[[paste0(name, '_edit_', i)]]
          }
        }
        tryCatch({
          callback.data <- callback.update(data = newdata,
                                           olddata = result$thedata,
                                           row = row)
          if(!is.null(callback.data) & is.data.frame(callback.data)) {
            result$thedata <- callback.data
          } else {
            result$thedata <- newdata
          }
          updateData(dt.proxy,
                     result$thedata[,view.cols],
                     rownames = FALSE)
          
          output[[PlotName]] <- renderPlot({
            inputargslist <- lapply(input.args, function(x)  input[[paste0(name, "_", x[[2]])]])
            
            names(inputargslist) <- foreach::foreach(a=input.args,.combine = "c") %do% a[[2]] 
            
            print(inputargslist)
            print(isolate(result$thedata))
            
            pltfunc(isolate(result$thedata),inputargslist)
            
            
          })
          
          shiny::removeModal()
          return(TRUE)
        }, error = function(e) {
          output[[paste0(name, '_message')]] <<- shiny::renderText(geterrmessage())
          return(FALSE)
        })
      }
    }
    return(FALSE)
  })
  
  editModal <- function(row) {
    output[[paste0(name, '_message')]] <- renderText('')
    fields <- getFields('_edit_', values = result$thedata[row,])
    shiny::modalDialog(title = title.edit,
                       shiny::div(shiny::textOutput(paste0(name, '_message')), style='color:red'),
                       fields,
                       footer = column(shiny::modalButton('Cancel'),
                                       shiny::actionButton(paste0(id, name, '_update'), 'Save'),
                                       width=12),
                       size = modal.size
    )
  }
  
  ##### Delete functions #####################################################
  
  observeEvent(input[[paste0(name, '_remove')]], {
    row <- input[[paste0(name, 'dt_rows_selected')]]
    if(!is.null(row)) {
      if(row > 0) {
        shiny::showModal(deleteModal(row))
      }
    }
    print(row)
  })
  
  observeEvent(input[[paste0(name, '_delete')]], {
    row <- input[[paste0(name, 'dt_rows_selected')]]
    if(!is.null(row)) {
      if(row > 0) {
        newdata <- callback.delete(data = result$thedata, row = row)
        if(!is.null(newdata) & is.data.frame(newdata)) {
          result$thedata <- newdata
        } else {
          result$thedata <- result$thedata[-row,]
        }
        updateData(dt.proxy,
                   result$thedata[,view.cols],
                   rownames = FALSE)
        
        output[[PlotName]] <- renderPlot({
          inputargslist <- lapply(input.args, function(x)  input[[paste0(name, "_", x[[2]])]])
          
          names(inputargslist) <- foreach::foreach(a=input.args,.combine = "c") %do% a[[2]] 
          
          print(inputargslist)
          print(isolate(result$thedata))
          
          pltfunc(isolate(result$thedata),inputargslist)
          
          
        })
        
        # output[[PlotName]] <- renderPlot({
        #   print(input[[paste0(name, '_type')]])
        #    print(isolate(result$thedata))
        #    pltfunc(isolate(result$thedata),tp=input[[paste0(name, '_type')]])
        #  })
        
        
        shiny::removeModal()
        return(TRUE)
      }
    }
    return(FALSE)
  })
  
  deleteModal <- function(row) {
    fields <- list()
    for(i in view.cols) {
      fields[[i]] <- div(paste0(i, ' = ', result$thedata[row,i]))
    }
    
    shiny::modalDialog(title = title.delete,
                       shiny::p('Are you sure you want to delete this record?'),
                       fields,
                       footer = shiny::column(modalButton('Cancel'),
                                              shiny::actionButton(paste0(id, name, '_delete'), 'Delete'),
                                              width=12),
                       size = modal.size
    )
  }
  
  ##### Build the UI for the DataTable and buttons ###########################
  
  plotinput <- function(widget, inputid, label, ..., ns = "shiny") {
    #print(widget)
    #print(paste("formatting shiny widgtes with", widget, "id =", inputid, "label =", label))
    return(call2(widget, 
                 call2("paste0", id, name, inputid), 
                 label,
                 ...,
                 .ns= ns)
    )
  }
  
  defalutwidgets <- list(
    list("actionButton", "_add", label.add),
    list("actionButton", "_edit", label.edit),
    list("actionButton", "_remove", label.delete),
    list("actionButton", "_copy", label.copy),
    list("actionButton", "_download", label.download))
  
  outputlist <- list(expr(DT::dataTableOutput(paste0(id, DataTableName))), 
                     expr(plotOutput(paste0(id, PlotName)))
  )
  
  argslist <- lapply(input.args, function(x) {
    x[[2]] <- paste0("_", x[[2]])
    x
  })
  
  
  argslist <- c(defalutwidgets, argslist)
  
  
  output[[name]] <- shiny::renderUI({
    do.call(div, c(lapply(argslist, function(x) do.call(plotinput, x)), outputlist))
  })
  
  
  
  return(result)
}


sbolcanvas <- function(gff) {
  enssentials = c("chr","gene","start","end","type")
  mat <- as.data.frame(matrix(data = rep(0,10), nrow = 2))
  colnames(mat) <- enssentials
  
  con <- dbConnect(RSQLite::SQLite(), ":memory:")
  
  if(!'gene_tab' %in% dbListTables(con)) {
    genetable <- .load.gff(gff)
    gene_tab <- genetable[, match(enssentials, colnames(genetable))]
    gene_tab$id <- 1:length(gene_tab$chr)
    gene_tab$activation <- " "
    gene_tab$repression <- " "
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
    query <- paste0("INSERT INTO gene_tab (id, chr, gene, start, end, type, activation, repression) VALUES (",
                    "", max(getgenetab()$id) + 1, ", ",
                    "'", as.character(data[row,]$gene), "', ",
                    "'", as.character(data[row,]$gene), "', ",
                    "", as.numeric(data[row,]$start), ", ",
                    "", as.numeric(data[row,]$end), ", ",
                    "'", as.character(data[row,]$type), "', ",
                    "'", as.character(data[row,]$activation), "', ",
                    "'", as.character(data[row,]$repression), "' ",
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
                    "end = '", as.numeric(data[row,]$end), "', ",
                    "activation = '", as.character(data[row,]$activation), "', ",
                    "repression = '", as.character(data[row,]$repression), "', ",
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
  
  ui <- dashboardPage(
    # skin = "midnight",
    header = dashboardHeader(
      title = "Genetic circuit canvas with rvis-sbol",
      dropdownMenu(type = "messages",
                   messageItem(
                     from = "Developer",
                     message = "Qian Li at Zhejiang Lab",
                     time = "2022"
                     
                   ),
                   messageItem(
                     from = "Contact",
                     message = "qianli@zhejianglab.com",
                     icon= icon("list-alt")
                     
                   ),
                   messageItem(
                     from = "Support",
                     message = "https://github.com/lynceuslq/rvis-sbol",
                     icon = icon("life-ring")
                   )
                   
      )
    ),
    body = dashboardBody(
      fluidRow(
        box(
          title = "table of plasmid features", 
          closable = TRUE, 
          width = 12,
          height = "800px",
          solidHeader = FALSE, 
          collapsible = TRUE,
          sidebar = boxSidebar(
            id = "mycardsidebar",
            checkboxGroupInput("coln", 
                               "please select columns to present:", 
                               #                        choices =  colnames(genetab)
            )),
          uiOutput("my_datatable")
        ),
      ),
      fluidRow(
        
        box(
          title = "draw plasmid", 
          closable = TRUE, 
          width = 6,
          height = "800px",
          solidHeader = FALSE, 
          collapsible = TRUE,
          actionButton("draw",label = "Plot genetic circuits"),
          plotOutput("my_plot")
        )
        
      )
    ),
    sidebar = dashboardSidebar(
      minified = FALSE, 
      collapsed = FALSE
      
    ),
    
    controlbar = dashboardControlbar(
      id = "controlbar",
      disable = TRUE
    ),
  )
  
  mywidgets <- list(list("selectInput", "type", "type to plot", c("Wiley", "Wadswort Internation Group", "Wadsworth & Brooks", 
                                                                  "Springer"), "Springer"), list("checkboxInput", "show", "To show plot or not"))
  
  server <- function(input, output) {

    
    gene_tab <- getgenetab()
    
    .dtedit(input, output,
           name = 'my_datatable',
           thedata = gene_tab,
           pltfunc = .circuitplot,
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

    observeEvent(input$update, {
      updateBoxSidebar("mycardsidebar")
    })
    
    v <- reactiveValues(data = { 
      gene_tab
    })
    
    observeEvent(input$my_datatable_cell_edit, {
      #get values
      info = input$my_datatable_cell_edit
      i = as.numeric(info$row)
      j = as.numeric(info$col)
      # k = as.numeric(info$value)
      k = info$value
      
      
      #write values to reactive
      v$data[i,j] <- k
    })
    #render plot
    output$my_plot2 <- renderPlot({
      req(input$draw) #require the input button to be non-0 (ie: don't load the plot when the app first loads)
      isolate(v$data) %>%  #don't react to any changes in the data
        ggplot(aes(x=start,y=end)) +
        geom_point() +
        geom_smooth(method = "lm")
    })
    
    output$my_plot <- renderPlot({
      req(input$draw)
      print(isolate(v$data))
      p <- .circuitplot(isolate(v$data), 1000, 5000)
      p
      
    })
    
  }
  shinyApp(ui = ui, server = server)
}

