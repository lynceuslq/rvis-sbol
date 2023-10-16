#!/bin/Rscript
#author: Qian Li
#contact: 402146079@qq.com
#An interactive app to load and plot genetic circuits based on SBOL3 standards

library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(ggpubr)
library(ggtree)
library(ggrepel)
library(plasmapR)
library(dplyr)

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
.circuitplot <- function(gene_tab_plot, startloc, endloc) {
  
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
  
  plot <- plot + 
    facet_wrap(~chr) +
    theme_gc()
  
  return(plot)
}

sbolcanvas <- function(gff) {
 # stopifnot(file.exists(gff))
 # genetab <- .load.gff(gff)
  
  
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
          width = 6,
          height = "800px",
          solidHeader = FALSE, 
          collapsible = TRUE,
          sidebar = boxSidebar(
            id = "mycardsidebar",
            checkboxGroupInput("coln", 
                                "please select columns to present:", 
        #                        choices =  colnames(genetab)
          )),
        div(style='height:500px;overflow-y: scroll;',
            uiOutput("my_datatable"))
        ),
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
  
  
  enssentials = c("chr","gene","start","end","type")
  mat <- as.data.frame(matrix(data = rep(0,10), nrow = 2))
  colnames(mat) <- enssentials
  server <- function(input, output) {
   # v <- reactiveValues(data = { 
   #   data.frame(x = numeric(0),y = numeric(0)) %>% 
  #      add_row(x = rep(0,10),y = rep(0,10))
   # })

    if(file.exists(gff)) {
      genetab <- .load.gff(gff)
      v <- reactiveValues(data = genetab[, match(enssentials, colnames(genetab))])
    }else{
      #initialize a blank dataframe
       v <- reactiveValues(data = { 
         mat
       })
    }
    output$my_datatable <- renderUI({
      output$tabtmp <- renderDT({
        v$show <- v$data
        DT::datatable(v$show, editable = TRUE)
      })
      DTOutput('tabtmp', height ='500px') 
    })
   # v <- reactiveValues(data = genetab[, match(enssentials, colnames(genetab))])
    #output the datatable based on the dataframe (and make it editable)
    
    #when there is any edit to a cell, write that edit to the initial dataframe
    #check to make sure it's positive, if not convert
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
    
    observeEvent(input$update, {
      updateBoxSidebar("mycardsidebar")
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

