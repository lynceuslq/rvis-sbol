#!/bin/Rscript
#author: Qian Li
#contact: 402146079@qq.com


library(ggplot2)
library(grid)
library(gridExtra)
library(ggtree)
library(ggrepel)

circuitplot <- function(file, startloc, endloc) {
  
  gene_tab_plot <- load.gff(file)
  
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