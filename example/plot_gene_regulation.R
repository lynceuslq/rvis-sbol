gunew <- read.table(file="plot.repression.tsv", sep="\t", header=T)

gff_file <- read.table(file="induction.cir.gff", header=F, sep="\t")
gene_tab_plot <- subset(gff.proc.func(gff_file))

ggplot() +
#  geom_terminator(data = gene_tab_plot[gene_tab_plot$type=="terminator",], aes(x=start, y=0, group = gene, label = gene , fill =gene), fontsize=5, h =0.01, show.legend = F)
  geom_repression(data = gunew, aes(xmin=donorloc, xmax=receptorloc , y=0, col =gene), fontsize=5, h =0.07, show.legend = F) + theme_gc() + facet_wrap(~chr) +
  geom_cds(data = subset(gene_tab_plot, type=="gene" & chr == "0x41v70"), aes(xmin=start, xmax=end, y=0, group = gene, label = gene , fill =gene), fontsize=5, h =0.01, show.legend = F) +
  geom_promoter(data = subset(gene_tab_plot, type=="promoter" & chr == "0x41v70"), aes(x=start, y=0, label = gene, col = gene), fontsize=4, r = 0.01, show.legend = F) +
  geom_text_repel(data = subset(gene_tab_plot, type=="promoter" & chr == "0x41v70"), aes(x=start, y=0, label = gene), check_overlap = F, show.legend = F, size=2, nudge_y = -0.005, angle = 0,
                  segment.color = 'grey'
              #    segment.color = 'transparent'
                  )
