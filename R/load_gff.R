gff.proc.func <- function(gfff) {
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
