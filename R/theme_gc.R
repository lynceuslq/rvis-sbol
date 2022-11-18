library(ggplot2)
library(grid)
library(gridExtra)

theme_gc <- function () 
{
    ggplot2::theme_grey() + ggplot2::theme(panel.background = ggplot2::element_blank(), 
        panel.grid.major.y = ggplot2::element_line(colour = "grey", 
            size = 1), 
        panel.grid.minor.y = ggplot2::element_blank(), 
        panel.grid.minor.x = ggplot2::element_blank(), 
        panel.grid.major.x = ggplot2::element_blank(), 
        
        axis.ticks.y = ggplot2::element_blank(), 
        axis.title.y=ggplot2::element_blank(),
        axis.text.y=ggplot2::element_blank(),
        
        axis.line.x = ggplot2::element_line(colour = "grey20", 
            size = 0.5), 
        axis.ticks.x = ggplot2::element_line(colour = "grey20", 
            size = 0.5), 
        strip.text = ggplot2::element_blank(), 
        strip.background = ggplot2::element_blank()
        )
}
