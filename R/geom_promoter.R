library(ggplot2)
library(grid)
library(gridExtra)

GeomPromoter <- ggproto("GeomPromoter", Geom,
        required_aes = c("x", "y"),
        default_aes = aes(h = 0.03, 
                          col = "yellow",
                          alpha = 1),
        draw_key = draw_key_point,
        draw_panel = function(data, panel_scales, coord) {
                ## Transform the data first
                coords <- coord$transform(data, panel_scales)
                
                ## Let's print out the structure of the 'coords' object
                str(coords)
                
              part1 <- segmentsGrob(x0 = coords$x, 
                                    x1 = coords$x, 
                                    y0 = coords$y, 
                                    y1 = coords$y + coords$h, 
                                    gp = gpar(col = alpha(coords$col, coords$alpha)))  
              
              part2 <- segmentsGrob(x0 = coords$x, 
                                    x1 = coords$x + coords$h * 2/3 , 
                                    y0 = coords$y + coords$h, 
                                    y1 = coords$y + coords$h, 
                                    arrow = arrow(angle = 30, length = unit(0.05, "inches"), ends = "last", type = "open") , 
                                    gp = gpar(col = alpha(coords$col, coords$alpha))
              )
              gTree(children = gList(part1, part2))
        })


geom_promoter <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", na.rm = FALSE, 
                         show.legend = NA, inherit.aes = TRUE, ...) {
        ggplot2::layer(
                geom = GeomPromoter, mapping = mapping,  
                data = data, stat = stat, position = position, 
                show.legend = show.legend, inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, ...)
        )
}
