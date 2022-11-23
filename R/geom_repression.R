library(ggplot2)
library(grid)
library(gridExtra)

GeomRepression <- ggproto("GeomRepression", Geom,
        required_aes = c("xmin", "xmax", "y"),
        default_aes = aes(h = 0.03, 
                          col = "black",
                          lwd = 5,
                          alpha =1),
        draw_key = draw_key_point,
        draw_panel = function(data, panel_scales, coord) {
                ## Transform the data first
                coords <- coord$transform(data, panel_scales)
                
                ## Let's print out the structure of the 'coords' object
                str(coords)
                
              part1 <-  segmentsGrob(x0 = coords$xmin, 
                                    x1 = coords$xmin, 
                                    y0 = coords$y, 
                                    y1 = coords$y + coords$h, 
                                    gp = gpar(col = alpha(coords$col, coords$alpha), lwd = coords$lwd)) 
              part2 <-  segmentsGrob(x0 = coords$xmin, 
                                     x1 = coords$xmax, 
                                     y0 = coords$y + coords$h, 
                                     y1 = coords$y + coords$h,
                                     gp = gpar(col = alpha(coords$col, coords$alpha), lwd = coords$lwd))
              
              part3 <-  segmentsGrob(x0 = coords$xmax, 
                                    x1 = coords$xmax, 
                                    y0 = coords$y+ coords$h * 2/3, 
                                    y1 = coords$y + coords$h, 
                                    gp = gpar(col = alpha(coords$col, coords$alpha), lwd = coords$lwd)) 
              
              part4 <-  segmentsGrob(x0 = coords$xmax - coords$h / 10, 
                                     x1 = coords$xmax + coords$h / 10, 
                                     y0 = coords$y + coords$h * 2/3, 
                                     y1 = coords$y + coords$h * 2/3,
                                     gp = gpar(col = alpha(coords$col, coords$alpha), lwd = coords$lwd))
              
              gTree(children = gList(part1, part2, part3, part4))  
              
        })

geom_repression <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", na.rm = FALSE, 
                         show.legend = NA, inherit.aes = TRUE, ...) {
        ggplot2::layer(
                geom = GeomRepression, mapping = mapping,  
                data = data, stat = stat, position = position, 
                show.legend = show.legend, inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, ...)
        )
}

