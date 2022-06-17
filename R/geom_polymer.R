library(ggplot2)
library(grid)
library(gridExtra)

GeomPolymer <- ggproto("GeomPolymer", Geom,
        required_aes = c("x", "y", "label"),
        default_aes = aes(r = 0.01, 
                          h = 0.02,
                          fill = "lightblue",
                          fontsize = 5,
                          alpha =1),
        draw_key = draw_key_point,
        draw_panel = function(data, panel_scales, coord) {
                ## Transform the data first
                coords <- coord$transform(data, panel_scales)
                
                ## Let's print out the structure of the 'coords' object
                str(coords)
                
              point <-  circleGrob(r = coords$r, 
                                   x = coords$x, 
                                   y = coords$y + coords$r + coords$h, 
                                   gp = gpar(fill = alpha(coords$fill, coords$alpha)))
              stick <-  segmentsGrob(x0 = coords$x, 
                                     x1 = coords$x, 
                                     y0 = coords$y, 
                                     y1 = coords$y + coords$h)
              text <- textGrob(label = coords$label, 
                                x = coords$x,
                                y = coords$y + coords$r + coords$h + 0.02,
                               gp=gpar(fontsize=coords$fontsize))
              gTree(children = gList(point, stick, text))  
              
        })

geom_polymer <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", na.rm = FALSE, 
                         show.legend = NA, inherit.aes = TRUE, ...) {
        ggplot2::layer(
                geom = GeomPolymer, mapping = mapping,  
                data = data, stat = stat, position = position, 
                show.legend = show.legend, inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, ...)
        )
}
