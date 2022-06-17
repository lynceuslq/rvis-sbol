library(ggplot2)
library(grid)
library(gridExtra)

GeomCDS <- ggproto("GeomCDS", Geom,
        required_aes = c("xmin", "xmax", "y", "label"),
        default_aes = aes(h = 0.01,  
                          fill = "lightblue",
                          fontsize = 5,
                          alpha =1),
        draw_key = draw_key_point,
        draw_panel = function(data, panel_scales, coord) {
                ## Transform the data first
                coords <- coord$transform(data, panel_scales)
                
                ## Let's print out the structure of the 'coords' object
                str(coords)
                
            
              p <-  polygonGrob(x = c(coords$xmin, 
                                      coords$xmin , 
                                      coords$xmin + ( coords$xmax - coords$xmin ) * 4/5, 
                                      coords$xmax, 
                                      coords$xmin + ( coords$xmax - coords$xmin ) * 4/5),
                                y = c( coords$y - coords$h / 2, 
                                       coords$y + coords$h / 2, 
                                       coords$y + coords$h /2 , 
                                       coords$y , 
                                       coords$y - coords$h / 2),
                                gp = gpar(fill = alpha(coords$fill, coords$alpha)),
                                id = rep(1:length(coords$xmin), 5)
                                )

              text <- textGrob(label = coords$label, 
                                x = (coords$xmin + coords$xmax) / 2,
                                y = coords$y + coords$h + 0.02,
                               gp=gpar(fontsize=coords$fontsize),
                               )
              gTree(children = gList(p, text))  
              
        })

geom_cds <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", na.rm = FALSE, 
                         show.legend = NA, inherit.aes = TRUE, ...) {
        ggplot2::layer(
                geom = GeomCDS, mapping = mapping,  
                data = data, stat = stat, position = position, 
                show.legend = show.legend, inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, ...)
        )
}
