library(ggplot2)
library(grid)
library(gridExtra)


semicir_x <- function(center=c(0,0), diameter=1, npoints=10){
    tt <- seq(0, pi, length.out=npoints)
    return(center[1] + diameter / 2 * cos(tt))
    
}

semicir_y <- function(center=c(0,0), diameter=1, npoints=10){
    tt <- seq(0, pi, length.out=npoints)
    return(center[2] + diameter / 2 * sin(tt))
    
}


GeomRBS <- ggproto("GeomRBS", Geom,
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
                
              point <-  curveGrob(x1 = coords$x - coords$r, 
                                   x2 = coords$x + coords$r, 
                                   y1 = coords$y,
                                   y2 = coords$y,
                                  curvature=1,
                                #   open = F, 
                                #   angle = 90,
                                   gp = gpar(fill = alpha(coords$fill, coords$alpha)))
              
              p <-  polygonGrob(x = foreach::foreach(a=coords$x, .combine = "c") %do% semicir_x(c(a,0), diameter=2*coords$r, npoints=20),
                                y = foreach::foreach(a=coords$y, .combine = "c") %do% semicir_y(c(0,a), diameter=2*coords$r, npoints=20),
                                gp = gpar(fill = alpha(coords$fill, coords$alpha)),
                                id = foreach::foreach(a=1:length(coords$x), .combine = "c") %do% rep(a, 20)
                            #    id = rep(1:length(coords$x), 10)
                                )
              
              text <- textGrob(label = coords$label, 
                                x = coords$x,
                                y = coords$y + coords$r + coords$h + 0.02,
                               gp=gpar(fontsize=coords$fontsize))
              
              gTree(children = gList(p, text))  
              
        })

geom_rbs <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", na.rm = FALSE, 
                         show.legend = NA, inherit.aes = TRUE, ...) {
        ggplot2::layer(
                geom = GeomRBS, mapping = mapping,  
                data = data, stat = stat, position = position, 
                show.legend = show.legend, inherit.aes = inherit.aes,
                params = list(na.rm = na.rm, ...)
        )
}
