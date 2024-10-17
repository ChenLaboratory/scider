#' Convert x,y nodes to sf polygons
#' 
#' @param spe A SpatialExperiment object with grid density calculated
#' @param x vector of x nodes of the polygons
#' @param y vector of y nodes of the polygons
#' 
#' @return List of sf polygons
#' 
#' @details
#' Default is to generate sf polygons for all grid.
#' For plotting, use sf::st_as_sfc(grid2sf2(spe)) to convert list into Geometry 
#' Set.
grid2sf <- function(spe,x=spe@metadata$grid_density$node_x,y=spe@metadata$grid_density$node_y) {
  if (is.null(spe@metadata$grid_info)) stop("Missing grid. Compute Density first")
  if (is.null(x) || is.null(y)) stop("Missing x or y")
  if (length(x) != length(x)) stop("x, y must be of equal length")
  
  x = as.numeric(x)
  y = as.numeric(y)
  
  nx = spe@metadata$grid_info$dims[1]
  ny = spe@metadata$grid_info$dims[2]
  
  `if`(spe@metadata$grid_info$grid_type=="hex",
       {# hexagon
         dx = spe@metadata$grid_info$xstep/2
         dy = spe@metadata$grid_info$ystep/3
         offset = c(spe@metadata$grid_info$xlim[1] - dx,
                    spe@metadata$grid_info$ylim[1] - dy * 2)
         
         xc = offset[1] + (0:(nx*2+1)) * dx
         yc = offset[2] + (0:(ny*3+1)) * dy
         
         make_poly = function(col,row) {
           is_right = !(row%%2)
           x_index = 2*col+c(0,1,1,0,-1,-1,0)+is_right
           y_index = 3*row+c(2,1,-1,-2,-1,1,2)
           sf::st_polygon(list(matrix(c(xc[x_index],yc[y_index]),7)))
         }
       },
       {# square
         cellsize = c(spe@metadata$grid_info$xstep,spe@metadata$grid_info$ystep)
         offset = c(spe@metadata$grid_info$xlim[1],spe@metadata$grid_info$ylim[1])
         
         xc = offset[1] + (0:nx) * cellsize[1]
         yc = offset[2] + (0:ny) * cellsize[2]
         
         make_poly = function(col,row) {
           x_index = col + c(0,1,1,0,0)
           y_index = row + c(0,0,1,1,0)
           sf::st_polygon(list(matrix(c(xc[x_index],yc[y_index]), 5)))
         }
       }
  )
  return(lapply(1:length(x), function(ii) {
    make_poly(x[ii],y[ii])
  }))
}
