#' Convert x,y nodes to sf polygons
#' 
#' @param spe A SpatialExperiment object with grid density calculated
#' @param x vector of x nodes of the polygons
#' @param y vector of y nodes of the polygons
#' @param reverseY Reverse y coordinates. Can be numeric to specify the value to 
#' subtract y coordinates from (reverseY - y coords).
#' 
#' @return List of sf polygons
#' 
#' @details
#' Default is to generate sf polygons for all grid.
#' For plotting with geom_sf, use sf::st_as_sfc(grid2sf(spe)) to convert list
#' into Geometry Set.
#' @export
grid2sf <- function(spe,
                    x=spe@metadata$grid_density$node_x,
                    y=spe@metadata$grid_density$node_y,
                    reverseY = FALSE) {
  if (is.null(spe@metadata$grid_info)) stop("Missing grid. Compute Density first")
  if (is.null(x) || is.null(y)) stop("Missing x or y")
  if (length(x) != length(x)) stop("x, y must be of equal length")
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  nx <- spe@metadata$grid_info$dims[1]
  ny <- spe@metadata$grid_info$dims[2]
  
  `if`(spe@metadata$grid_info$grid_type=="hex",
       {# hexagon
         dx <- spe@metadata$grid_info$xstep/2
         dy <- spe@metadata$grid_info$ystep/3
         offset <- c(spe@metadata$grid_info$xlim[1] - dx,
                     spe@metadata$grid_info$ylim[1] - dy * 2)
         
         xc <- offset[1] + (0:(nx*2+1)) * dx
         yc <- offset[2] + (0:(ny*3+1)) * dy
         is_right <- rep_len(c(0,1),ny)
         
         make_poly <- function(col,row) {
           x_index <- 2*col+c(0,1,1,0,-1,-1,0)+is_right[row]
           y_index <- 3*row+c(2,1,-1,-2,-1,1,2)
           sf::st_polygon(list(matrix(c(xc[x_index],yc[y_index]),7)))
         }
       },
       {# square
         cellsize <- c(spe@metadata$grid_info$xstep,spe@metadata$grid_info$ystep)
         offset <- c(spe@metadata$grid_info$xlim[1],spe@metadata$grid_info$ylim[1])
         
         xc <- offset[1] + (0:nx) * cellsize[1]
         yc <- offset[2] + (0:ny) * cellsize[2]
         
         make_poly <- function(col,row) {
           x_index <- col + c(0,1,1,0,0)
           y_index <- row + c(0,0,1,1,0)
           sf::st_polygon(list(matrix(c(xc[x_index],yc[y_index]), 5)))
         }
       }
  )
  if (is.numeric(reverseY)) {
    yc <- reverseY - yc
  } else if (reverseY == TRUE) {
    # yc <- sum(range(yc)) - yc
    yc <- range(yc)[2] - yc
  }
  
  return(lapply(seq_along(x), function(ii) {
    make_poly(x[ii],y[ii])
  }))
}

#' Convert x,y nodes to data.frame of polygons
#' 
#' @param spe A SpatialExperiment object with grid density calculated
#' @param x vector of x nodes of the polygons
#' @param y vector of y nodes of the polygons
#' @param reverseY Reverse y coordinates. Can be numeric to specify the value to 
#' subtract y coordinates from (reverseY - y coords).
#' @param ... other elements to be stored as columns of the data.frame. Each one
#' must be a vector same length as x.
#' 
#' @return data.frame with X, Y, and L2. Points with the same L2 belong to the 
#' same polygons
#' 
#' @details
#' Basically grid2sf() but returns a data.frame for plotting with geom_polygon(),
#' which allows for scale_*_transform(), unlike geom_sf().
#'
#' Column names are kept similar to sf::st_coordinates
#' @export
grid2df <- function(spe,
                    x=spe@metadata$grid_density$node_x,
                    y=spe@metadata$grid_density$node_y,
                    reverseY = FALSE,
                    ...) {
  if (is.null(spe@metadata$grid_info)) stop("Missing grid. Compute Density first")
  if (is.null(x) || is.null(y)) stop("Missing x or y")
  if (length(x) != length(x)) stop("x, y must be of equal length")
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  nx <- spe@metadata$grid_info$dims[1]
  ny <- spe@metadata$grid_info$dims[2]
  
  if (spe@metadata$grid_info$grid_type=="hex") {
    dx <- spe@metadata$grid_info$xstep/2
    dy <- spe@metadata$grid_info$ystep/3
    offset <- c(spe@metadata$grid_info$xlim[1] - dx,
                spe@metadata$grid_info$ylim[1] - dy * 2)
    
    xc <- offset[1] + (0:(nx*2+1)) * dx
    yc <- offset[2] + (0:(ny*3+1)) * dy
    
    n_vertices = 7
    x_pattern = c(0,1,1,0,-1,-1,0)
    y_pattern = c(2,1,-1,-2,-1,1,2)
    
    x <- 2*x + !y%%2
    y <- 3*y
  } else { # square
    dx <- c(spe@metadata$grid_info$xstep,spe@metadata$grid_info$ystep)
    offset <- c(spe@metadata$grid_info$xlim[1],spe@metadata$grid_info$ylim[1])
    
    xc <- offset[1] + (0:nx) * dx[1]
    yc <- offset[2] + (0:ny) * dx[2]
    
    n_vertices = 5
    x_pattern = c(0,1,1,0,0)
    y_pattern = c(0,0,1,1,0)
  }
  x <- xc[rep(x, each = n_vertices) + x_pattern]
  y <- yc[rep(y, each = n_vertices) + y_pattern]
  L2 <- rep(seq_len(length(x)/n_vertices),each=n_vertices)
  res <- data.frame(X=x,Y=y,L2=L2)
  
  lst <- lapply(list(...),function(ii) rep(ii,each=n_vertices))
  if (length(lst)) res <- cbind(res,lst)
  return(res)
}
