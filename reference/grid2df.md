# Convert x,y nodes to data.frame of polygons

Convert x,y nodes to data.frame of polygons

## Usage

``` r
grid2df(
  spe,
  x = spe@metadata$grid_density$node_x,
  y = spe@metadata$grid_density$node_y,
  reverseY = FALSE,
  ...
)
```

## Arguments

- spe:

  A SpatialExperiment object with grid density calculated

- x:

  vector of x nodes of the polygons

- y:

  vector of y nodes of the polygons

- reverseY:

  Reverse y coordinates. Can be numeric to specify the value to subtract
  y coordinates from (reverseY - y coords).

- ...:

  other elements to be stored as columns of the data.frame. Each one
  must be a vector same length as x.

## Value

data.frame with X, Y, and L2. Points with the same L2 belong to the same
polygons

## Details

Basically grid2sf() but returns a data.frame for plotting with
geom_polygon(), which allows for scale\_\*\_transform(), unlike
geom_sf().

Column names are kept similar to sf::st_coordinates
