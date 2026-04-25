# Convert x,y nodes to sf polygons

Convert x,y nodes to sf polygons

## Usage

``` r
grid2sf(
  spe,
  x = spe@metadata$grid_density$node_x,
  y = spe@metadata$grid_density$node_y,
  reverseY = FALSE
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

## Value

List of sf polygons

## Details

Default is to generate sf polygons for all grid. For plotting with
geom_sf, use sf::st_as_sfc(grid2sf(spe)) to convert list into Geometry
Set.
