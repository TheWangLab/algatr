# Remove islands from mapping

Remove islands from mapping

## Usage

``` r
rm_islands(input, shape, min_vertices = 10000)
```

## Arguments

- input:

  raster (SpatRaster or Raster\*) or coordinates (sf or two column
  dataframe) to mask or remove island points from

- shape:

  SpatialPolygons, sf, or sf object to create island mask

- min_vertices:

  minimum number of vertices in polygons to retain (defaults to 10000)

## Value

object (of input type) with islands removed
