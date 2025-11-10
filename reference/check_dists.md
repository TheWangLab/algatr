# Check geographic and environmental distances for collinearity

Check geographic and environmental distances for collinearity

## Usage

``` r
check_dists(envlayers, coords, type = "Euclidean", lyr = NULL, sig = 0.05)
```

## Arguments

- envlayers:

  SpatRaster or Raster\* object

- coords:

  dataframe with x and y sample coordinates

- type:

  the type of geographic distance to be calculated; options are
  "Euclidean" for direct distance, "topographic" for topographic
  distances, and "resistance" for resistance distances

- lyr:

  DEM raster for calculating topographic distances or resistance raster
  for calculating resistance distances

- sig:

  significance threshold for Mantel test

## Value

a list with (1) a dataframe of significantly correlated variables, (2) a
matrix of p-values, (3) a matrix of Mantel's r
