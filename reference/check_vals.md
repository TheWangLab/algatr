# Check extracted values for collinearity

Check extracted values for collinearity

## Usage

``` r
check_vals(envlayers, coords, threshold = 0.7)
```

## Arguments

- envlayers:

  SpatRaster or Raster\* object

- coords:

  dataframe with x and y sample coordinates

- threshold:

  the cutoff correlation coefficient for flagging variables as collinear
  (numeric)

## Value

a matrix of correlation coefficients
