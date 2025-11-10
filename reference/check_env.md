# Check environmental layers for collinearity

Check environmental layers for collinearity

## Usage

``` r
check_env(envlayers, threshold = 0.7)
```

## Arguments

- envlayers:

  SpatRaster or Raster\* object

- threshold:

  the cutoff correlation coefficient for flagging variables as collinear
  (numeric; defaults to 0.7)

## Value

a matrix of correlation coefficients
