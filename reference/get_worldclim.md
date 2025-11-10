# Download and merge WorldClim data for study area

Download and merge WorldClim data for study area

## Usage

``` r
get_worldclim(coords, res = 0.5, buff = 0.01, save_output = FALSE)
```

## Arguments

- coords:

  Dataframe with x and y sample coordinates.

- res:

  The resolution of WorldClim data to download; options are 0.5, 2.5, 5,
  and 10 arc-minutes (default = 0.5).

- buff:

  A buffer area around sample points for cropping the data layers,
  expressed as a proportion of the spatial extent for the coordinates
  (default = 0.01).

- save_output:

  Whether to save downloaded worldclim data in a tmp folder in the
  working directory (default = FALSE).

## Value

A SpatRaster of WorldClim layers.

## Details

If res = 0.5 then the individual WorldClim tiles that cover the sample
coordinates are downloaded and merged. If res \> 2.5 then global layers
are downloaded. The buffer area maintains a large extent for the final
cropped data layers around the sample coordinates. e.g. buff = 0.01
creates a 1% buffer area around the coordinates.
