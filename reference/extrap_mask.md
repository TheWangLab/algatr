# Create raster mask based on coordinates

Creates a raster that can be used to mask areas falling outside the
observation range of a dataset, as defined by coordinates and
corresponding raster values

## Usage

``` r
extrap_mask(coords, envlayers, method = "range", nsd = 2, buffer_width = NULL)

range_mask(coords, envlayers)

sd_mask(coords, envlayers, nsd)

buffer_mask(coords, envlayers, buffer_width = 0.8)

chull_mask(coords, envlayers, buffer_width = NULL)
```

## Arguments

- coords:

  data frame of coordinates (first column should be x and second should
  be y)

- envlayers:

  SpatRaster or Raster\* object with environmental values to base mask
  on

- method:

  method to create mask (can be "range", "sd", "buffer", defaults to
  "range"). See details for more information.

- nsd:

  number of standard deviations to use if using the "sd" method

- buffer_width:

  buffer width to supply to `gBuffer` if using "buffer" method

## Value

SpatRaster where values of 1 indicate areas that fall outside of
observation range

## Details

method can either be:

1.  range - uses `range_mask`, mask all areas with values outside of the
    range of any of the values of the coords

2.  sd - uses `sd_mask`, mask all areas outside the mean +/- stdev\*nsd
    of any of the values of the coords (`nsd` defaults to 2)

3.  buffer - uses `buffer_mask`, mask all areas outside of the
    buffer_width around the coords (`buffer_width` defaults to 0.8)

4.  chull - uses `chull_mask`, mask all areas outside a convex hull of
    the points

## Functions

- `range_mask()`: mask based on range of data

- `sd_mask()`: mask based on mean and standard deviation of data

- `buffer_mask()`: mask based on buffers around points

- `chull_mask()`: mask based on range of data
