# Helper function for kriging wingen moving window; checks raster resolution and runs kriging

Helper function for kriging wingen moving window; checks raster
resolution and runs kriging

## Usage

``` r
krig_helper(
  map,
  grd = NULL,
  index = 1,
  agg_grd = NULL,
  disagg_grd = NULL,
  agg_r = NULL,
  disagg_r = NULL
)
```

## Arguments

- map:

  RasterLayer or RasterStack produced by `window_gd()`

- grd:

  object to create grid for kriging; inherits from
  [krig_gd](https://rdrr.io/pkg/wingen/man/krig_gd.html) (defaults to
  NULL)

- index:

  integer indices of layers in raster stack to krige; inherits from
  [krig_gd](https://rdrr.io/pkg/wingen/man/krig_gd.html) (defaults to 1)

- agg_grd:

  factor to use for aggregation of `grd`; inherits from
  [krig_gd](https://rdrr.io/pkg/wingen/man/krig_gd.html) (defaults to
  NULL)

- disagg_grd:

  factor to use for disaggregation of `grd`; inherits from
  [krig_gd](https://rdrr.io/pkg/wingen/man/krig_gd.html) (defaults to
  NULL)

- agg_r:

  factor to use for aggregation of `r`; inherits from
  [krig_gd](https://rdrr.io/pkg/wingen/man/krig_gd.html) (defaults to
  NULL)

- disagg_r:

  factor to use for disaggregation of `r`; inherits from
  [krig_gd](https://rdrr.io/pkg/wingen/man/krig_gd.html) (defaults to
  NULL)

## Value

kriged map

## See also

Other wingen functions:
[`krig_agg_helper()`](https://thewanglab.github.io/algatr/reference/krig_agg_helper.md),
[`wingen_do_everything()`](https://thewanglab.github.io/algatr/reference/wingen_do_everything.md)
