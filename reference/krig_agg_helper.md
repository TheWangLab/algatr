# Helper function for krig_helper; calculates aggregated/disaggregated raster cell size

Helper function for krig_helper; calculates aggregated/disaggregated
raster cell size

## Usage

``` r
krig_agg_helper(to_krig, agg_disagg, agg_spec = "agg")
```

## Arguments

- to_krig:

  object to create grid for kriging, `grd` or `r` from
  [krig_gd](https://rdrr.io/pkg/wingen/man/krig_gd.html)

- agg_disagg:

  aggregation or disaggregation parameter, one of agg_grd, disagg_grd,
  agg_r, or disagg_r from
  [krig_gd](https://rdrr.io/pkg/wingen/man/krig_gd.html)

- agg_spec:

  whether aggregation or disaggregation is performed, options are
  `"disagg"` or `"agg"` (defaults to `"agg"`)

## Value

number of cells contained in final (aggregated or disaggregated) raster
layer

## See also

Other wingen functions:
[`krig_helper()`](https://thewanglab.github.io/algatr/reference/krig_helper.md),
[`wingen_do_everything()`](https://thewanglab.github.io/algatr/reference/wingen_do_everything.md)
