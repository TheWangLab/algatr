# wingen function to do everything (preview and generate moving window maps, krige, and mask)

wingen function to do everything (preview and generate moving window
maps, krige, and mask)

## Usage

``` r
wingen_do_everything(
  gen,
  lyr,
  coords,
  wdim = 3,
  fact = 0,
  sample_count = TRUE,
  min_n = 2,
  preview = FALSE,
  stat = "pi",
  rarify = FALSE,
  kriged = FALSE,
  grd = NULL,
  index = 1,
  agg_grd = NULL,
  disagg_grd = NULL,
  agg_r = NULL,
  disagg_r = NULL,
  masked = FALSE,
  mask = NULL,
  bkg = NULL,
  plot_count = FALSE,
  quiet = FALSE
)
```

## Arguments

- gen:

  Genetic data either as an object of type vcf or a path to a vcf file
  (*note:* order matters! The coordinate and genetic data should be in
  the same order; there are currently no checks for this).

- lyr:

  SpatRaster or RasterLayer to slide the window across (see Details for
  important information about projections). For `method = "resist"` this
  should also be the conductivity layer (see
  [resist_gd](https://rdrr.io/pkg/wingen/man/resist_gd.html)).

- coords:

  Coordinates of samples as sf points, a two-column matrix, or a
  data.frame representing x and y coordinates (see Details for important
  information about projections).

- wdim:

  If `method = "window"`, dimensions (height x width) of window; if only
  one value is provided, a square window is created (defaults to 3 x 3
  window).

- fact:

  Aggregation factor to apply to `lyr` (defaults to 0; *note:*
  increasing this value reduces computational time).

- sample_count:

  Whether to create plot of sample counts for each cell (defaults to
  TRUE).

- min_n:

  Minimum number of samples to use in calculations (any focal cell with
  a window containing less than this number of samples will be assigned
  a value of NA).

- preview:

  whether to produce preview of raster layer, window and focal cell size
  using [preview_gd](https://rdrr.io/pkg/wingen/man/preview_gd.html)
  (default = FALSE)

- stat:

  Genetic diversity statistic(s) to calculate (see Details, defaults to
  `"pi"`). Can be a single statistic or a vector of statistics.

- rarify:

  If rarify = TRUE, rarefaction is performed (defaults to FALSE).

- kriged:

  whether to smooth out mapped values using kriging using
  [krig_gd](https://rdrr.io/pkg/wingen/man/krig_gd.html) (default =
  FALSE)

- grd:

  Object to create grid for kriging; can be a SpatRaster or RasterLayer.
  If undefined, will use `r` to create a grid.

- index:

  Integer indices of layers in raster stack to krige (defaults to 1;
  i.e., the first layer).

- agg_grd:

  Factor to use for aggregation of `grd`, if provided (this will
  decrease the resolution of the final kriged raster; defaults to NULL).

- disagg_grd:

  Factor to use for disaggregation of `grd`, if provided (this will
  increase the resolution of the final kriged raster; defaults to NULL).

- agg_r:

  Factor to use for aggregation of `r`, if provided (this will decrease
  the number of points used in the kriging model; defaults to NULL).

- disagg_r:

  Factor to use for disaggregation of `r`, if provided (this will
  increase the number of points used in the kriging model; defaults to
  NULL).

- masked:

  whether to mask out areas outside region of interest using
  [mask_gd](https://rdrr.io/pkg/wingen/man/mask_gd.html) (default =
  FALSE)

- plot_count:

  if TRUE, whether to visualize sample counts using
  [plot_count](https://rdrr.io/pkg/wingen/man/plot_count.html) (default
  = FALSE)

- quiet:

  whether to operate quietly and suppress the output of tables and
  figures (defaults to FALSE)

## Value

RasterBrick object final raster

## Details

When using wingen, please cite the original citation: Bishop, A.P.,
Chambers, E.A., Wang, I.J. (2023). Generating continuous maps of genetic
diversity using moving windows. Methods Ecol. Evol. doi:
https://doi.org/10.1111/2041-210X.14090 N.B.: Be aware that this
function sets *many* of the wingen function arguments to defaults, which
my result in sub optimal results. We highly advise researchers to run
each wingen function separately for best results.

## See also

Other wingen functions:
[`krig_agg_helper()`](https://thewanglab.github.io/algatr/reference/krig_agg_helper.md),
[`krig_helper()`](https://thewanglab.github.io/algatr/reference/krig_helper.md)
