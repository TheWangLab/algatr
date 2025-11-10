# Run GDM and return model object

Run GDM and return model object

## Usage

``` r
gdm_run(
  gendist,
  coords,
  env,
  model = "full",
  sig = 0.05,
  nperm = 50,
  scale_gendist = FALSE,
  geodist_type = "Euclidean",
  distPreds = NULL,
  dist_lyr = NULL
)
```

## Arguments

- gendist:

  matrix of genetic distances (must range between 0 and 1 or set
  scale_gendist = TRUE)

- coords:

  dataframe with x (i.e., longitude) and y (i.e., latitude) coordinates;
  must be in this order

- env:

  dataframe or raster object with environmental values for each
  coordinate; if not provided, it will be calculated based on
  coords/envlayers

- model:

  whether to fit the model with all variables ("full") or to perform
  variable selection to determine the best set of variables ("best");
  defaults to "full"

- sig:

  alpha value for significance threshold (defaults to 0.05); only used
  if model = "best"

- nperm:

  number of permutations to use to calculate variable importance; only
  used if model = "best" (defaults to 50)

- scale_gendist:

  whether to scale genetic distance data from 0 to 1 (defaults to FALSE)

- geodist_type:

  the type of geographic distance to be calculated; options are
  "Euclidean" (default) for direct distance, "topographic" for
  topographic distances, and "resistance" for resistance distances.
  Note: creation and plotting of the GDM raster is only possible for
  "Euclidean" distances

- dist_lyr:

  DEM raster for calculating topographic distances or resistance raster
  for calculating resistance distances

## Value

GDM model

## See also

Other GDM functions:
[`gdm_coeffs()`](https://thewanglab.github.io/algatr/reference/gdm_coeffs.md),
[`gdm_df()`](https://thewanglab.github.io/algatr/reference/gdm_df.md),
[`gdm_do_everything()`](https://thewanglab.github.io/algatr/reference/gdm_do_everything.md),
[`gdm_format()`](https://thewanglab.github.io/algatr/reference/gdm_format.md),
[`gdm_map()`](https://thewanglab.github.io/algatr/reference/gdm_map.md),
[`gdm_plot_diss()`](https://thewanglab.github.io/algatr/reference/gdm_plot_diss.md),
[`gdm_plot_isplines()`](https://thewanglab.github.io/algatr/reference/gdm_plot_isplines.md),
[`gdm_plot_vars()`](https://thewanglab.github.io/algatr/reference/gdm_plot_vars.md),
[`gdm_table()`](https://thewanglab.github.io/algatr/reference/gdm_table.md),
[`gdm_var_sel()`](https://thewanglab.github.io/algatr/reference/gdm_var_sel.md),
[`scale01()`](https://thewanglab.github.io/algatr/reference/scale01.md)
