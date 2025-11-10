# MMRR function to do everything

MMRR function to do everything

## Usage

``` r
mmrr_do_everything(
  gendist,
  coords,
  env,
  geo = TRUE,
  model = "full",
  geodist_type = "Euclidean",
  dist_lyr = NULL,
  nperm = 999,
  stdz = TRUE,
  quiet = FALSE,
  plot_type = "all"
)
```

## Arguments

- gendist:

  matrix of genetic distances

- coords:

  dataframe with x and y coordinates

- env:

  dataframe with environmental data or a Raster\* type object from which
  environmental values for the coordinates can be extracted

- geo:

  whether to include geographic, topographic, or resistance distance as
  an independent variable (defaults to TRUE)

- model:

  whether to fit the model with all variables (`"full"`) or to perform
  variable selection to determine the best set of variables (`"best"`);
  defaults to "full"

- geodist_type:

  if `geo = TRUE`, the type of geographic distance to be calculated;
  options are "Euclidean" (default) for direct distance, "topographic"
  for topographic distances, and "resistance" for resistance distances

- dist_lyr:

  if `geodist_type = "topographic"`, DEM raster for calculating
  topographic distances or if `geodist_type = "resistance"`, resistance
  raster for calculating resistance distances

- nperm:

  number of permutations to be used in significance tests (default =
  999)

- stdz:

  if TRUE then matrices will be standardized (default = TRUE)

- quiet:

  whether to operate quietly and suppress the output of tables and
  figures (defaults to FALSE)

- plot_type:

  if `quiet = FALSE`, which plots to produce (options: (1) "vars" to
  plot single variable relationships, (2) "fitted" to plot the fitted
  relationship, (3) "cov" to plot covariances between the predictor
  variables, (4) "all" to produce all plots (default))

## Value

list with final model results and regression coefficients

## Details

The MMRR method is described here: Wang, I.J. (2013). Examining the full
effects of landscape heterogeneity on spatial genetic variation: a
multiple matrix regression approach for quantifying geographic and
ecological isolation. Evolution 67(12):3403-3411. doi:
https://doi.org/10.1111/evo.12134

## See also

Other MMRR functions:
[`MMRR()`](https://thewanglab.github.io/algatr/reference/MMRR.md),
[`mmrr_df()`](https://thewanglab.github.io/algatr/reference/mmrr_df.md),
[`mmrr_plot()`](https://thewanglab.github.io/algatr/reference/mmrr_plot.md),
[`mmrr_run()`](https://thewanglab.github.io/algatr/reference/mmrr_run.md),
[`mmrr_table()`](https://thewanglab.github.io/algatr/reference/mmrr_table.md),
[`mmrr_var_sel()`](https://thewanglab.github.io/algatr/reference/mmrr_var_sel.md),
[`unfold()`](https://thewanglab.github.io/algatr/reference/unfold.md)
