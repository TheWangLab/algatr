# Make map from model

Make map from model

## Usage

``` r
gdm_map(
  gdm_model,
  envlayers,
  coords,
  plot_vars = TRUE,
  scl = 1,
  display_axes = FALSE,
  quiet = FALSE
)
```

## Arguments

- gdm_model:

  GDM model

- envlayers:

  SpatRaster or Raster\* object (LAYER NAMES MUST CORRESPOND WITH GDM
  MODEL)

- coords:

  data frame with x and y coordinates

- plot_vars:

  whether to create PCA plot to help in variable and map interpretation
  (defaults to TRUE)

- scl:

  constant for rescaling variable vectors for plotting (defaults to 1)

- display_axes:

  display PC axes text, labels, and ticks (defaults to FALSE)

- quiet:

  whether to operate quietly and suppress the output of tables and
  figures (defaults to FALSE)

## Value

GDM RGB map

## See also

Other GDM functions:
[`gdm_coeffs()`](https://thewanglab.github.io/algatr/reference/gdm_coeffs.md),
[`gdm_df()`](https://thewanglab.github.io/algatr/reference/gdm_df.md),
[`gdm_do_everything()`](https://thewanglab.github.io/algatr/reference/gdm_do_everything.md),
[`gdm_format()`](https://thewanglab.github.io/algatr/reference/gdm_format.md),
[`gdm_plot_diss()`](https://thewanglab.github.io/algatr/reference/gdm_plot_diss.md),
[`gdm_plot_isplines()`](https://thewanglab.github.io/algatr/reference/gdm_plot_isplines.md),
[`gdm_plot_vars()`](https://thewanglab.github.io/algatr/reference/gdm_plot_vars.md),
[`gdm_run()`](https://thewanglab.github.io/algatr/reference/gdm_run.md),
[`gdm_table()`](https://thewanglab.github.io/algatr/reference/gdm_table.md),
[`gdm_var_sel()`](https://thewanglab.github.io/algatr/reference/gdm_var_sel.md),
[`scale01()`](https://thewanglab.github.io/algatr/reference/scale01.md)
