# Create a PCA plot for GDM

Create a PCA plot for GDM

## Usage

``` r
gdm_plot_vars(
  pcaSamp,
  pcaRast,
  pcaRastRGB,
  coords,
  x = "PC1",
  y = "PC2",
  scl = 1,
  display_axes = FALSE
)
```

## Arguments

- pcaSamp:

  PCA results from running prcomp()

- pcaRast:

  raster PCA

- pcaRastRGB:

  raster PCA rescaled to RGB

- coords:

  dataframe with x and y coordinates

- x:

  x-axis PC

- y:

  y-axis PC

- scl:

  constant for rescaling variable vectors for plotting

- display_axes:

  whether to display axes

## Value

GDM PCA plot

## See also

Other GDM functions:
[`gdm_coeffs()`](https://thewanglab.github.io/algatr/reference/gdm_coeffs.md),
[`gdm_df()`](https://thewanglab.github.io/algatr/reference/gdm_df.md),
[`gdm_do_everything()`](https://thewanglab.github.io/algatr/reference/gdm_do_everything.md),
[`gdm_format()`](https://thewanglab.github.io/algatr/reference/gdm_format.md),
[`gdm_map()`](https://thewanglab.github.io/algatr/reference/gdm_map.md),
[`gdm_plot_diss()`](https://thewanglab.github.io/algatr/reference/gdm_plot_diss.md),
[`gdm_plot_isplines()`](https://thewanglab.github.io/algatr/reference/gdm_plot_isplines.md),
[`gdm_run()`](https://thewanglab.github.io/algatr/reference/gdm_run.md),
[`gdm_table()`](https://thewanglab.github.io/algatr/reference/gdm_table.md),
[`gdm_var_sel()`](https://thewanglab.github.io/algatr/reference/gdm_var_sel.md),
[`scale01()`](https://thewanglab.github.io/algatr/reference/scale01.md)
