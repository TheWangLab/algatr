# Plot all kriged Q values for each K

Plot all kriged Q values for each K

## Usage

``` r
tess_plot_allK(
  krig_admix,
  coords = NULL,
  col_pal = algatr_col_default("base"),
  col_breaks = 20,
  ...
)
```

## Arguments

- krig_admix:

  RasterStack returned by
  [tess_krig](https://thewanglab.github.io/algatr/reference/tess_krig.md)

- coords:

  dataframe with x and y coordinates for plotting (optional)

- col_breaks:

  number of breaks for plotting (defaults to 20)

- ...:

  Graphical parameters. Any argument that can be passed to image.plot
  and to base plot

## See also

Other TESS functions:
[`allK_plot_helper()`](https://thewanglab.github.io/algatr/reference/allK_plot_helper.md),
[`bestK()`](https://thewanglab.github.io/algatr/reference/bestK.md),
[`geom_tess()`](https://thewanglab.github.io/algatr/reference/geom_tess.md),
[`ggbarplot_helper()`](https://thewanglab.github.io/algatr/reference/ggbarplot_helper.md),
[`pops_helper()`](https://thewanglab.github.io/algatr/reference/pops_helper.md),
[`tess_barplot()`](https://thewanglab.github.io/algatr/reference/tess_barplot.md),
[`tess_col_default()`](https://thewanglab.github.io/algatr/reference/tess_col_default.md),
[`tess_do_everything()`](https://thewanglab.github.io/algatr/reference/tess_do_everything.md),
[`tess_ggbarplot()`](https://thewanglab.github.io/algatr/reference/tess_ggbarplot.md),
[`tess_ggplot()`](https://thewanglab.github.io/algatr/reference/tess_ggplot.md),
[`tess_krig()`](https://thewanglab.github.io/algatr/reference/tess_krig.md),
[`tess_ktest()`](https://thewanglab.github.io/algatr/reference/tess_ktest.md),
[`tess_legend()`](https://thewanglab.github.io/algatr/reference/tess_legend.md)
