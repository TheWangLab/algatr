# ggplot of TESS results

ggplot of TESS results

## Usage

``` r
tess_ggplot(
  krig_admix,
  coords = NULL,
  plot_method = "maxQ",
  ggplot_fill = algatr_col_default("ggplot"),
  minQ = 0.1,
  plot_axes = FALSE,
  rel_widths = c(3, 1),
  list = FALSE
)
```

## Arguments

- krig_admix:

  SpatRaster returned by
  [tess_krig](https://thewanglab.github.io/algatr/reference/tess_krig.md)

- coords:

  dataframe with x and y coordinates for plotting (optional)

- plot_method:

  method for making rainbow map of kriged layers (options: "maxQ" to
  only plot the max Q value for each cell (default), "allQ" to plot all
  Q values greater than `minQ`, "maxQ_poly" or "allQ_poly" to create the
  plots as previously described, but as polygons for each K instead of
  continuous Q values)

- ggplot_fill:

  any ggplot2 scale fill discrete function (default:
  scale_fill_viridis_d, `option = "turbo"`)

- minQ:

  threshold for minimum Q-value for rainbow plotting if
  `plot_method = "allQ"` or `plot_method = "allQ_poly"` is used
  (defaults to 0.10)

- plot_axes:

  whether to plot axes or not (defaults to FALSE)

- rel_widths:

  if `plot_method = "maxQ"` or `plot_method = "allQ"` is used, sets
  relative widths of kriged TESS map and legend (defaults to 3:1), from
  [plot_grid](https://wilkelab.org/cowplot/reference/plot_grid.html)

- list:

  if `plot_method = "maxQ"` or `"allQ"`, whether to output list of
  ggplots (i.e., the base plot and legend; TRUE) or a single plot
  (FALSE; default). Use `list = TRUE` if you would like to add
  additional ggplot2 functions to the plot

## Value

ggplot object of TESS results

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
[`tess_krig()`](https://thewanglab.github.io/algatr/reference/tess_krig.md),
[`tess_ktest()`](https://thewanglab.github.io/algatr/reference/tess_ktest.md),
[`tess_legend()`](https://thewanglab.github.io/algatr/reference/tess_legend.md),
[`tess_plot_allK()`](https://thewanglab.github.io/algatr/reference/tess_plot_allK.md)
