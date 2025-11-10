# Create geom of TESS results that can be added to a ggplot object

This function creates a ggplot2 geom object for visualizing TESS plots
based on kriging admixture data.

## Usage

``` r
geom_tess(krig_admix, plot_method = "maxQ", minQ = 0.1)
```

## Arguments

- krig_admix:

  SpatRaster returned by
  [tess_krig](https://thewanglab.github.io/algatr/reference/tess_krig.md)

- plot_method:

  method for making rainbow map of kriged layers (options: "maxQ" to
  only plot the max Q value for each cell (default), "allQ" to plot all
  Q values greater than `minQ`, "maxQ_poly" or "allQ_poly" to create the
  plots as previously described, but as polygons for each K instead of
  continuous Q values)

- minQ:

  threshold for minimum Q-value for rainbow plotting if
  `plot_method = "allQ"` or `plot_method = "allQ_poly"` is used
  (defaults to 0.10)

## Value

A list containing ggplot2 geom objects for plotting.

## See also

Other TESS functions:
[`allK_plot_helper()`](https://thewanglab.github.io/algatr/reference/allK_plot_helper.md),
[`bestK()`](https://thewanglab.github.io/algatr/reference/bestK.md),
[`ggbarplot_helper()`](https://thewanglab.github.io/algatr/reference/ggbarplot_helper.md),
[`pops_helper()`](https://thewanglab.github.io/algatr/reference/pops_helper.md),
[`tess_barplot()`](https://thewanglab.github.io/algatr/reference/tess_barplot.md),
[`tess_col_default()`](https://thewanglab.github.io/algatr/reference/tess_col_default.md),
[`tess_do_everything()`](https://thewanglab.github.io/algatr/reference/tess_do_everything.md),
[`tess_ggbarplot()`](https://thewanglab.github.io/algatr/reference/tess_ggbarplot.md),
[`tess_ggplot()`](https://thewanglab.github.io/algatr/reference/tess_ggplot.md),
[`tess_krig()`](https://thewanglab.github.io/algatr/reference/tess_krig.md),
[`tess_ktest()`](https://thewanglab.github.io/algatr/reference/tess_ktest.md),
[`tess_legend()`](https://thewanglab.github.io/algatr/reference/tess_legend.md),
[`tess_plot_allK()`](https://thewanglab.github.io/algatr/reference/tess_plot_allK.md)
