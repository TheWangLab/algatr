# Krige admixture values

Krige admixture values

## Usage

``` r
tess_krig(qmat, coords, grid = NULL, correct_kriged_Q = TRUE)
```

## Arguments

- qmat:

  qmatrix

- coords:

  coordinates of samples as sf points, a two-column matrix, or a
  data.frame representing x and y coordinates (see Details for important
  information about projections)

- grid:

  SpatRaster for kriging

- correct_kriged_Q:

  whether to correct kriged Q values so values greater than 1 are set to
  1 and values less than 0 are set to 0 (defaults to TRUE)

## Value

Raster\\ type object of kriged Q values

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
[`tess_ktest()`](https://thewanglab.github.io/algatr/reference/tess_ktest.md),
[`tess_legend()`](https://thewanglab.github.io/algatr/reference/tess_legend.md),
[`tess_plot_allK()`](https://thewanglab.github.io/algatr/reference/tess_plot_allK.md)
