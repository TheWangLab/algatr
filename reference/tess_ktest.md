# Test multiple K values

Test multiple K values

## Usage

``` r
tess_ktest(
  gen,
  coords,
  Kvals = 1:10,
  grid = NULL,
  tess_method = "projected.ls",
  lambda = 1,
  K_selection = "manual",
  ploidy = 2,
  quiet = FALSE
)
```

## Arguments

- gen:

  genotype dosage matrix (rows = individuals & columns = snps) or `vcfR`
  object

- coords:

  coordinates of samples as sf points, a two-column matrix, or a
  data.frame representing x and y coordinates (see Details for important
  information about projections)

- Kvals:

  vector of K values to test

- grid:

  SpatRaster for kriging

- tess_method:

  the type of TESS method to be run ("projected.ls" for projected least
  squares algorithm (default) or "qp" for quadratic programming
  algorithm)

- lambda:

  numeric value for the spatial regularization parameter. The default
  value lambda = 1 attributes equal weights to the loss function and to
  the penalty function.

- K_selection:

  how to perform K selection ("manual" to enter into console (default)
  or "auto" for automatic selection based on
  [bestK](https://thewanglab.github.io/algatr/reference/bestK.md))

- ploidy:

  ploidy of data (defaults to 2)

- quiet:

  whether to operate quietly and suppress the output of tables and
  figures (defaults to FALSE)

## Value

list with results from testing different K-values

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
[`tess_legend()`](https://thewanglab.github.io/algatr/reference/tess_legend.md),
[`tess_plot_allK()`](https://thewanglab.github.io/algatr/reference/tess_plot_allK.md)
