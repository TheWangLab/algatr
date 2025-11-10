# TESS function to do everything

TESS function to do everything

## Usage

``` r
tess_do_everything(
  gen,
  coords,
  grid = NULL,
  Kvals = 1:10,
  K_selection = "manual",
  plot_method = "maxQ",
  col_breaks = 20,
  minQ = 0.1,
  tess_method = "projected.ls",
  lambda = 1,
  ploidy = 2,
  correct_kriged_Q = TRUE,
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

- grid:

  SpatRaster for kriging

- Kvals:

  vector of K values to test

- K_selection:

  how to perform K selection ("manual" to enter into console (default)
  or "auto" for automatic selection based on
  [bestK](https://thewanglab.github.io/algatr/reference/bestK.md))

- plot_method:

  method for making rainbow map of kriged layers (options: "maxQ" to
  only plot the max Q value for each cell (default), "allQ" to plot all
  Q values greater than `minQ`, "maxQ_poly" or "allQ_poly" to create the
  plots as previously described, but as polygons for each K instead of
  continuous Q values)

- col_breaks:

  number of breaks for plotting (defaults to 20)

- minQ:

  threshold for minimum Q value for rainbow plotting if `method = "all"`
  is used (defaults to 0.10)

- tess_method:

  the type of TESS method to be run ("projected.ls" for projected least
  squares algorithm (default) or "qp" for quadratic programming
  algorithm)

- lambda:

  numeric value for the spatial regularization parameter. The default
  value lambda = 1 attributes equal weights to the loss function and to
  the penalty function.

- ploidy:

  ploidy of data (defaults to 2)

- correct_kriged_Q:

  whether to correct kriged Q values so values greater than 1 are set to
  1 and values less than 0 are set to 0 (defaults to TRUE)

- quiet:

  whether to operate quietly and suppress the output of tables and
  figures (defaults to FALSE)

## Value

list with all TESS results, final K value, and final kriged raster

## Details

TESS is run using the tess3r package: Caye, K., François, O. (2016).
tess3r: Inference of Spatial Population Genetic Structure. R package
version 1.1.0. See also: Caye, K., Deist, T.M., Martins, H., Michel, O.,
François, O. (2016). TESS3: fast inference of spatial population
structure and genome scans for selection. Mol. Ecol. Res. 16(2):540-548.
https://doi.org/10.1111/1755-0998.12471

Coordinates and rasters should be in a projected (planar) coordinate
system. Therefore, spherical systems (including latitute-longitude
coordinate systems) should be projected prior to use. Transformation can
be performed using
[st_set_crs](https://r-spatial.github.io/sf/reference/st_crs.html) for
coordinates or
[project](https://rspatial.github.io/terra/reference/project.html) for
rasters (see vignette for more details).

## See also

Other TESS functions:
[`allK_plot_helper()`](https://thewanglab.github.io/algatr/reference/allK_plot_helper.md),
[`bestK()`](https://thewanglab.github.io/algatr/reference/bestK.md),
[`geom_tess()`](https://thewanglab.github.io/algatr/reference/geom_tess.md),
[`ggbarplot_helper()`](https://thewanglab.github.io/algatr/reference/ggbarplot_helper.md),
[`pops_helper()`](https://thewanglab.github.io/algatr/reference/pops_helper.md),
[`tess_barplot()`](https://thewanglab.github.io/algatr/reference/tess_barplot.md),
[`tess_col_default()`](https://thewanglab.github.io/algatr/reference/tess_col_default.md),
[`tess_ggbarplot()`](https://thewanglab.github.io/algatr/reference/tess_ggbarplot.md),
[`tess_ggplot()`](https://thewanglab.github.io/algatr/reference/tess_ggplot.md),
[`tess_krig()`](https://thewanglab.github.io/algatr/reference/tess_krig.md),
[`tess_ktest()`](https://thewanglab.github.io/algatr/reference/tess_ktest.md),
[`tess_legend()`](https://thewanglab.github.io/algatr/reference/tess_legend.md),
[`tess_plot_allK()`](https://thewanglab.github.io/algatr/reference/tess_plot_allK.md)
