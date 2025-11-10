# Plot RDA results

Plot RDA results

## Usage

``` r
rda_plot(
  mod,
  rda_snps = NULL,
  pvalues = NULL,
  axes = "all",
  biplot_axes = NULL,
  sig = 0.05,
  manhattan = NULL,
  rdaplot = NULL,
  binwidth = NULL
)
```

## Arguments

- mod:

  model object of class `rda`; if this is all that's provided,
  histograms with loadings will be generated

- rda_snps:

  vector of outlier SNPs (defaults to NULL)

- pvalues:

  if creating a Manhattan plot (i.e., `manhattan = TRUE`), a matrix of
  p-values (defaults to NULL)

- axes:

  which RDA axes to include while plotting (defaults to `"all"`)

- biplot_axes:

  if creating an RDA biplot (i.e., `rdaplot = TRUE`), which pairs of
  axes to plot. Defaults to plotting all pairs of axes possible,
  otherwise can be set to a single pair of axes (e.g., c(1,2)) or a list
  of axes pairs (e.g., list(c(1,2), c(2,3))))

- sig:

  if creating a Manhattan plot, significance threshold for y axis
  (defaults to 0.05)

- manhattan:

  whether to produce Manhattan plot (defaults to `TRUE`)

- rdaplot:

  whether to produce an RDA biplot (defaults to `TRUE`). If only one
  axis is provided, instead of a biplot, a histogram will be created

- binwidth:

  width of bins for histograms (defaults to NULL)

## See also

Other RDA functions:
[`rda_cor()`](https://thewanglab.github.io/algatr/reference/rda_cor.md),
[`rda_do_everything()`](https://thewanglab.github.io/algatr/reference/rda_do_everything.md),
[`rda_getoutliers()`](https://thewanglab.github.io/algatr/reference/rda_getoutliers.md),
[`rda_run()`](https://thewanglab.github.io/algatr/reference/rda_run.md),
[`rda_table()`](https://thewanglab.github.io/algatr/reference/rda_table.md),
[`rda_varpart()`](https://thewanglab.github.io/algatr/reference/rda_varpart.md),
[`rda_varpart_table()`](https://thewanglab.github.io/algatr/reference/rda_varpart_table.md)
