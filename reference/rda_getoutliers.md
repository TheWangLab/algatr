# Get significant outliers from RDA model

Get significant outliers from RDA model

## Usage

``` r
rda_getoutliers(
  mod,
  naxes = "all",
  outlier_method = "p",
  p_adj = "fdr",
  sig = 0.05,
  z = 3,
  plot = TRUE
)
```

## Arguments

- naxes:

  number of RDA axes to use (defaults to "all" to use all axes), if set
  to "manual" a selection option with a terminal prompt will be given,
  otherwise can be any integer that is less than or equal to the total
  number of axes

- outlier_method:

  method to determine outliers. Can either be "p" to use the p-value
  method from
  [here](https://github.com/Capblancq/RDA-landscape-genomics) or "z" to
  use the z-score based method from
  [here](https://popgen.nescent.org/2018-03-27_RDA_GEA.html)

- p_adj:

  if `outlier_method = "p"`, method to use for p-value correction
  (defaults to "fdr"); other options can be found in
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)

- sig:

  if `outlier_method = "p"`, the significance level to use to identify
  SNPs (defaults to 0.05)

- z:

  if `outlier_method = "z"`, the number of standard deviations to use to
  identify SNPs (defaults to 3)

- plot:

  whether to produce scree plot of RDA axes (defaults to TRUE)

## Value

results from outlier tests. If `outlier_method = "p"`, a list of outlier
SNPs, p-values, and results from rdadapt (see [Capblancq et al.
2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12906)).
If `outlier_method = "z"`, a dataframe with outlier SNP Z-scores for
each axis

## See also

Other RDA functions:
[`rda_cor()`](https://thewanglab.github.io/algatr/reference/rda_cor.md),
[`rda_do_everything()`](https://thewanglab.github.io/algatr/reference/rda_do_everything.md),
[`rda_plot()`](https://thewanglab.github.io/algatr/reference/rda_plot.md),
[`rda_run()`](https://thewanglab.github.io/algatr/reference/rda_run.md),
[`rda_table()`](https://thewanglab.github.io/algatr/reference/rda_table.md),
[`rda_varpart()`](https://thewanglab.github.io/algatr/reference/rda_varpart.md),
[`rda_varpart_table()`](https://thewanglab.github.io/algatr/reference/rda_varpart_table.md)
