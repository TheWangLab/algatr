# Run LFMM

Run LFMM

## Usage

``` r
lfmm_run(
  gen,
  env,
  K,
  lfmm_method = "ridge",
  p_adj = "fdr",
  sig = 0.05,
  calibrate = "gif"
)
```

## Arguments

- gen:

  genotype dosage matrix (rows = individuals & columns = SNPs) or `vcfR`
  object

- env:

  dataframe with environmental data or a Raster\* type object from which
  environmental values for the coordinates can be extracted

- K:

  number of latent factors (if left as NULL (default), K value selection
  will be conducted)

- lfmm_method:

  lfmm method (either `"ridge"` (default) or `"lasso"`)

- p_adj:

  method to use for p-value correction (defaults to "fdr"); other
  options can be found in
  [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html)

- sig:

  alpha level for determining candidate SNPs (defaults to 0.05)

- calibrate:

  a character string, "gif" or "median+MAD". If the "gif" option is set
  (default), significance values are calibrated by using the genomic
  control method. Genomic control uses a robust estimate of the variance
  of z-scores called "genomic inflation factor". If the "median+MAD"
  option is set, the pvalues are calibrated by computing the median and
  MAD of the zscores. If `NULL`, the pvalues are not calibrated.

## See also

Other LFMM functions:
[`lfmm_df()`](https://thewanglab.github.io/algatr/reference/lfmm_df.md),
[`lfmm_do_everything()`](https://thewanglab.github.io/algatr/reference/lfmm_do_everything.md),
[`lfmm_manhattanplot()`](https://thewanglab.github.io/algatr/reference/lfmm_manhattanplot.md),
[`lfmm_qqplot()`](https://thewanglab.github.io/algatr/reference/lfmm_qqplot.md),
[`lfmm_table()`](https://thewanglab.github.io/algatr/reference/lfmm_table.md),
[`lfmm_test_tidy()`](https://thewanglab.github.io/algatr/reference/lfmm_test_tidy.md),
[`quick_elbow()`](https://thewanglab.github.io/algatr/reference/quick_elbow.md),
[`select_K()`](https://thewanglab.github.io/algatr/reference/select_K.md)
