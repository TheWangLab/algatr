# Create `gt` table of RDA results

Create `gt` table of RDA results

## Usage

``` r
rda_table(
  cor_df,
  sig = 0.05,
  sig_only = TRUE,
  top = FALSE,
  order = FALSE,
  var = NULL,
  nrow = NULL,
  digits = 2
)
```

## Arguments

- cor_df:

  dataframe of correlation results output from
  [rda_cor](https://thewanglab.github.io/algatr/reference/rda_cor.md)

- sig:

  if `outlier_method = "p"`, the significance level to use to identify
  SNPs (defaults to 0.05)

- sig_only:

  whether to only include loci with p-values less than `sig` (defaults
  to TRUE)

- top:

  whether to only include only keep the top variable for each snp in the
  table by the strength of the correlation (defaults to FALSE)

- order:

  whether to order by the magnitude of the correlation (defaults to
  FALSE)

- var:

  which variables to include (defaults to including all variables)

- nrow:

  number of rows to display (defaults to displaying all rows)

- digits:

  number of digits to include (defaults to 2)

## Value

An object of class `gt_tbl`

## See also

Other RDA functions:
[`rda_cor()`](https://thewanglab.github.io/algatr/reference/rda_cor.md),
[`rda_do_everything()`](https://thewanglab.github.io/algatr/reference/rda_do_everything.md),
[`rda_getoutliers()`](https://thewanglab.github.io/algatr/reference/rda_getoutliers.md),
[`rda_plot()`](https://thewanglab.github.io/algatr/reference/rda_plot.md),
[`rda_run()`](https://thewanglab.github.io/algatr/reference/rda_run.md),
[`rda_varpart()`](https://thewanglab.github.io/algatr/reference/rda_varpart.md),
[`rda_varpart_table()`](https://thewanglab.github.io/algatr/reference/rda_varpart_table.md)
