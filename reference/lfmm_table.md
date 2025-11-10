# Create `gt` table of LFMM results

Create `gt` table of LFMM results

## Usage

``` r
lfmm_table(
  df,
  sig = 0.05,
  sig_only = TRUE,
  top = FALSE,
  order = FALSE,
  var = NULL,
  rows = NULL,
  digits = 2,
  footnotes = TRUE
)
```

## Arguments

- df:

  df element from
  [`lfmm_run`](https://thewanglab.github.io/algatr/reference/lfmm_run.md)
  results

- sig:

  alpha level for determining candidate snps (defaults to 0.5)

- sig_only:

  only include SNPs that exceeded the significance threshold in the
  table (defaults to TRUE)

- top:

  if there are SNPs that are significantly associated with multiple
  environmental variables, only display the top association (i.e.,
  variable with the maximum B value; defaults to FALSE)

- order:

  if TRUE, will order rows by decreasing B value (defaults to FALSE and
  orders rows based on variable)

- var:

  display significant SNPs associated with particular environmental
  variable (defaults to NULL)

- rows:

  number of rows to include in table (defaults to NULL)

- digits:

  number of decimal points to include (defaults to 2)

## Value

table of LFMM results

## See also

Other LFMM functions:
[`lfmm_df()`](https://thewanglab.github.io/algatr/reference/lfmm_df.md),
[`lfmm_do_everything()`](https://thewanglab.github.io/algatr/reference/lfmm_do_everything.md),
[`lfmm_manhattanplot()`](https://thewanglab.github.io/algatr/reference/lfmm_manhattanplot.md),
[`lfmm_qqplot()`](https://thewanglab.github.io/algatr/reference/lfmm_qqplot.md),
[`lfmm_run()`](https://thewanglab.github.io/algatr/reference/lfmm_run.md),
[`lfmm_test_tidy()`](https://thewanglab.github.io/algatr/reference/lfmm_test_tidy.md),
[`quick_elbow()`](https://thewanglab.github.io/algatr/reference/quick_elbow.md),
[`select_K()`](https://thewanglab.github.io/algatr/reference/select_K.md)
