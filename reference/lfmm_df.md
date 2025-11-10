# Convert LFMM results into a tidy dataframe for downstream processing

Convert LFMM results into a tidy dataframe for downstream processing

## Usage

``` r
lfmm_df(x)
```

## Arguments

- x:

  lfmm_test_result element from
  [`lfmm_run`](https://thewanglab.github.io/algatr/reference/lfmm_run.md)
  results

## Value

tidy dataframe with LFMM results with each SNP, its p-value, association
with env var and other relevant statistics

## See also

Other LFMM functions:
[`lfmm_do_everything()`](https://thewanglab.github.io/algatr/reference/lfmm_do_everything.md),
[`lfmm_manhattanplot()`](https://thewanglab.github.io/algatr/reference/lfmm_manhattanplot.md),
[`lfmm_qqplot()`](https://thewanglab.github.io/algatr/reference/lfmm_qqplot.md),
[`lfmm_run()`](https://thewanglab.github.io/algatr/reference/lfmm_run.md),
[`lfmm_table()`](https://thewanglab.github.io/algatr/reference/lfmm_table.md),
[`lfmm_test_tidy()`](https://thewanglab.github.io/algatr/reference/lfmm_test_tidy.md),
[`quick_elbow()`](https://thewanglab.github.io/algatr/reference/quick_elbow.md),
[`select_K()`](https://thewanglab.github.io/algatr/reference/select_K.md)
