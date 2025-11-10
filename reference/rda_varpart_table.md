# Create `gt` table with RDA variance partitioning results

Create `gt` table with RDA variance partitioning results

## Usage

``` r
rda_varpart_table(df, digits = 2, call_col = FALSE)
```

## Arguments

- df:

  dataframe of variance partitioning results output from
  [rda_varpart](https://thewanglab.github.io/algatr/reference/rda_varpart.md)

- digits:

  number of digits to include (defaults to 2)

- call_col:

  whether to include column with RDA call (defaults to FALSE)

## Value

object of class `gt` with RDA variance partitioning results

## See also

Other RDA functions:
[`rda_cor()`](https://thewanglab.github.io/algatr/reference/rda_cor.md),
[`rda_do_everything()`](https://thewanglab.github.io/algatr/reference/rda_do_everything.md),
[`rda_getoutliers()`](https://thewanglab.github.io/algatr/reference/rda_getoutliers.md),
[`rda_plot()`](https://thewanglab.github.io/algatr/reference/rda_plot.md),
[`rda_run()`](https://thewanglab.github.io/algatr/reference/rda_run.md),
[`rda_table()`](https://thewanglab.github.io/algatr/reference/rda_table.md),
[`rda_varpart()`](https://thewanglab.github.io/algatr/reference/rda_varpart.md)
