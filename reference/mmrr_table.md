# Create `gt` table of MMRR results

Create `gt` table of MMRR results

## Usage

``` r
mmrr_table(mmrr_results, digits = 2, summary_stats = TRUE)
```

## Arguments

- mmrr_results:

  results from MMRR

- digits:

  the number of decimal places to round to

- summary_stats:

  whether to add summary statistics (R-squared, F-statistic, F p-value)
  to bottom of table (defaults to TRUE)

## Value

an object of class `gt_tbl`

## See also

Other MMRR functions:
[`MMRR()`](https://thewanglab.github.io/algatr/reference/MMRR.md),
[`mmrr_df()`](https://thewanglab.github.io/algatr/reference/mmrr_df.md),
[`mmrr_do_everything()`](https://thewanglab.github.io/algatr/reference/mmrr_do_everything.md),
[`mmrr_plot()`](https://thewanglab.github.io/algatr/reference/mmrr_plot.md),
[`mmrr_run()`](https://thewanglab.github.io/algatr/reference/mmrr_run.md),
[`mmrr_var_sel()`](https://thewanglab.github.io/algatr/reference/mmrr_var_sel.md),
[`unfold()`](https://thewanglab.github.io/algatr/reference/unfold.md)
