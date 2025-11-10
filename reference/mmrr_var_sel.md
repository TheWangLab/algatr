# mmrr_var_sel performs MMRR with backward elimination variable selection

mmrr_var_sel performs MMRR with backward elimination variable selection

## Usage

``` r
mmrr_var_sel(Y, X, nperm = 999, stdz = TRUE)
```

## Arguments

- Y:

  is a dependent distance matrix

- X:

  is a list of independent distance matrices (with optional names)

- nperm:

  number of permutations to be used in significance tests (default =
  999)

- stdz:

  if TRUE then matrices will be standardized (default = TRUE)

## See also

Other MMRR functions:
[`MMRR()`](https://thewanglab.github.io/algatr/reference/MMRR.md),
[`mmrr_df()`](https://thewanglab.github.io/algatr/reference/mmrr_df.md),
[`mmrr_do_everything()`](https://thewanglab.github.io/algatr/reference/mmrr_do_everything.md),
[`mmrr_plot()`](https://thewanglab.github.io/algatr/reference/mmrr_plot.md),
[`mmrr_run()`](https://thewanglab.github.io/algatr/reference/mmrr_run.md),
[`mmrr_table()`](https://thewanglab.github.io/algatr/reference/mmrr_table.md),
[`unfold()`](https://thewanglab.github.io/algatr/reference/unfold.md)
