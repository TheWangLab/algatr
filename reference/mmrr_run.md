# Run MMRR and return model object

Run MMRR and return model object

## Usage

``` r
mmrr_run(Y, X, nperm = 999, stdz = TRUE, model = "full")
```

## Arguments

- Y:

  dependent distance matrix

- X:

  list of independent distance matrices (with optional names)

- nperm:

  number of permutations to be used in significance tests (default =
  999)

- stdz:

  if TRUE then matrices will be standardized (default = TRUE)

- model:

  whether to fit the model with all variables (`"full"`) or to perform
  variable selection to determine the best set of variables (`"best"`);
  defaults to "full"

## Value

list with final model results and regression coefficients

## See also

Other MMRR functions:
[`MMRR()`](https://thewanglab.github.io/algatr/reference/MMRR.md),
[`mmrr_df()`](https://thewanglab.github.io/algatr/reference/mmrr_df.md),
[`mmrr_do_everything()`](https://thewanglab.github.io/algatr/reference/mmrr_do_everything.md),
[`mmrr_plot()`](https://thewanglab.github.io/algatr/reference/mmrr_plot.md),
[`mmrr_table()`](https://thewanglab.github.io/algatr/reference/mmrr_table.md),
[`mmrr_var_sel()`](https://thewanglab.github.io/algatr/reference/mmrr_var_sel.md),
[`unfold()`](https://thewanglab.github.io/algatr/reference/unfold.md)
