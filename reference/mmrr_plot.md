# Plot MMRR results

Plot MMRR results

## Usage

``` r
mmrr_plot(
  Y = NULL,
  X,
  mod = NULL,
  plot_type = "all",
  stdz = TRUE,
  var_names = NULL
)
```

## Arguments

- Y:

  the dependent variable in the form of a distance matrix

- X:

  a list of independent variables in the form of distance matrices
  (required if `plot_type = "fitted"`, `vars` or `"all"`)

- mod:

  the fitted MMRR model (required if `plot_type = "fitted"` or `"all"`)

- plot_type:

  which plots to produce (options: (1) "vars" to plot single variable
  relationships, (2) "fitted" to plot the fitted relationship, (3) "cov"
  to plot covariances between the predictor variables, (4) "all" to
  produce all plots (default))

- stdz:

  if TRUE then matrices will be standardized (default = TRUE)

- var_names:

  add variable names to plot (defaults to NULL)

## Value

plots of MMRR single variable relationships, the fitted relationship,
and the covariances between predictor variables

## See also

Other MMRR functions:
[`MMRR()`](https://thewanglab.github.io/algatr/reference/MMRR.md),
[`mmrr_df()`](https://thewanglab.github.io/algatr/reference/mmrr_df.md),
[`mmrr_do_everything()`](https://thewanglab.github.io/algatr/reference/mmrr_do_everything.md),
[`mmrr_run()`](https://thewanglab.github.io/algatr/reference/mmrr_run.md),
[`mmrr_table()`](https://thewanglab.github.io/algatr/reference/mmrr_table.md),
[`mmrr_var_sel()`](https://thewanglab.github.io/algatr/reference/mmrr_var_sel.md),
[`unfold()`](https://thewanglab.github.io/algatr/reference/unfold.md)
