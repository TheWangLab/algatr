# Plot I-splines for each variable

Plot I-splines for each variable

## Usage

``` r
gdm_plot_isplines(gdm_model, scales = "free", nrow = NULL, ncol = NULL)
```

## Arguments

- gdm_model:

  GDM model

- scales:

  Whether scales should be free ("free"; default), free in one dimension
  ("free_x", "free_y") or fixed ("fixed"). We recommend setting this to
  "free_x" to allow the x-axis to vary while keeping the y-axis fixed
  across all plots such that relative importance can be visualized.

- nrow:

  Number of rows

- ncol:

  Number of cols

## Value

plot for each I-spline

## See also

Other GDM functions:
[`gdm_coeffs()`](https://thewanglab.github.io/algatr/reference/gdm_coeffs.md),
[`gdm_df()`](https://thewanglab.github.io/algatr/reference/gdm_df.md),
[`gdm_do_everything()`](https://thewanglab.github.io/algatr/reference/gdm_do_everything.md),
[`gdm_format()`](https://thewanglab.github.io/algatr/reference/gdm_format.md),
[`gdm_map()`](https://thewanglab.github.io/algatr/reference/gdm_map.md),
[`gdm_plot_diss()`](https://thewanglab.github.io/algatr/reference/gdm_plot_diss.md),
[`gdm_plot_vars()`](https://thewanglab.github.io/algatr/reference/gdm_plot_vars.md),
[`gdm_run()`](https://thewanglab.github.io/algatr/reference/gdm_run.md),
[`gdm_table()`](https://thewanglab.github.io/algatr/reference/gdm_table.md),
[`gdm_var_sel()`](https://thewanglab.github.io/algatr/reference/gdm_var_sel.md),
[`scale01()`](https://thewanglab.github.io/algatr/reference/scale01.md)
