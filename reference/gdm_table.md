# Create `gt` table of GDM results

Create `gt` table of GDM results

## Usage

``` r
gdm_table(gdm_result, digits = 2, summary_stats = TRUE, footnote = TRUE)
```

## Arguments

- gdm_result:

  output of
  [gdm_run](https://thewanglab.github.io/algatr/reference/gdm_run.md) or
  [gdm_do_everything](https://thewanglab.github.io/algatr/reference/gdm_do_everything.md)
  or a GDM model object

- digits:

  number of digits to include (defaults to 2)

## Value

An object of class `gt_tbl`

## See also

Other GDM functions:
[`gdm_coeffs()`](https://thewanglab.github.io/algatr/reference/gdm_coeffs.md),
[`gdm_df()`](https://thewanglab.github.io/algatr/reference/gdm_df.md),
[`gdm_do_everything()`](https://thewanglab.github.io/algatr/reference/gdm_do_everything.md),
[`gdm_format()`](https://thewanglab.github.io/algatr/reference/gdm_format.md),
[`gdm_map()`](https://thewanglab.github.io/algatr/reference/gdm_map.md),
[`gdm_plot_diss()`](https://thewanglab.github.io/algatr/reference/gdm_plot_diss.md),
[`gdm_plot_isplines()`](https://thewanglab.github.io/algatr/reference/gdm_plot_isplines.md),
[`gdm_plot_vars()`](https://thewanglab.github.io/algatr/reference/gdm_plot_vars.md),
[`gdm_run()`](https://thewanglab.github.io/algatr/reference/gdm_run.md),
[`gdm_var_sel()`](https://thewanglab.github.io/algatr/reference/gdm_var_sel.md),
[`scale01()`](https://thewanglab.github.io/algatr/reference/scale01.md)
