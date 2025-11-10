# Create TESS barplot using ggplot2

Create TESS barplot using ggplot2

## Usage

``` r
tess_ggbarplot(
  qmat,
  ggplot_fill = algatr_col_default("ggplot"),
  sort_by_Q = TRUE,
  legend = TRUE
)
```

## Arguments

- qmat:

  Q matrix

- ggplot_fill:

  any ggplot2 scale fill discrete function (default:
  scale_fill_viridis_d, `option = "turbo"`)

- sort_by_Q:

  whether to sort bars by Q value (equivalent to barplot sort.by.Q)

- legend:

  whether to display legend (defaults to TRUE)

## Value

ggplot object of TESS results as a barplot

## See also

Other TESS functions:
[`allK_plot_helper()`](https://thewanglab.github.io/algatr/reference/allK_plot_helper.md),
[`bestK()`](https://thewanglab.github.io/algatr/reference/bestK.md),
[`geom_tess()`](https://thewanglab.github.io/algatr/reference/geom_tess.md),
[`ggbarplot_helper()`](https://thewanglab.github.io/algatr/reference/ggbarplot_helper.md),
[`pops_helper()`](https://thewanglab.github.io/algatr/reference/pops_helper.md),
[`tess_barplot()`](https://thewanglab.github.io/algatr/reference/tess_barplot.md),
[`tess_col_default()`](https://thewanglab.github.io/algatr/reference/tess_col_default.md),
[`tess_do_everything()`](https://thewanglab.github.io/algatr/reference/tess_do_everything.md),
[`tess_ggplot()`](https://thewanglab.github.io/algatr/reference/tess_ggplot.md),
[`tess_krig()`](https://thewanglab.github.io/algatr/reference/tess_krig.md),
[`tess_ktest()`](https://thewanglab.github.io/algatr/reference/tess_ktest.md),
[`tess_legend()`](https://thewanglab.github.io/algatr/reference/tess_legend.md),
[`tess_plot_allK()`](https://thewanglab.github.io/algatr/reference/tess_plot_allK.md)
