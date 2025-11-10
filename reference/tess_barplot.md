# Create TESS barplot

Based on code from:
https://github.com/bcm-uga/TESS3_encho_sen/blob/master/R/plotQ.R

## Usage

``` r
tess_barplot(
  qmat,
  col_pal = algatr_col_default("base"),
  sort_by_Q = TRUE,
  legend = TRUE,
  legend_position = "bottomright",
  border = NA,
  space = 0,
  ...
)
```

## Arguments

- qmat:

  Q matrix

- sort_by_Q:

  whether to sort bars by Q value (equivalent to barplot sort.by.Q)

- legend:

  whether to display legend (defaults to TRUE)

- legend_position:

  the x and y coordinates or keyword to determine legend position
  (defaults to bottom right)

- border:

  the color to be used for the border of the bars. Use `border = NA` to
  omit borders. If there are shading lines, `border = TRUE` means use
  the same colour for the border as for the shading lines.

- space:

  the amount of space (as a fraction of the average bar width) left
  before each bar. May be given as a single number or one number per
  bar. If `height` is a matrix and `beside` is `TRUE`, `space` may be
  specified by two numbers, where the first is the space between bars in
  the same group, and the second the space between the groups. If not
  given explicitly, it defaults to `c(0,1)` if `height` is a matrix and
  `beside` is `TRUE`, and to 0.2 otherwise.

- ...:

  other parameters of the function
  [`barplot.default`](https://rdrr.io/r/graphics/barplot.html).

## Value

STRUCTURE-style bar plot of TESS results

## See also

Other TESS functions:
[`allK_plot_helper()`](https://thewanglab.github.io/algatr/reference/allK_plot_helper.md),
[`bestK()`](https://thewanglab.github.io/algatr/reference/bestK.md),
[`geom_tess()`](https://thewanglab.github.io/algatr/reference/geom_tess.md),
[`ggbarplot_helper()`](https://thewanglab.github.io/algatr/reference/ggbarplot_helper.md),
[`pops_helper()`](https://thewanglab.github.io/algatr/reference/pops_helper.md),
[`tess_col_default()`](https://thewanglab.github.io/algatr/reference/tess_col_default.md),
[`tess_do_everything()`](https://thewanglab.github.io/algatr/reference/tess_do_everything.md),
[`tess_ggbarplot()`](https://thewanglab.github.io/algatr/reference/tess_ggbarplot.md),
[`tess_ggplot()`](https://thewanglab.github.io/algatr/reference/tess_ggplot.md),
[`tess_krig()`](https://thewanglab.github.io/algatr/reference/tess_krig.md),
[`tess_ktest()`](https://thewanglab.github.io/algatr/reference/tess_ktest.md),
[`tess_legend()`](https://thewanglab.github.io/algatr/reference/tess_legend.md),
[`tess_plot_allK()`](https://thewanglab.github.io/algatr/reference/tess_plot_allK.md)
