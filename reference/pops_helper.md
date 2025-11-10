# Helper function to get population assignments for best K

Helper function to get population assignments for best K

## Usage

``` r
pops_helper(gen, tess3_obj, K)
```

## Arguments

- gen:

  genotype dosage matrix (rows = individuals & columns = SNPs) or `vcfR`
  object

- tess3_obj:

  list produced by `tess3`

- K:

  K value

## Value

population assignments for each individual based on max Q values

## See also

Other TESS functions:
[`allK_plot_helper()`](https://thewanglab.github.io/algatr/reference/allK_plot_helper.md),
[`bestK()`](https://thewanglab.github.io/algatr/reference/bestK.md),
[`geom_tess()`](https://thewanglab.github.io/algatr/reference/geom_tess.md),
[`ggbarplot_helper()`](https://thewanglab.github.io/algatr/reference/ggbarplot_helper.md),
[`tess_barplot()`](https://thewanglab.github.io/algatr/reference/tess_barplot.md),
[`tess_col_default()`](https://thewanglab.github.io/algatr/reference/tess_col_default.md),
[`tess_do_everything()`](https://thewanglab.github.io/algatr/reference/tess_do_everything.md),
[`tess_ggbarplot()`](https://thewanglab.github.io/algatr/reference/tess_ggbarplot.md),
[`tess_ggplot()`](https://thewanglab.github.io/algatr/reference/tess_ggplot.md),
[`tess_krig()`](https://thewanglab.github.io/algatr/reference/tess_krig.md),
[`tess_ktest()`](https://thewanglab.github.io/algatr/reference/tess_ktest.md),
[`tess_legend()`](https://thewanglab.github.io/algatr/reference/tess_legend.md),
[`tess_plot_allK()`](https://thewanglab.github.io/algatr/reference/tess_plot_allK.md)
