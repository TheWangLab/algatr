# MMRR performs Multiple Matrix Regression with Randomization analysis

MMRR performs Multiple Matrix Regression with Randomization analysis

## Usage

``` r
MMRR(Y, X, nperm = 999, scale = TRUE)
```

## Arguments

- Y:

  is a dependent distance matrix

- X:

  is a list of independent distance matrices (with optional names)

- nperm:

  is the number of permutations to be used in significance tests.
  Default = 999.

- scale:

  if TRUE then matrices will be standardized. Default = TRUE.

## Details

When using MMRR, please cite the original citation: Wang I.J. (2013)
Examining the full effects of landscape heterogeneity on spatial genetic
variation: a multiple matrix regression approach for quantifying
geographic and ecological isolation. Evolution, 67: 3403-3411.

## See also

Other MMRR functions:
[`mmrr_df()`](https://thewanglab.github.io/algatr/reference/mmrr_df.md),
[`mmrr_do_everything()`](https://thewanglab.github.io/algatr/reference/mmrr_do_everything.md),
[`mmrr_plot()`](https://thewanglab.github.io/algatr/reference/mmrr_plot.md),
[`mmrr_run()`](https://thewanglab.github.io/algatr/reference/mmrr_run.md),
[`mmrr_table()`](https://thewanglab.github.io/algatr/reference/mmrr_table.md),
[`mmrr_var_sel()`](https://thewanglab.github.io/algatr/reference/mmrr_var_sel.md),
[`unfold()`](https://thewanglab.github.io/algatr/reference/unfold.md)
