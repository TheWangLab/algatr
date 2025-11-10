# Helper function to plot cross entropy scores from SNMF

Helper function to plot cross entropy scores from SNMF

## Usage

``` r
plot_crossent(ce_values)
```

## Arguments

- ce_values:

  df with run, K-value, and cross entropy created in
  [snmf_bestK](https://thewanglab.github.io/algatr/reference/snmf_bestK.md)

## Value

ggplots of cross entropy values compared to K-values (and across runs)
