# Plot the relationship between two distance metrics

Plot the relationship between two distance metrics

## Usage

``` r
gen_dist_corr(dist_x, dist_y, metric_name_x, metric_name_y)
```

## Arguments

- dist_x:

  df containing square distance matrix for x axis

- dist_y:

  df containing square distance matrix for y axis

- metric_name_x:

  name of distance metric for x axis; if DPS used, must be `"dps"`

- metric_name_y:

  name of distance metric for y axis; if DPS used, must be `"dps"`

## Value

scatterplot comparing two user-defined genetic distance metrics
