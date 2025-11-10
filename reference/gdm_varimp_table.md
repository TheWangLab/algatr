# Generate a Variable Importance Table for GDM Models

This function generates a table displaying the variable importance for
Generalized Dissimilarity Models (GDM). It can take either a `gdmData`
object or a precomputed variable importance object and outputs a
formatted table.

## Usage

``` r
gdm_varimp_table(
  varimp,
  digits = 2,
  summary_stats = TRUE,
  nPerm = 50,
  geo = TRUE
)
```

## Arguments

- varimp:

  a `gdmData` object or a variable importance object created by running
  [gdm.varImp](https://mfitzpatrick.al.umces.edu/gdm/reference/gdm.varImp.html).

- digits:

  number of digits to include (defaults to 2)

- summary_stats:

  whether to add summary statistics to bottom of table (defaults to
  TRUE).

- nPerm:

  number of permutations to use if `varimp` is a `gdmData` object.
  Default is 50.

- geo:

  whether to include geographic distance in the GDM model. Default is
  TRUE.

## Value

A `gt` table object displaying the variable importance.
