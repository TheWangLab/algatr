# Partial RDA variance partitioning

Partial RDA variance partitioning

## Usage

``` r
rda_varpart(gen, env, coords, Pin, R2permutations, R2scope, nPC)
```

## Arguments

- gen:

  genotype dosage matrix (rows = individuals & columns = SNPs) or `vcfR`
  object

- env:

  dataframe with environmental data or a Raster\* type object from which
  environmental values for the coordinates can be extracted

- coords:

  dataframe with coordinates (only needed if correctGEO = TRUE) or if
  env is a Raster\* from which values should be extracted

- Pin:

  if `model = "best"`, limits of permutation P-values for adding (`Pin`)
  a term to the model, or dropping (`Pout`) from the model. Term is
  added if` P <= Pin`, and removed if `P > Pout` (see
  [ordiR2step](https://vegandevs.github.io/vegan/reference/ordistep.html))
  (defaults to 0.05)

- R2permutations:

  if `model = "best"`, number of permutations used in the estimation of
  adjusted R2 for cca using RsquareAdj (see
  [ordiR2step](https://vegandevs.github.io/vegan/reference/ordistep.html))
  (defaults to 1000)

- R2scope:

  if `model = "best"` and set to TRUE (default), use adjusted R2 as the
  stopping criterion: only models with lower adjusted R2 than scope are
  accepted (see
  [ordiR2step](https://vegandevs.github.io/vegan/reference/ordistep.html))

- nPC:

  number of PCs to use if correctPC = TRUE (defaults to 3); if set to
  "manual" a selection option with a terminal prompt will be provided

## Value

df with relevant statistics from variance partitioning analysis

## See also

Other RDA functions:
[`rda_cor()`](https://thewanglab.github.io/algatr/reference/rda_cor.md),
[`rda_do_everything()`](https://thewanglab.github.io/algatr/reference/rda_do_everything.md),
[`rda_getoutliers()`](https://thewanglab.github.io/algatr/reference/rda_getoutliers.md),
[`rda_plot()`](https://thewanglab.github.io/algatr/reference/rda_plot.md),
[`rda_run()`](https://thewanglab.github.io/algatr/reference/rda_run.md),
[`rda_table()`](https://thewanglab.github.io/algatr/reference/rda_table.md),
[`rda_varpart_table()`](https://thewanglab.github.io/algatr/reference/rda_varpart_table.md)
