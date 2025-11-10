# Quickly choose an elbow for a PC

At variance below 5% per component, choose the largest % drop. Designed
for variance percentages, but will also work given a full set of
Evalues. Quickly estimate the 'elbow' of a scree plot (PCA)

## Usage

``` r
quick_elbow(varpc, low = 0.08, max.pc = 0.9)
```

## Arguments

- varpc:

  numeric, vector of eigenvalues, or 'percentage of variance' explained
  by datapoints for each principal component. If only using a partial
  set of components, should first pass to `estimate.eig.vpcs()` to
  estimate any missing eigenvalues

- low:

  numeric (between zero and one); the threshold that defines whether a
  principal component explains 'much' of the variance

- max.pc:

  maximum percentage of the variance to capture before the elbow
  (cumulative sum to PC 'n')

## Value

the number of principal components to keep, prior to the determined
elbow cutoff

## Details

This function uses a rough algorithm to estimate a sensible 'elbow' to
choose for a PCA screeplot of eigenvalues. The function looks at an
initial arbitrarily 'low' level of variance and looks for the first
eigenvalue lower than this. If the very first eigenvalue is actually
lower than this (i.e, when the PCs are not very explanatory) then this
'low' value is iteratively halved until this is no longer the case.
After starting below this arbitrary threshold the drop in variance
explained by each pair of consecutive PCs is standardized by dividing
over the larger of the pair. The largest percentage drop in the series
below 'low' % is selected as the 'elbow'.

## See also

`estimate.eig.vpcs`

Other LFMM functions:
[`lfmm_df()`](https://thewanglab.github.io/algatr/reference/lfmm_df.md),
[`lfmm_do_everything()`](https://thewanglab.github.io/algatr/reference/lfmm_do_everything.md),
[`lfmm_manhattanplot()`](https://thewanglab.github.io/algatr/reference/lfmm_manhattanplot.md),
[`lfmm_qqplot()`](https://thewanglab.github.io/algatr/reference/lfmm_qqplot.md),
[`lfmm_run()`](https://thewanglab.github.io/algatr/reference/lfmm_run.md),
[`lfmm_table()`](https://thewanglab.github.io/algatr/reference/lfmm_table.md),
[`lfmm_test_tidy()`](https://thewanglab.github.io/algatr/reference/lfmm_test_tidy.md),
[`select_K()`](https://thewanglab.github.io/algatr/reference/select_K.md)

## Author

Nicholas Cooper

## Examples

``` r
# correlated data
mat <- sim.cor(100, 50)
#> Error in sim.cor(100, 50): could not find function "sim.cor"
result <- princomp(mat)
#> Error: object 'mat' not found
eig <- result$sdev^2
#> Error: object 'result' not found
elb.a <- quick_elbow(eig)
#> Error: object 'eig' not found
pca.scree.plot(eig, elbow = elb.a, M = mat)
#> Error in pca.scree.plot(eig, elbow = elb.a, M = mat): could not find function "pca.scree.plot"
elb.b <- quick_elbow(eig, low = .05) # decrease 'low' to select more components
#> Error: object 'eig' not found
pca.scree.plot(eig, elbow = elb.b, M = mat)
#> Error in pca.scree.plot(eig, elbow = elb.b, M = mat): could not find function "pca.scree.plot"
# random (largely independent) data, usually higher elbow #
mat2 <- generate.test.matrix(5, 3)
#> Error in generate.test.matrix(5, 3): could not find function "generate.test.matrix"
result2 <- princomp(mat2)
#> Error: object 'mat2' not found
eig2 <- result2$sdev^2
#> Error: object 'result2' not found
elb2 <- quick_elbow(result2$sdev^2)
#> Error: object 'result2' not found
pca.scree.plot(eig2, elbow = elb2, M = mat2)
#> Error in pca.scree.plot(eig2, elbow = elb2, M = mat2): could not find function "pca.scree.plot"
```
