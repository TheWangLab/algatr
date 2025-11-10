# K selection

K selection

## Usage

``` r
select_K(
  gen,
  K_selection = "tracy_widom",
  coords = NULL,
  Kvals = 1:10,
  criticalpoint = 2.023,
  low = 0.08,
  max.pc = 0.9,
  perc.pca = 90,
  max.n.clust = 10
)

select_K_tw(gen, criticalpoint = 2.0234)

select_K_elbow(gen, low = 0.08, max.pc = 0.9)

select_K_tess(
  gen,
  coords,
  Kvals = 1:10,
  tess_method = "projected.ls",
  ploidy = 2
)

select_K_fc(gen, perc.pca, max.n.clust)
```

## Arguments

- gen:

  a genotype matrix

- K_selection:

  method for performing K selection (options: "tracy_widom" (default),
  "quick_elbow", or "tess")

- coords:

  coordinates for "tess"

- Kvals:

  values of K to test for "tess"

- criticalpoint:

  if `K_selection = "tracy_widom"`, a numeric value corresponding to the
  significance level. If the significance level is 0.05, 0.01, 0.005, or
  0.001, the criticalpoint should be set to be 0.9793, 2.0234, 2.4224,
  or 3.2724, respectively (defaults to 2.0234)

- low:

  if `K_selection = "quick_elbow"`, numeric, between zero and one, the
  threshold that defines whether a principal component explains 'much'
  of the variance (defaults to 0.08).

- max.pc:

  if `K_selection = "quick_elbow"`, maximum percentage of the variance
  to capture before the elbow (cumulative sum to PC 'n'; defaults to
  0.90).

- perc.pca:

  if `K_selection = "find_clusters"`, a numeric value between 0 and 100
  indicating the minimal percentage of the total variance of the data to
  be expressed by the retained axes of PCA (defaults to 90).

- max.n.clust:

  if `K_selection = "find_clusters"`, an integer indicating the maximum
  number of clusters to try. Values of 'k' will be picked up between 1
  and max.n.clust (defaults to 10)

- tess_method:

  method to use for "tess"

- ploidy:

  ploidy for "tess"

## Value

prints the best K value given the specified K selection procedure

## Functions

- `select_K_tw()`: select K using Tracy-Widom Test

- `select_K_elbow()`: select K using PCA and `quick_elbow` method

- `select_K_tess()`: select K using TESS and `bestK` method

- `select_K_fc()`: select K using find.clusters method

## Note

uses the \`tw\` function (originally from the archived CRAN package
AssocTests)

uses the
[find.clusters](https://rdrr.io/pkg/adegenet/man/find.clusters.html)
function

## See also

Other LFMM functions:
[`lfmm_df()`](https://thewanglab.github.io/algatr/reference/lfmm_df.md),
[`lfmm_do_everything()`](https://thewanglab.github.io/algatr/reference/lfmm_do_everything.md),
[`lfmm_manhattanplot()`](https://thewanglab.github.io/algatr/reference/lfmm_manhattanplot.md),
[`lfmm_qqplot()`](https://thewanglab.github.io/algatr/reference/lfmm_qqplot.md),
[`lfmm_run()`](https://thewanglab.github.io/algatr/reference/lfmm_run.md),
[`lfmm_table()`](https://thewanglab.github.io/algatr/reference/lfmm_table.md),
[`lfmm_test_tidy()`](https://thewanglab.github.io/algatr/reference/lfmm_test_tidy.md),
[`quick_elbow()`](https://thewanglab.github.io/algatr/reference/quick_elbow.md)

Other LFMM functions:
[`lfmm_df()`](https://thewanglab.github.io/algatr/reference/lfmm_df.md),
[`lfmm_do_everything()`](https://thewanglab.github.io/algatr/reference/lfmm_do_everything.md),
[`lfmm_manhattanplot()`](https://thewanglab.github.io/algatr/reference/lfmm_manhattanplot.md),
[`lfmm_qqplot()`](https://thewanglab.github.io/algatr/reference/lfmm_qqplot.md),
[`lfmm_run()`](https://thewanglab.github.io/algatr/reference/lfmm_run.md),
[`lfmm_table()`](https://thewanglab.github.io/algatr/reference/lfmm_table.md),
[`lfmm_test_tidy()`](https://thewanglab.github.io/algatr/reference/lfmm_test_tidy.md),
[`quick_elbow()`](https://thewanglab.github.io/algatr/reference/quick_elbow.md)

Other LFMM functions:
[`lfmm_df()`](https://thewanglab.github.io/algatr/reference/lfmm_df.md),
[`lfmm_do_everything()`](https://thewanglab.github.io/algatr/reference/lfmm_do_everything.md),
[`lfmm_manhattanplot()`](https://thewanglab.github.io/algatr/reference/lfmm_manhattanplot.md),
[`lfmm_qqplot()`](https://thewanglab.github.io/algatr/reference/lfmm_qqplot.md),
[`lfmm_run()`](https://thewanglab.github.io/algatr/reference/lfmm_run.md),
[`lfmm_table()`](https://thewanglab.github.io/algatr/reference/lfmm_table.md),
[`lfmm_test_tidy()`](https://thewanglab.github.io/algatr/reference/lfmm_test_tidy.md),
[`quick_elbow()`](https://thewanglab.github.io/algatr/reference/quick_elbow.md)

Other LFMM functions:
[`lfmm_df()`](https://thewanglab.github.io/algatr/reference/lfmm_df.md),
[`lfmm_do_everything()`](https://thewanglab.github.io/algatr/reference/lfmm_do_everything.md),
[`lfmm_manhattanplot()`](https://thewanglab.github.io/algatr/reference/lfmm_manhattanplot.md),
[`lfmm_qqplot()`](https://thewanglab.github.io/algatr/reference/lfmm_qqplot.md),
[`lfmm_run()`](https://thewanglab.github.io/algatr/reference/lfmm_run.md),
[`lfmm_table()`](https://thewanglab.github.io/algatr/reference/lfmm_table.md),
[`lfmm_test_tidy()`](https://thewanglab.github.io/algatr/reference/lfmm_test_tidy.md),
[`quick_elbow()`](https://thewanglab.github.io/algatr/reference/quick_elbow.md)

Other LFMM functions:
[`lfmm_df()`](https://thewanglab.github.io/algatr/reference/lfmm_df.md),
[`lfmm_do_everything()`](https://thewanglab.github.io/algatr/reference/lfmm_do_everything.md),
[`lfmm_manhattanplot()`](https://thewanglab.github.io/algatr/reference/lfmm_manhattanplot.md),
[`lfmm_qqplot()`](https://thewanglab.github.io/algatr/reference/lfmm_qqplot.md),
[`lfmm_run()`](https://thewanglab.github.io/algatr/reference/lfmm_run.md),
[`lfmm_table()`](https://thewanglab.github.io/algatr/reference/lfmm_table.md),
[`lfmm_test_tidy()`](https://thewanglab.github.io/algatr/reference/lfmm_test_tidy.md),
[`quick_elbow()`](https://thewanglab.github.io/algatr/reference/quick_elbow.md)
