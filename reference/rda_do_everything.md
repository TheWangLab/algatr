# RDA function to do everything

RDA function to do everything

## Usage

``` r
rda_do_everything(
  gen,
  env,
  coords = NULL,
  impute = "structure",
  K_impute = 3,
  entropy = TRUE,
  repetitions = 10,
  project = "new",
  quiet_impute = TRUE,
  save_output = FALSE,
  output_filename = NULL,
  model = "full",
  correctGEO = FALSE,
  correctPC = FALSE,
  outlier_method = "p",
  sig = 0.05,
  z = 3,
  p_adj = "fdr",
  cortest = TRUE,
  nPC = 3,
  varpart = FALSE,
  naxes = "all",
  Pin = 0.05,
  R2permutations = 1000,
  R2scope = T,
  stdz = TRUE,
  quiet = FALSE
)
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

- impute:

  if NAs in `gen`, imputation will be performed on missing values;
  options are "structure" which uses the
  [`str_impute()`](https://thewanglab.github.io/algatr/reference/str_impute.md)
  function to impute based on population structure inferred with
  [`LEA::snmf`](https://rdrr.io/pkg/LEA/man/main_sNMF.html) (default);
  other option is "simple" based on
  [`simple_impute()`](https://thewanglab.github.io/algatr/reference/simple_impute.md)
  which imputes to the median

- K_impute:

  if `impute = "structure"`, an integer vector (range or single value)
  corresponding to the number of ancestral populations for which the
  sNMF algorithm estimates have to be calculated (defaults to 3)

- entropy:

  A boolean value. If true, the cross-entropy criterion is calculated
  (see
  [`create.dataset`](https://rdrr.io/pkg/LEA/man/main_createDataSet.html)
  and
  [`cross.entropy.estimation`](https://rdrr.io/pkg/LEA/man/main_crossEntropyEstimation.html)).

- repetitions:

  An integer corresponding with the number of repetitions for each value
  of `K`.

- project:

  A character string among "continue", "new", and "force". If
  "continue", the results are stored in the current project. If "new",
  the current project is removed and a new one is created to store the
  result. If "force", the results are stored in the current project even
  if the input file has been modified since the creation of the project.

- quiet_impute:

  if `impute = "structure"`, whether to suppress results of
  cross-entropy scores (defaults to TRUE; only does so if K is range of
  values); only displays run with minimum cross-entropy

- save_output:

  if `impute = "structure"`, if TRUE, saves SNP GDS and ped (plink)
  files with retained SNPs in new directory; if FALSE returns object
  (defaults to FALSE)

- output_filename:

  if `impute = "structure"` and `save_output = TRUE`, name prefix for
  saved .geno file, SNMF project file, and SNMF output file results
  (defaults to FALSE, in which no files are saved)

- model:

  whether to fit the model with all variables ("full") or to perform
  variable selection to determine the best set of variables ("best");
  defaults to "full"

- correctGEO:

  whether to condition on geographic coordinates

- correctPC:

  whether to condition on PCs from PCA of genotypes

- outlier_method:

  method to determine outliers. Can either be "p" to use the p-value
  method from
  [here](https://github.com/Capblancq/RDA-landscape-genomics) or "z" to
  use the z-score based method from
  [here](https://popgen.nescent.org/2018-03-27_RDA_GEA.html)

- sig:

  if `outlier_method = "p"`, the significance level to use to identify
  SNPs (defaults to 0.05)

- z:

  if `outlier_method = "z"`, the number of standard deviations to use to
  identify SNPs (defaults to 3)

- p_adj:

  if `outlier_method = "p"`, method to use for p-value correction
  (defaults to "fdr"); other options can be found in
  [`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)

- cortest:

  whether to create table of correlations for SNPs and environmental
  variable (defaults to TRUE)

- nPC:

  number of PCs to use if correctPC = TRUE (defaults to 3); if set to
  "manual" a selection option with a terminal prompt will be provided

- varpart:

  whether to perform variance partitioning (defaults to FALSE)

- naxes:

  number of RDA axes to use (defaults to "all" to use all axes), if set
  to "manual" a selection option with a terminal prompt will be given,
  otherwise can be any integer that is less than or equal to the total
  number of axes

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

- stdz:

  whether to center and scale environmental data (defaults to TRUE)

- quiet:

  whether to operate quietly and suppress the output of tables and
  figures (defaults to FALSE)

## Value

list containing (1) outlier SNPs, (2) dataframe with correlation test
results, if `cortest = TRUE`, (3) the RDA model, (4) results from
outlier analysis (output from
[rda_getoutliers](https://thewanglab.github.io/algatr/reference/rda_getoutliers.md)),
(5) RDA R-Squared, (6) RDA ANOVA, (7) p-values if
`outlier_method = "p"`, and (8) results from variance partitioning
analysis, if `varpart = TRUE`

## Details

Much of algatr's code is adapted from Capblancq T., Forester B.R. 2021.
Redundancy analysis: A swiss army knife for landscape genomics. Methods
Ecol. Evol. 12:2298-2309. doi: https://doi.org/10.1111/2041-210X.13722.

## See also

Other RDA functions:
[`rda_cor()`](https://thewanglab.github.io/algatr/reference/rda_cor.md),
[`rda_getoutliers()`](https://thewanglab.github.io/algatr/reference/rda_getoutliers.md),
[`rda_plot()`](https://thewanglab.github.io/algatr/reference/rda_plot.md),
[`rda_run()`](https://thewanglab.github.io/algatr/reference/rda_run.md),
[`rda_table()`](https://thewanglab.github.io/algatr/reference/rda_table.md),
[`rda_varpart()`](https://thewanglab.github.io/algatr/reference/rda_varpart.md),
[`rda_varpart_table()`](https://thewanglab.github.io/algatr/reference/rda_varpart_table.md)
