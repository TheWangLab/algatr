# Imputation of missing values using population structure inferred with `LEA::snmf`

Imputation of missing values using population structure inferred with
[`LEA::snmf`](https://rdrr.io/pkg/LEA/man/main_sNMF.html)

## Usage

``` r
str_impute(
  gen,
  K,
  entropy = TRUE,
  repetitions = 10,
  project = "new",
  quiet = TRUE,
  save_output = FALSE,
  output_filename = NULL
)
```

## Arguments

- gen:

  a dosage matrix, an object of class 'vcfR', or an object of type
  snmfProject

- K:

  An integer vector corresponding to the number of ancestral populations
  for which the snmf algorithm estimates have to be calculated.

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

- quiet:

  whether to operate quietly and suppress the results of cross-entropy
  scores (defaults to TRUE; only does so if more than one K-value); only
  displays run with minimum cross-entropy

- save_output:

  if TRUE, saves SNP GDS and ped (plink) files with retained SNPs in new
  directory; if FALSE returns object (defaults to FALSE)

- output_filename:

  if `save_output = TRUE`, name prefix for saved .geno file, sNMF
  project file, and sNMF output file results (defaults to FALSE, in
  which no files are saved)

## Value

dosage matrix with imputed missing values

## See also

Other Imputation functions:
[`gen_to_geno()`](https://thewanglab.github.io/algatr/reference/gen_to_geno.md),
[`geno_to_dosage()`](https://thewanglab.github.io/algatr/reference/geno_to_dosage.md),
[`simple_impute()`](https://thewanglab.github.io/algatr/reference/simple_impute.md),
[`snmf_bestK()`](https://thewanglab.github.io/algatr/reference/snmf_bestK.md),
[`snmf_crossent_helper()`](https://thewanglab.github.io/algatr/reference/snmf_crossent_helper.md)
