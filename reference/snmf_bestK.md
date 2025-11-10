# Helper function to select "best" K based on minimizing cross-entropy criteria from sNMF results

Helper function to select "best" K based on minimizing cross-entropy
criteria from sNMF results

## Usage

``` r
snmf_bestK(snmf_proj, K, quiet)
```

## Arguments

- snmf_proj:

  object of type snmfProject

- K:

  integer corresponding to K-value

- quiet:

  whether to operate quietly and suppress the results of cross-entropy
  scores (defaults to TRUE; only does so if more than one K-value); only
  displays run with minimum cross-entropy

## Value

list with best K-value and run number and all cross-entropy scores

## See also

Other Imputation functions:
[`gen_to_geno()`](https://thewanglab.github.io/algatr/reference/gen_to_geno.md),
[`geno_to_dosage()`](https://thewanglab.github.io/algatr/reference/geno_to_dosage.md),
[`simple_impute()`](https://thewanglab.github.io/algatr/reference/simple_impute.md),
[`snmf_crossent_helper()`](https://thewanglab.github.io/algatr/reference/snmf_crossent_helper.md),
[`str_impute()`](https://thewanglab.github.io/algatr/reference/str_impute.md)
