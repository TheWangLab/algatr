# Helper function to retrieve cross entropy scores from SNMF project

Helper function to retrieve cross entropy scores from SNMF project

## Usage

``` r
snmf_crossent_helper(snmf_proj, K, select_min = TRUE)
```

## Arguments

- snmf_proj:

  object of type snmfProject

- K:

  K-value(s)

- select_min:

  whether to return minimum

## Value

cross entropy scores for given K

## See also

Other Imputation functions:
[`gen_to_geno()`](https://thewanglab.github.io/algatr/reference/gen_to_geno.md),
[`geno_to_dosage()`](https://thewanglab.github.io/algatr/reference/geno_to_dosage.md),
[`simple_impute()`](https://thewanglab.github.io/algatr/reference/simple_impute.md),
[`snmf_bestK()`](https://thewanglab.github.io/algatr/reference/snmf_bestK.md),
[`str_impute()`](https://thewanglab.github.io/algatr/reference/str_impute.md)
