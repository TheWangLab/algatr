# Impute NA values NOTE: use extreme caution when using this form of simplistic imputation. We mainly provide this code for creating test datasets and highly discourage its use in analyses.

Impute NA values NOTE: use extreme caution when using this form of
simplistic imputation. We mainly provide this code for creating test
datasets and highly discourage its use in analyses.

## Usage

``` r
simple_impute(x, FUN = median)
```

## Arguments

- x:

  matrix

- f:

  function to use for imputation (defaults to median)

## Value

matrix of values with missing values imputed

## See also

Other Imputation functions:
[`gen_to_geno()`](https://thewanglab.github.io/algatr/reference/gen_to_geno.md),
[`geno_to_dosage()`](https://thewanglab.github.io/algatr/reference/geno_to_dosage.md),
[`snmf_bestK()`](https://thewanglab.github.io/algatr/reference/snmf_bestK.md),
[`snmf_crossent_helper()`](https://thewanglab.github.io/algatr/reference/snmf_crossent_helper.md),
[`str_impute()`](https://thewanglab.github.io/algatr/reference/str_impute.md)
