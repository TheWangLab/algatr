# Convert lfmm/geno matrix to dosage matrix (N.B.: this only works for diploids!)

Convert lfmm/geno matrix to dosage matrix (N.B.: this only works for
diploids!)

## Usage

``` r
geno_to_dosage(geno)
```

## Arguments

- geno:

  matrix of LEA geno or lfmm format (i.e., 0 corresponds to zero
  reference alleles)

## Value

matrix encoded as dosage type object (0 corresponds to two reference
alleles)

## See also

Other Imputation functions:
[`gen_to_geno()`](https://thewanglab.github.io/algatr/reference/gen_to_geno.md),
[`simple_impute()`](https://thewanglab.github.io/algatr/reference/simple_impute.md),
[`snmf_bestK()`](https://thewanglab.github.io/algatr/reference/snmf_bestK.md),
[`snmf_crossent_helper()`](https://thewanglab.github.io/algatr/reference/snmf_crossent_helper.md),
[`str_impute()`](https://thewanglab.github.io/algatr/reference/str_impute.md)
