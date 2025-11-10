# Convert dosage matrix or vcf to geno type object (N.B.: this only works for diploids!)

Convert dosage matrix or vcf to geno type object (N.B.: this only works
for diploids!)

## Usage

``` r
gen_to_geno(gen)
```

## Arguments

- gen:

  a dosage matrix, an object of class 'vcfR', or an object of type
  snmfProject

## Value

matrix encoded as geno type object

## See also

Other Imputation functions:
[`geno_to_dosage()`](https://thewanglab.github.io/algatr/reference/geno_to_dosage.md),
[`simple_impute()`](https://thewanglab.github.io/algatr/reference/simple_impute.md),
[`snmf_bestK()`](https://thewanglab.github.io/algatr/reference/snmf_bestK.md),
[`snmf_crossent_helper()`](https://thewanglab.github.io/algatr/reference/snmf_crossent_helper.md),
[`str_impute()`](https://thewanglab.github.io/algatr/reference/str_impute.md)
