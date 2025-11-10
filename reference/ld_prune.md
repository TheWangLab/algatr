# ld_prune prunes SNPs based on linkage disequilibrium using `SNPRelate` and `SeqArray` packages

ld_prune prunes SNPs based on linkage disequilibrium using `SNPRelate`
and `SeqArray` packages

## Usage

``` r
ld_prune(
  vcf,
  out_name,
  out_format,
  nodes = 1,
  ld.threshold = 0.6,
  slide.max.n = 100,
  maf = 0.05,
  seed = 1234,
  method = c("corr"),
  save_output = TRUE
)
```

## Arguments

- vcf:

  is the path to the vcf file containing all SNPs

- out_name:

  prefix name of output files (will append with param settings)

- out_format:

  output file format ("plink" will produce ped and map files while "vcf"
  will produce a vcf and a GDS)

- nodes:

  is the number of parallel processors (numeric)

- ld.threshold:

  is the threshold for LD pruning (numeric; 0 - 1; defaults to 0.6)

- slide.max.n:

  is the maximum number of SNPs in a sliding window (numeric; defaults
  to 100)

- maf:

  is the minor allele frequency cutoff (numeric; defaults to 0.05)

- seed:

  is the random starting seed (defaults to 1234)

- method:

  is the LD threshold method; default to corr which is r2 correlation
  coefficient

- save_output:

  if TRUE, saves SNP GDS and ped (plink) files with retained SNPs in new
  directory; if FALSE returns object (defaults to TRUE)

## Value

LD-pruned vcf-type object

## Details

`SNPRelate` package citation: Zheng et al. (2012):
https://doi.org/10.1093/bioinformatics/bts606
