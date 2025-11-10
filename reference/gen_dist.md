# Calculate genetic distances

Calculate genetic distances

## Usage

``` r
gen_dist(
  gen = NULL,
  dist_type = "euclidean",
  plink_file = NULL,
  plink_id_file = NULL,
  npc_selection = "auto",
  criticalpoint = 2.0234
)
```

## Arguments

- gen:

  path to vcf file, a `vcfR` type object, or a dosage matrix

- dist_type:

  the type of genetic distance to calculate (options: `"euclidean"`
  (default), `"bray_curtis"`, `"dps"` for proportion of shared alleles
  (requires vcf), `"plink"`, or `"pc"` for PC-based)

- plink_file:

  if `"plink"` dist_type is used, path to plink distance file (typically
  ".dist"; required only for calculating plink distance). File must be a
  **square** distance matrix.

- plink_id_file:

  if `"plink"` dist_type is used, path to plink id file (typically
  ".dist.id"; required only for calculating plink distance)

- npc_selection:

  if `dist_type = "pc"`, how to perform K selection (options: `"auto"`
  for automatic selection based on significant eigenvalues from
  Tracy-Widom test (default), or `"manual"` to examine PC screeplot and
  enter no. PCs into console)

- criticalpoint:

  if `dist_type = "pc"` used with `npc_selection = "auto"`, the critical
  point for the significance threshold for the Tracy-Widom test within
  the PCA (defaults to 2.0234 which corresponds to an alpha of 0.01)

## Value

pairwise distance matrix for given distance metric

## Details

Euclidean and Bray-Curtis distances calculated using the ecodist
package: Goslee, S.C. and Urban, D.L. 2007. The ecodist package for
dissimilarity-based analysis of ecological data. Journal of Statistical
Software 22(7):1-19. DOI:10.18637/jss.v022.i07. Proportions of shared
alleles calculated using the adegenet package: Jombart T. and Ahmed I.
(2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP
data. Bioinformatics. doi:10.1093/bioinformatics/btr521. For calculating
proportions of shared alleles, missing values are ignored (i.e., prop
shared alleles calculated from present values; no scaling performed)
