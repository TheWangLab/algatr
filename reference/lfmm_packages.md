# Install LFMM packages

Checks for the presence of packages required for LFMM. If any of these
packages are not already installed, it will automatically install them.

## Usage

``` r
lfmm_packages()
```

## Value

None

## Details

The following packages will be installed if not already present:

- "adegenet"

- "lfmm"

- "tess3r" (from GitHub repository bcm-uga/TESS3_encho_sen)

- "LEA" (from Bioconductor repository)

## Examples

``` r
if (FALSE) lfmm_packages() # \dontrun{}
```
