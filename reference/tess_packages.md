# Install TESS packages

Checks for the presence of packages required for TESS. If any of these
packages are not already installed, it will automatically install them.

## Usage

``` r
tess_packages()
```

## Value

None

## Details

The following packages will be installed if not already present:

- "automap"

- "graphics"

- "LEA" (from Bioconductor repository)

- "tess3r" (from GitHub repository bcm-uga/TESS3_encho_sen)

- "fields"

- "rworldmap"

- "cowplot"

## Examples

``` r
if (FALSE) tess_packages() # \dontrun{}
```
