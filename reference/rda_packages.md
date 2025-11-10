# Install RDA packages

Checks for the presence of packages required for RDA. If any of these
packages are not already installed, it will automatically install them.

## Usage

``` r
rda_packages()
```

## Value

None

## Details

The following packages will be installed if not already present:

- "ggrepel"

- "qvalue" (from Bioconductor repository)

- "robust"

- "tibble"

- "vegan"

## Examples

``` r
if (FALSE) rda_packages() # \dontrun{}
```
