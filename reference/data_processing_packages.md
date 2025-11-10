# Install data processing packages

Checks for the presence of packages required for genetic data
processing. If any of these packages are not already installed, it will
automatically install them.

## Usage

``` r
data_processing_packages()
```

## Value

None

## Details

The following packages will be installed if not already present:

- "gdsfmt" (from Bioconductor repository)

- "SeqArray" (from Bioconductor repository)

- "SNPRelate" (from Bioconductor repository)

## Examples

``` r
if (FALSE) data_processing_packages() # \dontrun{}
```
