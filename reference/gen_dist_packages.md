# Install genetic distance packages

Checks for the presence of packages required for genetic distance
calculations. If any of these packages are not already installed, it
will automatically install them.

## Usage

``` r
gen_dist_packages()
```

## Value

None

## Details

The following packages will be installed if not already present:

- "adegenet"

- "readr"

- "tibble"

- "ecodist"

- "cowplot"

## Examples

``` r
if (FALSE) gen_dist_packages() # \dontrun{}
```
