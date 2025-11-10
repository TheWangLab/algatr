# Install environmental and geographic data processing packages

Checks for the presence of packages required for environmental and
geographic data processing. If any of these packages are not already
installed, it will automatically install them.

## Usage

``` r
envirodata_packages()
```

## Value

None

## Details

The following packages will be installed if not already present:

- "RStoolbox" (from GitHub repository bleutner/RStoolbox)

- "geodata"

- "corrplot"

- "vegan"

- "gdistance"

- "topoDistance"

- "rmapshaper"

- "wingen"

## Examples

``` r
if (FALSE) envirodata_packages() # \dontrun{}
```
