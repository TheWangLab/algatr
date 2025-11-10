# Lazy run of all landscape genomic analyses contained within `algatr`

Disclaimer: this is probably a bad idea...

## Usage

``` r
do_everything_for_me(gen, coords, envlayers, quiet = FALSE, gators = FALSE)
```

## Arguments

- gen:

  path to vcf file, a `vcfR` type object, or a dosage matrix

- coords:

  dataframe with x (i.e., longitude) and y (i.e., latitude) coordinates;
  must be in this order

- envlayers:

  SpatRaster or Raster\* for mapping (if env is provided, the dataframe
  column names and envlayers layer names should be the same)

- quiet:

  whether to operate quietly and suppress the output of tables and
  figures (defaults to FALSE)

- gators:

  set to TRUE to see some gators...

## Value

results from all six analyses contained within algatr
