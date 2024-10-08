
#' Lazy run of all landscape genomic analyses contained within `algatr`
#'
#' Disclaimer: this is probably a bad idea...
#'
#' @param gen path to vcf file, a `vcfR` type object, or a dosage matrix
#' @param coords dataframe with x (i.e., longitude) and y (i.e., latitude) coordinates; must be in this order
#' @param envlayers SpatRaster or Raster* for mapping (if env is provided, the dataframe column names and envlayers layer names should be the same)
#' @param quiet whether to operate quietly and suppress the output of tables and figures (defaults to FALSE)
#' @param gators set to TRUE to see some gators...
#'
#' @return results from all six analyses contained within algatr
#' @export
do_everything_for_me <- function(gen, coords, envlayers, quiet = FALSE, gators = FALSE) {

  message("Please be aware: the do_everything functions are meant to be exploratory. We do not recommend their use for final analyses unless certain they are properly parameterized.")

  # Data processing ---------------------------------------------------------

  gendist <- gen_dist(gen, dist_type = "euclidean")

  lyr <- wingen::coords_to_raster(coords, res = 0.25, buffer = 5)

  # Run algatr methods ------------------------------------------------------

  wingen <- wingen_do_everything(gen, lyr, coords, kriged = TRUE, grd = lyr, masked = TRUE, mask = envlayers, quiet = quiet)

  if (gators & !quiet) ascii_alligator("wingen")

  tess <- tess_do_everything(gen, coords, raster::aggregate(envlayers[[1]], 10), Kvals = 1:10, K_selection = "auto", quiet = quiet)

  if (gators & !quiet) ascii_alligator("TESS")

  mmrr <- mmrr_do_everything(gendist, coords, envlayers, model = "best", nperm = 100, quiet = quiet)

  if (is.null(mmrr)) {
    warning("MMRR model = \"best\" did not find a significant model, running a full model instead")
    mmrr <- quiet_mmrr(gendist, coords, envlayers, model = "full", nperm = 100, quiet = quiet)
  }

  if (gators & !quiet) ascii_alligator("MMRR")

  gdm <- gdm_do_everything(gendist, coords, envlayers, model = "best", scale_gendist = TRUE, quiet = quiet)

  if (is.null(gdm)) {
    warning("GDM model = \"best\" did not find a significant model, running a full model instead")
    gdm <- gdm_do_everything(gendist, coords, envlayers, model = "full", scale_gendist = TRUE, quiet = quiet)
  }

  if (gators & !quiet) ascii_alligator("GDM")

  rda <- rda_do_everything(gen, envlayers, coords, quiet = quiet)

  if (gators & !quiet) ascii_alligator("RDA")

  lfmm <- lfmm_do_everything(gen, envlayers, coords, quiet = quiet)

  if (gators & !quiet) ascii_alligator("LFMM")

  return(list(wingen = wingen, tess = tess, gdm = gdm, mmrr = mmrr, rda = rda, lfmm = lfmm))
}
