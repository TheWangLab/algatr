
#' Lazy run of all landscape genomic analyses contained within `algatr`
#'
#' Disclaimer: this is probably a bad idea...
#'
#' @param gen path to vcf file, a `vcfR` type object, or a dosage matrix
#' @param coords dataframe with x (i.e., longitude) and y (i.e., latitude) coordinates; must be in this order
#' @param envlayers envlayers for mapping (if env is provided, the dataframe column names and envlayers layer names should be the same)
#' @param quiet whether to print output tables and figures (defaults to FALSE)
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

  quiet_wingen <- purrr::quietly(wingen_do_everything)
  wingen <- quiet_wingen(gen, lyr, coords, kriged = TRUE, grd = lyr, masked = TRUE, mask = envlayers, quiet = quiet)$result

  if (gators & !quiet) ascii_alligator("wingen")

  quiet_tess <- purrr::quietly(tess_do_everything)
  tess <- quiet_tess(gen, coords, raster::aggregate(envlayers[[1]], 10), Kvals = 1:10, K_selection = "auto", quiet = quiet)$result

  if (gators & !quiet) ascii_alligator("TESS")

  quiet_mmrr <- purrr::quietly(mmrr_do_everything)
  mmrr <- quiet_mmrr(gendist, coords, envlayers, model = "best", nperm = 100, quiet = quiet)$result

  if (is.null(mmrr)) {
    warning("MMRR model = \"best\" did not find a significant model, running a full model instead")
    mmrr <- quiet_mmrr(gendist, coords, envlayers, model = "full", nperm = 100, quiet = quiet)$result
  }

  if (gators & !quiet) ascii_alligator("MMRR")

  quiet_gdm <- purrr::quietly(gdm_do_everything)
  gdm <- quiet_gdm(gendist, coords, envlayers, model = "best", scale_gendist = TRUE, quiet = quiet)$result

  if (is.null(gdm)) {
    warning("GDM model = \"best\" did not find a significant model, running a full model instead")
    gdm <- quiet_gdm(gendist, coords, envlayers, model = "full", scale_gendist = TRUE, quiet = quiet)$result
  }

  if (gators & !quiet) ascii_alligator("GDM")

  quiet_rda <- purrr::quietly(rda_do_everything)
  rda <- quiet_rda(gen, envlayers, coords, quiet = quiet)$result

  if (gators & !quiet) ascii_alligator("RDA")

  quiet_lfmm <- purrr::quietly(lfmm_do_everything)
  lfmm <- quiet_lfmm(gen, envlayers, coords, quiet = quiet)$result

  if (gators & !quiet) ascii_alligator("LFMM")

  return(list(wingen = wingen, tess = tess, gdm = gdm, mmrr = mmrr, rda = rda, lfmm = lfmm))
}
