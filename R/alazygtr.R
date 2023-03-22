
#' Lazy run of all landscape genomic analyses contained within `algatr`
#'
#' Disclaimer: this is probably a bad idea...
#'
#' @param vcf path to vcf file or a `vcfR` type object
#' @param coords dataframe with x (i.e., longitude) and y (i.e., latitude) coordinates; must be in this order
#' @param envlayers envlayers for mapping (if env is provided, the dataframe column names and envlayers layer names should be the same)
#'
#' @return
#' @export
#'
#' @examples
do_everything_for_me <- function(vcf, coords, envlayers){

  gendist <- gen_dist(vcf = vcf, dist_type = "dps")

  tess <- tess_do_everything(vcf, coords, raster::aggregate(envlayers[[1]], 10), Kvals = 1:10, K_selection = "auto")

  ascii_alligator("TESS")

  mmrr <- mmrr_do_everything(gendist, coords, envlayers, model = "best", nperm = 100)

  ascii_alligator("MMRR")

  if(is.null(mmrr)) {
    warning("MMRR model = \"best\" did not find a significant model, running a full model instead")
    mmrr <- gdm_do_everything(gendist, liz_coords, CA_env, model = "full", scale_gendist = TRUE, nperm = 100)
    }

  gdm <- gdm_do_everything(gendist, coords, envlayers, model = "best", scale_gendist = TRUE)

  if(is.null(gdm)) {
    warning("GDM model = \"best\" did not find a significant model, running a full model instead")
    gdm <- gdm_do_everything(gendist, coords, envlayers, model = "full", scale_gendist = TRUE)
    }

  ascii_alligator("GDM")

  rda <- rda_do_everything(vcf, envlayers, coords)

  ascii_alligator("RDA")

  lfmm <- lfmm_do_everything(vcf, envlayers, coords)

  ascii_alligator("LFMM")

  return(list(gdm = gdm, mmrr = mmrr, rda = rda, lfmm = lfmm))
}
