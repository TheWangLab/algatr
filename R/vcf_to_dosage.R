
#' Convert a vcf to a dosage matrix
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#'
#' @return dosage matrix
#' @export
vcf_to_dosage <- function(x) {
  genlight <- vcfR::vcfR2genlight(x)
  gen <- as.matrix(genlight)
  return(gen)
}
