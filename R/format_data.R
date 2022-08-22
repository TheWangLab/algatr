

#' Convert a vcf to a dosage matrix
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#'
#' @return returns dosage matrix
#' @export
#'
#'
vcf_to_dosage <- function(x) {
  # check vcf
  vcf <- vcf_check(x)

  # convert to genlight
  genlight <- vcfR::vcfR2genlight(vcf)

  # convert to dosage matrix
  gen <- as.matrix(genlight)

  return(gen)
}


#' Check if an object is a vcf or a path to a vcf
#'
#' @param x vcfR object or path to vcf
#'
#' @return vcf object
#' @export
#'
#' @noRd
vcf_check <- function(x) {
  if (class(x)[1] == "vcfR")

    vcf <- x

  else if (is.character(x)) {

    if (file.exists(x))
      vcf <- vcfR::read.vcfR(x)

    else
      stop("Cannot open file: No such file or directory")

  } else stop("Input is expected to be an object of class 'vcfR' or a path to a .vcf file")

  return(vcf)
}


#' Remove islands from mapping
#' TODO: fill in param details
#' @param input
#' @param shape
#' @param min_vertices
#'
#' @return
#' @export
#'
#' @examples
rm_islands <- function(input, shape, min_vertices = 10000){
  no_island <- rmapshaper::ms_filter_islands(shape, min_vertices = min_vertices)

  if(class(input)[1] == "RasterLayer" | class(input)[1] == "RasterStack" ){
    raster_noIsland <- raster::mask(input, no_island)
    return(raster_noIsland)
  }

  if(class(input)[1] == "data.frame" & all(colnames(input) %in% c("ID", "x", "y"))){
    sp <- coords
    coordinates(sp) <- ~x+y
    crs(sp) <- raster::crs("+proj=longlat +datum=WGS84 +no_defs")

    sp_sub <- over(no_island, sp, returnList = TRUE)
    IDs <- sp_sub[[1]]$ID
    coords_noIsland <- coords[coords$ID %in% IDs,]
    return(coords_noIsland)
  }

}
