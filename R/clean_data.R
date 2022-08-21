
#' Convert vcf to dosage matrix
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#'
#' @return returns dosage matrix
#' @export
#'
#' @examples
#' data("ex_vcf")
#' vcf_to_dosage(ex_vcf)
#'
vcf_to_dosage <- function(x){

  if(class(x)[1] == "vcfR"){
    vcf <- x
  } else if(file.exists(x)){
    vcf <- vcfR::read.vcfR(x)
  }

  # convert to genlight
  genlight <- vcfR::vcfR2genlight(vcf)

  # convert to dosage matrix
  gen <- as.matrix(genlight)

  return(gen)
}

#' Wrapper for vcfR2genind function that assigns pops
#'
#' @param x can either be an object of class 'vcfR' or a path to a .vcf file
#' @param pops if NULL (default), and there are no pops detected from the vcf, each individual is assigned its own pop. If FALSE then genind$pop is left NULL. Alternatively, a vector of population assignments for each individual can be provided
#'
#' @return returns genind object
#' @export
#'
#' @examples
#'
#' data("ex_vcf")
#' vcf_to_genind(ex_vcf)
#'
vcf_to_genind <- function(x, pops = NULL){

  if(class(x)[1] == "vcfR"){
    vcf <- x
  } else if(file.exists(x)){
    vcf <- vcfR::read.vcfR(x)
  }

  # Convert to genind object
  genind <- vcfR::vcfR2genind(vcf)

  # Assign pops if null or pop vector provided
  if(is.null(genind$pop) | is.vector(pops)){

    if(is.null(pops)){
      genind$pop <- as.factor(1:nrow(genind@tab))

      warning("no pops were provided, assigning a pop to each individual (to stop this, set pop = FALSE)")
    }

    if(is.vector(pops)){

      if(!is.null(genind$pop) & is.vector(pops)){

        warning("overwriting genind pops with vector of pops provided")

      }

      if(length(pops) != nrow(genind@tab)){

        stop("length of pops does not match number of individuals in genind")

      } else {

        genind$pop <- as.factor(pops)

      }
    }

  }

  return(genind)
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
