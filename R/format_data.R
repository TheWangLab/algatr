
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
#'
#' @param input RasterLayer or RasterStack object with islands to be removed; also accepts coords
#' @param shape SpatialPolygons, sf, sfc, or other polygon object to filter; also accepts SpatVector object
#' @param min_vertices minimum number of vertices in polygons to retain (defaults to 10000)
#'
#' @return object (of input type) with islands removed
#' @export
#'
#' @examples
rm_islands <- function(input, shape, min_vertices = 10000){

  # Convert if SpatVector provided ------------------------------------------
  if(inherits(shape, "SpatVector")) shape <- sf::st_as_sf(shape)

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


#' Impute NA values
#' NOTE: use extreme caution when using this form of simplistic imputation. We mainly provide this code for creating test datasets and highly discourage its use in analyses.
#' @param x matrix
#' @param f function to use for imputation (defaults to median)
#'
#' @return
#' @export
#'
#' @examples
simple_impute <- function(x, FUN = median){
  x_noNA <- apply(x, 2, impute_helper, FUN)
  return(x_noNA)
}

#' Helper function for imputation
#' @export
#' @noRd
impute_helper <- function(i, FUN = median){
  i[which(is.na(i))] <- FUN(i[-which(is.na(i))], na.rm = TRUE)
  return(i)
}

#' Scale three layers of environmental data to R, G, and B for mapping
#'
#' @param env object of type SpatRaster or RasterStack
#'
#' @return
#' @export
#'
#' @examples
scaleRGB <- function(env){

  # Convert to SpatRaster if RasterStack provided ---------------------------
  if (!inherits(env, "SpatRaster")) env <- terra::rast(env)

  # Assign RGB values to each layer -----------------------------------------
  for(layer in 1:3){
    minval <- terra::minmax(env[[layer]])[1,]
    maxval <- terra::minmax(env[[layer]])[2,]
    env[[layer]] <- ((env[[layer]] - minval) / (maxval - minval))*255
  }

  return(env)
}
