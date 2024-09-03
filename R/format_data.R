#' Check if an object is a vcf or a path to a vcf
#'
#' @param x vcfR object or path to vcf
#'
#' @return vcf object
#' @export
#'
#' @noRd
vcf_check <- function(x) {
  if (class(x)[1] == "vcfR") {
    vcf <- x
  } else if (is.character(x)) {
    if (file.exists(x)) {
      vcf <- vcfR::read.vcfR(x)
    } else {
      stop("Cannot open file: No such file or directory")
    }
  } else {
    stop("Input is expected to be an object of class 'vcfR' or a path to a .vcf file")
  }

  return(vcf)
}

#' Remove islands from mapping
#'
#' @param input raster (SpatRaster or Raster*) or coordinates (sf or two column dataframe) to mask or remove island points from
#' @param shape SpatialPolygons, sf, or sf object to create island mask
#' @param min_vertices minimum number of vertices in polygons to retain (defaults to 10000)
#'
#' @return object (of input type) with islands removed
#' @export

rm_islands <- function(input, shape, min_vertices = 10000) {
  # Convert if SpatVector provided ------------------------------------------
  if (inherits(shape, "SpatVector")) shape <- sf::st_as_sf(shape)

  no_island <- rmapshaper::ms_filter_islands(shape, min_vertices = min_vertices)

  # convert SpatRaster to Raster
  if (inherits(input, "Raster")) input <- terra::rast(input)

  if (inherits(input, "SpatRaster")) {
    raster_noIsland <- terra::mask(input, no_island)
    return(raster_noIsland)
  }

  if (inherits(input, "data.frame") & !inherits(input, "sf")) {
    if (ncol(input) != 2) stop("Expected two columns (x and y)")
    colnames(input) <- c("x", "y")
    input <- sf::st_as_sf(input, coords = c("x", "y"), crs = sf::st_crs(shape))
  }

  if (inherits(input, "sf")) {
    input <- sf::st_transform(input, sf::st_crs(shape))
    coords_noIsland <- sf::st_intersection(input, no_island)
    return(coords_noIsland)
  }

}

#' Scale three layers of environmental data to R, G, and B for mapping
#'
#' @param env SpatRaster or Raster* with three layers
#'
#' @return RGB-scaled values
#' @export
scaleRGB <- function(env) {
  # Convert to SpatRaster if RasterStack provided ---------------------------
  if (!inherits(env, "SpatRaster")) env <- terra::rast(env)

  # Assign RGB values to each layer -----------------------------------------
  for (layer in 1:3) {
    minval <- terra::minmax(env[[layer]])[1, ]
    maxval <- terra::minmax(env[[layer]])[2, ]
    env[[layer]] <- ((env[[layer]] - minval) / (maxval - minval)) * 255
  }

  return(env)
}
