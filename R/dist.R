#' Calculate distance between environmental vars
#'
#' @param env dataframe or vector of environmental variables for locations
#' @param stdz if TRUE then environmental values will be standardized (default = TRUE)
#'
#' @return list of environmental distances between samples (for each environmental variable)
#' @export
env_dist <- function(env, stdz = TRUE) {
  if (!is.null(dim)) distmat <- dplyr::as_tibble(env) %>% purrr::map(env_dist_helper, stdz = stdz)
  if (is.null(dim)) distmat <- env_dist_helper(env, stdz)
  return(distmat)
}

#' Helper function to convert an environmental vector to a distance matrix
#'
#' @inheritParams env_dist
#'
#' @export
#' @noRd
env_dist_helper <- function(env, stdz = TRUE) {
  # Standardize environmental variables
  if (stdz) env <- scale(env, center = TRUE, scale = TRUE)

  distmat <- as.matrix(dist(env, diag = TRUE, upper = TRUE))

  return(distmat)
}

#' Calculate geographic distance between coordinates
#'
#' @param coords dataframe with x and y coordinates
#' @param type the type of geographic distance to be calculated; options are "Euclidean" for direct distance, "topographic" for topographic distances, and "resistance" for resistance distances.
#' @param lyr SpatRaster or Raster* DEM for calculating topographic distances or resistance raster for calculating resistance distances (RasterLayer or SpatRaster object)
#' @details
#' Euclidean, or linear, distances are calculated using the geodist package: Padgham M, Sumner M (2021). geodist: Fast, Dependency-Free Geodesic Distance Calculations. R package version 0.0.7, Available: https://CRAN.R-project.org/package=geodist.
#' Topographic distances are calculated using the topoDistance package: Wang I.J. (2020) Topographic path analysis for modeling dispersal and functional connectivity: calculating topographic distances using the TOPODISTANCE R package. Methods in Ecology and Evolution, 11: 265-272.
#' Resistance distances are calculated using the gdistance package: van Etten, J. (2017). R package gdistance: Distances and routes on geographical grids. Journal of Statistical Software, 76(1), 1–21.
#'
#' @return geographic distance matrix
#' @export
geo_dist <- function(coords, type = "Euclidean", lyr = NULL) {
  if (type == "Euclidean" | type == "euclidean" | type == "linear") {
    # Format coordinates
    coords <- coords_to_sf(coords)
    # Calculate geodesic distance between points
    distmat <- sf::st_distance(coords)
  } else if (type == "topo" | type == "topographic") {
    # Format coordinates
    coords <- coords_to_df(coords)

    if (is.null(lyr)) stop("Calculating topographic distances requires a DEM layer for argument lyr.")
    message("Calculating topo distances... This can be time consuming with many points and large rasters.")

    # Convert to RasterLayer if SpatRaster object
    if (inherits(lyr, "SpatRaster")) lyr <- raster::raster(lyr)

    distmat <- topoDistance::topoDist(lyr, coords, paths = FALSE)
  } else if (type == "resistance" | type == "cost" | type == "res") {
    if (is.null(lyr)) stop("Calculating resistance distances requires a resistance surface for argument lyr.")
    message("Calculating resistance distances... This can be time consuming with many points and large rasters.")

    # Format coordinates
    coords <- coords_to_df(coords)

    # Convert to RasterLayer if SpatRaster object
    if (inherits(lyr, "SpatRaster")) lyr <- raster::raster(lyr)

    # Convert resistance surface to conductance surface
    cond.r <- 1 / lyr
    trSurface <- gdistance::transition(cond.r, transitionFunction = mean, directions = 8) # Create transition surface
    trSurface <- gdistance::geoCorrection(trSurface, type = "c", scl = FALSE)
    sp <- sp::SpatialPoints(coords = coords)
    distmat <- as.matrix(gdistance::commuteDistance(trSurface, sp)) # Calculate circuit distances
  }

  return(distmat)
}

#' convert coordinates to sf
#' @noRd
coords_to_sf <- function(coords) {
  if (inherits(coords, "sf")) {
    return(coords)
  }
  if (inherits(coords, "SpatVector")) {
    return(sf::st_as_sf(coords))
  }
  if (is.matrix(coords)) coords <- data.frame(coords)
  if (is.data.frame(coords)) colnames(coords) <- c("x", "y")
  return(sf::st_as_sf(coords, coords = c("x", "y")))
}
