
#' Perform Rastar PCA
#'
#' Runs Principal Component Analysis (PCA) on a Raster stack.
#'
#' @param x A Raster* or SpatRaster object
#' @param scale a logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place. The default is TRUE, unlike for \link[stats]{prcomp}, since in most cases scaling should be performed.
#' @inheritParams stats::prcomp
#'
#' @return A SpatRaster of principal component scores
#'
#' @examples
#' library(terra)
#'
#' # Create a Raster stack with random data
#' set.seed(123)
#' r <- c(rast(matrix(rnorm(100), nrow = 10)),  rast(matrix(rnorm(100), nrow = 10)))
#'
#' pca_result <- raster_pca(r)
#'
#' @export
raster_pca <- function(x, center = TRUE, scale = TRUE, tol = NULL, rank. = NULL, ...) {
  # Check and convert x to a SpatRaster
  if (!inherits(x, "SpatRaster")) x <- terra::rast(x)

  # Convert to a data frame, removing NA values if present
  df <- terra::as.data.frame(x, na.rm = TRUE)

  # Perform PCA
  pca <- stats::prcomp(df, center = center, scale = scale, tol = tol, rank. = rank., ...)

  # Apply the PCA model to the raster to produce raster PCs
  env_pcs <- terra::predict(x, pca)

  return(env_pcs)
}


