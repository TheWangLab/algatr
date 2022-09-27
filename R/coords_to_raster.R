
#' Create a raster from coordinates
#'
#' @param coords x and y coordinates (two columns, the first should be x and the second should be y)
#' @param buff buffer to add to edge of raster (default = 0)
#' @param res desired resolution of raster (default = NULL). Can be a single value for square cells or a vector with two values representing x and y resolutions.
#' @param agg aggregation factor to apply to raster (default = NULL)
#' @param disagg disaggregation factor to apply to raster (default = NULL)
#' @param plot whether to plot resulting raster with coords (default = FALSE)
#'
#' @return RasterLayer
#' @export
#'
#' @examples
#' load_mini_ex()
#' coords_to_raster(mini_coords, buffer = 1, plot = TRUE)
coords_to_raster <- function(coords, buff = 0, res = NULL, agg = NULL, disagg = NULL, plot = FALSE) {

  # Make a matrix
  r <- make_raster(coords, buff = buff, res = res)

  # Aggregate or disaggregate
  if (!is.null(agg) & !is.null(disagg)) {
    warning("both agg and disagg were provided. Did you mean to do this? (if so, note that aggregation will occur first and then disaggregation second")
  }
  if (!is.null(agg)) r <- raster::aggregate(r, agg)
  if (!is.null(disagg)) r <- raster::disaggregate(r, disagg)

  # Assign values to make it easier to visualize the resolution
  r <- raster::init(r)
  r[] <- 1:raster::ncell(r)

  # Plot raster
  if (plot) {
    raster::plot(r, legend = FALSE, col = viridis::mako(raster::ncell(r)))
    graphics::points(coords, col = viridis::magma(1, begin = 0.7), pch = 3, lwd = 2)
  }

  return(r)
}

#' coords to raster converter
#'
#' @inheritParams coords_to_raster
#'
#' @export
#' @noRd
make_raster <- function(coords, buff = 0, res = NULL){
  # Format coords
  coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")

  # Get x and y min and max and round up to nearest integer
  # (Note: must be an integer for assigning nrow and ncol of a matrix)
  xmin <- ceiling(min(coords$x, na.rm = TRUE) - buff)
  xmax <- ceiling(max(coords$x, na.rm = TRUE) + buff)
  ymin <- ceiling(min(coords$y, na.rm = TRUE) - buff)
  ymax <- ceiling(max(coords$y, na.rm = TRUE) + buff)

  # Make matrix
  m <- matrix(nrow = (ymax - ymin), ncol = (xmax - xmin))

  # Turn into raster
  r <- raster::raster(m)

  # Set extent
  raster::extent(r) <- c(xmin, xmax, ymin, ymax)

  # Set resolution
  if(length(res) > 2) stop("invalid res provided")
  if(!is.null(res)) raster::res(r) <- res

  return(r)
}
