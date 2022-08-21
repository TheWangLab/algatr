#' Create a raster from coordinates
#'
#' Generate a raster layer from coordinates which can be used in \link[wingen]{window_gd} as the RasterLayer to move the window across
#'
#' @param coords coordinates (two columns, the first should be x and the second should be y)
#' @param buffer buffer to add to edge of raster (defaults to 0)
#' @param res desired resolution of raster (defaults to NULL). Can be a single value for square cells or a vector with two values representing x and y resolutions.
#' @param agg aggregation factor to apply to raster (defaults to NULL)
#' @param disagg disaggregation factor to apply to raster (defaults to NULL)
#' @param plot whether to plot resulting raster with coords (defaults to FALSE)
#'
#' @return RasterLayer
#' @export
#'
#' @examples
#' load_mini_ex()
#' coords_to_raster(mini_coords, buffer = 1, plot = TRUE)
coords_to_raster <- function(coords, buffer = 0, res = NULL, agg = NULL, disagg = NULL, plot = FALSE) {

  # make a matrix
  r <- make_raster(coords, buffer = buffer, res = res)

  # aggregate or disaggregate
  if (!is.null(agg) & !is.null(disagg)) {
    warning("both agg and disagg were provided. Did you mean to do this? (if so, note that aggregation will occur first and then disaggregation second")
  }
  if (!is.null(agg)) r <- raster::aggregate(r, agg)
  if (!is.null(disagg)) r <- raster::disaggregate(r, disagg)

  # assign values to make it easier to visualize the resolution
  r <- raster::init(r)
  r[] <- 1:raster::ncell(r)

  # plot raster
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
make_raster <- function(coords, buffer = 0, res = NULL){
  # format coords
  coords <- data.frame(coords)
  colnames(coords) <- c("x", "y")

  # get x and y min and max and round up to nearest integer
  # (note: must be an integer for assigning nrow and ncol of a matrix)
  xmin <- ceiling(min(coords$x, na.rm = TRUE) - buffer)
  xmax <- ceiling(max(coords$x, na.rm = TRUE) + buffer)
  ymin <- ceiling(min(coords$y, na.rm = TRUE) - buffer)
  ymax <- ceiling(max(coords$y, na.rm = TRUE) + buffer)

  # make matrix
  m <- matrix(nrow = (ymax - ymin), ncol = (xmax - xmin))

  # turn into raster
  r <- raster::raster(m)

  # set extent
  raster::extent(r) <- c(xmin, xmax, ymin, ymax)

  # set resolution
  if(length(res) > 2) stop("invalid res provided")
  if(!is.null(res)) raster::res(r) <- res


  return(r)

}
