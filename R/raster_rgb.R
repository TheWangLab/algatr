#' RGB plot of a raster
#'
#' Wrapper of the \link[terra]{plotRGB} function with added steps to perform
#' re-scaling of raster layers (e.g., from 0 to 255 for color layers).
#'
#' @param x RasterStack or SpatRaster
#' @inheritParams terra::plotRGB
#'
#' @export
plot_rgb <- function(x, r=1, g=2, b=3, a=NULL, scale=NULL, mar=0,
                     stretch=NULL, smooth=TRUE, colNA="white", alpha=NULL, bgalpha=NULL,
                     zlim=NULL, zcol=FALSE, axes=FALSE, ...){
  if (!inherits(x, "SpatRaster")) x <- terra::rast(x)

  # Scale rasters to get colors (each layer will correspond with R, G, or B in the final plot)
  pcaRastRGB <- stack_to_rgb(x)

  # If there is an alpha layer, rescale that layer from 0 to 1
  if (!is.null(a)) {
    mm <- terra::minmax(x[[a]])
    pcaRastRGB[[a]] <- (x[[a]] - mm[1,]) / (mm[2,] - mm[1,])
  }

  # Plot RGB
  terra::plotRGB(pcaRastRGB, r=r, g=g, b=b, a=a, scale=scale, mar=mar,
                 stretch=stretch, smooth=smooth, colNA=colNA, alpha=alpha, bgalpha=bgalpha,
                 zlim=zlim, zcol=zcol, axes=axes, ...)
}

#' Scale a raster stack from 0 to 255
#'
#' @param s RasterStack
#'
#' @noRd
#' @export
stack_to_rgb <- function(s) {
  stack_list <- as.list(s)
  new_stack <- terra::rast(purrr::map(stack_list, raster_to_rgb))
  return(new_stack)
}
