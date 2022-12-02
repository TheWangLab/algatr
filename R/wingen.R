
#' wingen function to do everything  (preview and generate moving window maps, krige, and mask)
#'
#' @param preview whether to produce preview of raster layer, window and focal cell size using \link[wingen]{preview_gd} (default = FALSE)
#' @param kriged whether to smooth out mapped values using kriging using \link[wingen]{krig_gd} (default = FALSE)
#' @param masked whether to mask out areas outside region of interest using \link[wingen]{mask_gd} (default = FALSE)
#' @param plot_count if TRUE, whether to visualize sample counts using \link[wingen]{plot_count} (default = FALSE)
#' @inheritParams wingen::preview_gd
#' @inheritParams wingen::window_gd
#' @inheritParams wingen::krig_gd
#' @inheritParams wingen::mask_gd
#' @inheritParams wingen::plot_count
#'
#' @return RasterBrick object final raster
#' @family wingen function
#'
#' @details
#' When using wingen, please cite the original citation: Bishop, A.P., Chambers, E.A., Wang, I.J. (202X) TODO: citation
#' N.B.: Be aware that this function sets *many* of the wingen function arguments to defaults, which may result in sub optimal results. We highly advise researchers to run each wingen function separately for best results.
#'
wingen_do_everything <- function(preview = FALSE, lyr, coords, wdim = 3, fact = 0, sample_count = TRUE, min_n = 2,
                                 vcf, stat, rarify = FALSE,
                                 kriged = FALSE, grd = NULL, index = 1, agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL,
                                 masked = FALSE, mask = NULL, bkg = NULL, plot_count = FALSE){

  if(preview == TRUE){
    if (fact == 0) lyr <- lyr * 0 else lyr <- raster::aggregate(lyr, fact, fun = mean) * 0
    if (ncell(lyr) > 10000) warning("The number of cells exceeds 10,000; you may want to increase the aggregation factor using the `fact` argument to decrease computational time!")
    print(wingen::preview_gd(lyr = lyr, coords = coords, wdim = wdim, fact = fact, sample_count = sample_count, min_n = min_n))
    input <- utils::menu(c("Y", "N"), title="Would you like to continue running wingen with these parameters?")
    if (input == 1) {print("OK, running wingen...")}
    if (input == 2) {stop("Stopping the run")}
  }

  if (fact == 0) lyr <- lyr * 0 else lyr <- raster::aggregate(lyr, fact, fun = mean) * 0
  if (ncell(lyr) > 10000) warning("The number of cells exceeds 10,000; you may want to increase the aggregation factor using the `fact` argument to decrease computational time!")
  map <- wingen::window_gd(vcf = vcf, coords = coords, lyr = lyr, stat = stat,
                           wdim = wdim, fact = fact, rarify = rarify)


  # KRIGING -----------------------------------------------------------------

  if (kriged == TRUE) map <- krig_helper(map, grd = grd, index = index, agg_grd = agg_grd, disagg_grd = disagg_grd, agg_r = agg_r, disagg_r = disagg_r)

  # MASKING -----------------------------------------------------------------

  if(masked == TRUE){map <- wingen::mask_gd(x = map, mask = mask)}

  # RESULTS -----------------------------------------------------------------
  # Plot genetic diversity
  print(wingen::plot_gd(map, bkg = bkg, index = index))

  # Plot sample counts
  if(plot_count == TRUE) print(wingen::plot_count(map))

  # Return object with results
  return(map)
}

#' Helper function for kriging wingen moving window; checks raster resolution and runs kriging
#'
#' @param map RasterLayer or RasterStack produced by `window_gd()`
#' @param grd object to create grid for kriging; inherits from \link[wingen]{krig_gd} (defaults to NULL)
#' @param index integer indices of layers in raster stack to krige; inherits from \link[wingen]{krig_gd} (defaults to 1)
#' @param agg_grd factor to use for aggregation of `grd`; inherits from \link[wingen]{krig_gd} (defaults to NULL)
#' @param disagg_grd factor to use for disaggregation of `grd`; inherits from \link[wingen]{krig_gd} (defaults to NULL)
#' @param agg_r factor to use for aggregation of `r`; inherits from \link[wingen]{krig_gd} (defaults to NULL)
#' @param disagg_r factor to use for disaggregation of `r`; inherits from \link[wingen]{krig_gd} (defaults to NULL)
#'
#' @return kriged map
#' @export
#' @family wingen functions
#'
#' @examples
krig_helper <- function(map, grd = NULL, index = 1, agg_grd = NULL, disagg_grd = NULL, agg_r = NULL, disagg_r = NULL){
  # Perform checks ----------------------------------------------------------
  if (!is.null(agg_grd)) grd <- krig_agg_helper(to_krig = grd, agg_disagg = agg_grd, agg_spec = "agg")
  if (!is.null(disagg_grd)) grd <- krig_agg_helper(to_krig = grd, agg_disagg = disagg_grd, agg_spec = "disagg")
  if (!is.null(agg_r)) grd <- krig_agg_helper(to_krig = r, agg_disagg = agg_r, agg_spec = "agg")
  if (!is.null(disagg_r)) grd <- krig_agg_helper(to_krig = r, agg_disagg = disagg_r, agg_spec = "disagg")

  # Warning if too many cells in agg/disagg raster
  if (raster::ncell(grd) > 10000) {
    warning("The resolution of your kriging raster layer is very high and no aggregation is being performed; you should perform aggregation to reduce computational time!")
  }

  map <- wingen::krig_gd(map, grd = grd, index = index)

  return(map)
}

#' Helper function for krig_helper; calculates aggregated/disaggregated raster cell size
#'
#' @param to_krig object to create grid for kriging, `grd` or `r` from \link[wingen]{krig_gd}
#' @param agg_disagg aggregation or disaggregation parameter, one of agg_grd, disagg_grd, agg_r, or disagg_r from \link[wingen]{krig_gd}
#' @param agg_spec whether aggregation or disaggregation is performed, options are `"disagg"` or `"agg"` (defaults to `"agg"`)
#'
#' @return number of cells contained in final (aggregated or disaggregated) raster layer
#' @export
#' @family wingen function
#'
#' @examples
krig_agg_helper <- function(to_krig, agg_disagg, agg_spec = "agg"){
  if (agg_spec == "agg") raster::aggregate(to_krig, agg_disagg)
  if (agg_spec == "disagg") raster::disaggregate(to_krig, agg_disagg)
}
