
#' TESS function to do everything
#'
#' @param gen genotype dosage matrix (rows = individuals & columns = snps) or `vcfR` object
#' @param coords dataframe with x and y coordinates
#' @param grid SpatRaster or other gridded spatial object for kriging
#' @param Kvals vector of K values to test
#' @param K_selection how to perform K selection ("manual" to enter into console (default) or "auto" for automatic selection based on \link[algatr]{bestK})
#' @param plot_method method for making rainbow map of kriged layers (options: "maxQ" to only plot the max Q value for each cell (default), "allQ" to plot all Qvalues greater than \code{minQ}, "maxQ_poly" or "allQ_poly" to create the plots as previously described, but as polygons for each K instead of continuous Q values)
#' @param col_breaks number of breaks for plotting (defaults to 20)
#' @param col_alpha
#' @param minQ threshold for minimum Q-value for rainbow plotting if \code{method = "all"} is used (defaults to 0.10)
#' @param tess_method the type of TESS method to be run ("projected.ls" for projected least squares algorithm (default) or "qp" for quadratic programming algorithm)
#' @param ploidy ploidy of data (defaults to 2)
#' @param correct_kriged_Q whether to correct kriged Q values so values greater than 1 are set to 1 and values less than 0 are set to 0 (defaults to TRUE)
#' @inheritParams tess3r::tess3
#'
#' @family TESS functions
#'
#' @return list with all TESS results, final K value, and final kriged raster
#' @export
#'
#' @examples
tess_do_everything <- function(gen, coords, grid, Kvals = 1:10, K_selection = "manual",
                      plot_method = "maxQ", col_breaks = 20, col_alpha = 0.5, minQ = 0.10,
                      tess_method = "projected.ls", ploidy = 2, correct_kriged_Q = TRUE){

  # RUN TESS ---------------------------------------------------------------------------------------------------

  # Convert vcf to dosage
  if(inherits(gen, "vcfR")) gen <- vcf_to_dosage(gen)

  # Convert coords to matrix
  coords <- as.matrix(coords)

  # Test different k values, if more than one provided
  if(length(Kvals) > 1){

    # Run TESS K test
    tess_results <- tess_ktest(gen, coords, Kvals = Kvals, tess_method = tess_method, K_selection = K_selection, ploidy = ploidy)

    # Get K
    K <- tess_results[["K"]]

    # Get tessobj
    tess3_obj <- tess_results[["tess3_obj"]]

    }

  # If only one K value is provided, just use that
  if(length(Kvals) == 1){

    # K is just Kvals if there is only one value
    K <- Kvals

    # Run tess for given K value
    tess3_obj <- tess3r::tess3(X = gen, coord = coords, K = Kvals, method = tess_method, ploidy = ploidy)

  }

  # KRIGE QMATRIX  -----------------------------------------------------------------------------------------------

  # Get Qmatrix
  qmat <- tess3r::qmatrix(tess3_obj, K = K)

  # Give warning if K = 1
  if(K == 1) warning("K = 1, skipping kriging and plotting")

  # Krige Qmatrix
  if(K != 1) krig_admix <- tess_krig(qmat = qmat, coords = coords, grid = grid, correct_kriged_Q = correct_kriged_Q) else krig_admix <- NULL


  # PLOTS --------------------------------------------------------------------------------------------------------

  # Plot max Q-values
  if(K != 1) print(tess_ggplot(krig_admix, coords, plot_method = "maxQ", ggplot_fill = algatr_col_default("ggplot")))

  # Make barplot
  if(K != 1) print(tess_barplot(qmat = qmat, col_pal = algatr_col_default("base")))

  # OUTPUTS ------------------------------------------------------------------------------------------------------

  # Create list with all outputs
  tess_results <- list(K = K,
                       Qmatrix = qmat,
                       krig_admix = krig_admix,
                       tess_results = tess3_obj,
                       coords = coords,
                       Kvals = Kvals,
                       grid = grid)

  return(tess_results)
}


#' Test multiple K values
#'
#' @inheritParams tess_do_everything
#' @return
#' @export
#'
#' @family TESS functions
#'
#' @examples
tess_ktest <- function(gen, coords, Kvals = 1:10, grid = NULL, tess_method = "projected.ls", K_selection = "manual", ploidy = 2){

  # Format coordinates
  coords <- as.matrix(coords)

  # Run tess for all K values
  tess3_obj <- tess3r::tess3(X = gen, coord = coords, K = Kvals, method = tess_method, ploidy = ploidy)

  # Plot CV results
  plot(tess3_obj, pch = 19, col = "blue",
       xlab = "Number of ancestral populations",
       ylab = "Cross-validation score")

  # Get best K value
  if(K_selection == "auto"){K <-  bestK(tess3_obj, Kvals)}
  if(K_selection == "manual"){K <- as.numeric(readline(prompt = "Enter K Value: "))}

  # Mark the K-value selected
  abline(v = K, col = "red", lty = "dashed")

  # Create list with tess3 object and K value
  tess_results <- list(K = K,
                    tess3_obj = tess3_obj,
                    coords = coords,
                    Kvals = Kvals,
                    grid = grid)

  return(tess_results)
}


#' Krige admixture values
#'
#' @param qmat qmatrix
#' @inheritParams tess_do_everything
#'
#' @return
#' @export
#'
#' @family TESS functions
#'
#' @examples
tess_krig <- function(qmat, coords, grid, correct_kriged_Q = TRUE){

  # Convert RasterLayer/Stack/Brick to SpatRaster
  if(!inherits(grid, "SpatRaster")) grid <- terra::rast(grid)

  # Get K
  K <- ncol(qmat)

  # Make grid for kriging
  if (inherits(grid, "SpatRaster")) {
    krig_grid <- raster_to_grid(grid)
  } else if (sp::gridded(grid)) {
    krig_grid <- grid
  } else {
    stop(" unable to find an inherited method for type of grid provided")
  }

  # Remove any CRS values before kriging (autoKrige doesn't like lonlat projection systems)
  terra::crs(krig_grid) <- NA

  # Make coords into spatial object
  krig_df <- data.frame(coords)
  sp::coordinates(krig_df) <- ~x+y

  # Krige each K value
  # krig_admix <- raster::stack(purrr::map(1:K, krig_K, qmat, krig_grid, krig_df))
  krig_admix_terra = terra::rast(purrr::map(1:K, krig_K, qmat, krig_grid, krig_df))

  krig_admix <- purrr::map(1:K, krig_K, qmat, krig_grid, krig_df)

  # If NULL, return NULL
  if(is.null(krig_admix)) return(NULL)

  # Convert all values in raster greater than 1 to 1 and all values less than 0 to 0
  if(correct_kriged_Q){
    krig_admix[krig_admix < 0] <- 0
    krig_admix[krig_admix > 1] <- 1
  }

  # Rename layers
  names(krig_admix) <- paste0("K", 1:K)

  # Return stack
  return(krig_admix)
}

#' Krige one K value
#'
#' @param K K value
#' @param qmat Q matrix
#'
#' @export
#' @noRd
#' @family TESS functions
#'
krig_K <- function(K, qmat, krig_grid, krig_df){

  # Add Q values to spatial dataframe
  krig_df$Q <- qmat[,K]

  # Skip if all of the Q values are identical (kriging not possible)
  if(length(unique(krig_df$Q)) == 1){
    warning(paste0("Only one unique Q value for K = ", K, ", returning NA raster (note: may want to consider different K value)")); NULL
    blank <- terra::rast(krig_grid, type = "xyz")
    blank[] <- NA
    return(blank)
  }

  # Krige (capture output so it is not printed automatically)
  co <- capture.output(krig_res <- automap::autoKrige(Q ~ 1, krig_df, krig_grid))

  # Get Krige output
  krig_spdf <- krig_res$krige_output

  # Turn into raster
  krig_raster_terra <- terra::rast(krig_spdf, type = "xyz")
  # krig_raster <- raster::rasterFromXYZ(krig_spdf)

  return(krig_raster)
}


#' Convert a raster to a grid
#'
#' @param x SpatRaster
#'
#' @return gridded SpatialPixelsDataFrame
#'
#' @family TESS functions
#'
#' @export
#' @noRd
raster_to_grid <- function(x) {

  # Convert raster to dataframe
  grd <- terra::as.data.frame(x, xy = TRUE, na.rm = FALSE)

  # Convert dataframe to spatial dataframe
  sp::coordinates(grd) <- ~ x + y

  # Convert into gridded object
  sp::gridded(grd) <- TRUE

  return(grd)
}

#' ggplot of TESS results
#'
#' @param krig_admix SpatRaster returned by \link[algatr]{tess_krig}
#' @param coords dataframe with x and y coordinates for plotting (optional)
#' @param plot_method method for making rainbow map of kriged layers (options: "maxQ" to only plot the max Q value for each cell (default), "allQ" to plot all Qvalues greater than \code{minQ}, "maxQ_poly" or "allQ_poly" to create the plots as previously described, but as polygons for each K instead of continuous Q values)
#' @param ggplot_fill any ggplot2 scale fill discrete function (default: \link[algatr]{scale_fill_viridis_d}, \code{option = "turbo"})
#' @param minQ threshold for minimum Q-value for rainbow plotting if \code{method = "all"} is used (defaults to 0.10)
#' @inheritParams tess_do_everything
#'
#' @family TESS functions
#'
#' @return
#' @export
#'
#' @examples
tess_ggplot <- function(krig_admix, coords = NULL, plot_method = "maxQ", ggplot_fill = algatr_col_default("ggplot"), minQ = 0.10, plot_axes = FALSE){

  # Set up ggplot df
  gg_df <- krig_admix %>%
    terra::as.data.frame(x, xy = TRUE, na.rm = FALSE) %>%
    tidyr::as_tibble() %>%
    tidyr::gather("K", "Q", -c(x, y)) %>%
    dplyr::mutate(K = as.factor(gsub("K", "", K))) %>%
    dplyr::group_by(x, y)

  # Use max or all Q
  if(plot_method == "maxQ" | plot_method == "maxQ_poly") gg_df <- gg_df %>% dplyr::top_n(1, Q)
  if(plot_method == "allQ" | plot_method == "allQ_poly") gg_df <- gg_df %>% dplyr::filter(Q > 0.20)

  # Set up base plot
  plt <- ggplot2::ggplot()

  # Plot as polygon or continuous Q
  if(plot_method == "maxQ_poly" | plot_method == "allQ_poly"){
    plt <- plt + ggplot2::geom_tile(data = gg_df, ggplot2::aes(x = x, y = y, fill = K), alpha = 0.5)
  } else {
    plt <- plt +
      ggplot2::geom_tile(data = gg_df, ggplot2::aes(x = x, y = y, fill = K, alpha = Q)) +
      ggplot2::scale_alpha_binned(breaks = round(seq(0, 1, by = 0.10), 1),
                                guide = ggplot2::guide_legend())
  }

  # Add color
  plt <- plt + ggplot_fill

  # Add themes and coord controls
  plt <-  plt + ggplot2::coord_equal() + ggplot2::theme_bw()

  # Add axes
  if(plot_axes) plt <- plt + ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(),
                                            panel.grid.major.y = ggplot2::element_blank(),
                                            panel.grid.minor.x = ggplot2::element_blank(),
                                            panel.grid.major.x = ggplot2::element_blank(),
                                            aspect.ratio = 1)

  if(!plot_axes) plt <- plt + ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank(),
                                            panel.grid.major.y = ggplot2::element_blank(),
                                            panel.grid.minor.x = ggplot2::element_blank(),
                                            panel.grid.major.x = ggplot2::element_blank(),
                                            axis.title.x = ggplot2::element_blank(),
                                            axis.text.x = ggplot2::element_blank(),
                                            axis.ticks.x = ggplot2::element_blank(),
                                            axis.title.y = ggplot2::element_blank(),
                                            axis.text.y = ggplot2::element_blank(),
                                            axis.ticks.y = ggplot2::element_blank(),
                                            panel.border = ggplot2::element_blank(),
                                            aspect.ratio = 1)

  # Add coords
  if(!is.null(coords)) plt <- plt + ggplot2::geom_point(data = data.frame(coords), ggplot2::aes(x = x, y = y))

  return(plt)
}

#' Make rainbow TESS plot from kriged admixture rasters
#'
#' @param krig_admix SpatRaster returned by \link[algatr]{tess_krig}
#' @param coords dataframe with x and y coordinates for plotting (optional)
#' @param plot_method method for making rainbow map of kriged layers (options: \code{"maxQ"} to only plot the max Q value for each cell (default), \code{"allQ"} to plot all Qvalues greater than \code{minQ}, \code{"maxQ_poly"} or \code{"allQ_poly"} to create the plots as previously described, but as polygons for each K instead of continuous Q values)
#' @param col_pal function that creates a vector of contiguous colors (defaults to \link[viridis]{turbo}). If using \code{plot_method = "allQ"} or \code{plot_method = "allQ_poly"}, the function must accept \code{alpha} as an argument
#' @param col_breaks if using \code{"maxQ"} and \code{"allQ"} plot methods, the number of breaks to use when plotting kriged maps
#' @param col_alpha if using the \code{"allQ"} plot method, an alpha-transparency level in the range [0,1] (0 means transparent and 1 means opaque) provided to \code{col_pal} function
#' @param minQ threshold for minimum Q-value for rainbow plotting if \code{plot_method = "allQ"} is used (defaults to 0.10)
#' @param legend whether to include legend (defaults to TRUE)
#' @inheritParams tess_do_everything
#'
#' @return
#' @export
#' @family TESS functions
#'
#' @examples
tess_plot <- function(krig_admix, coords = NULL, plot_method = "maxQ", col_pal = algatr_col_default("base"), col_breaks = 20, col_alpha = 0.50, minQ = 0.10, legend = TRUE){

  # Get K based on the number of layers
  K <- terra::nlyr(krig_admix)

  # Select method and options and suppress irrelevant plot warnings
  suppressWarnings({
  if(plot_method == "maxQ") tess_plot_max(krig_admix, K = K, coords = coords, poly = FALSE, col_pal = col_pal, col_breaks = col_breaks, legend = TRUE)
  if(plot_method == "allQ") tess_plot_all(krig_admix, K = K, coords = coords, poly = FALSE, col_pal = col_pal, col_breaks = col_breaks, col_alpha = col_alpha, minQ = minQ, legend = TRUE)
  if(plot_method == "maxQ_poly") tess_plot_max(krig_admix, K = K, coords = coords, poly = TRUE, col_pal = col_pal, legend = TRUE)
  if(plot_method == "allQ_poly") tess_plot_all(krig_admix, K = K, coords = coords, poly = TRUE, col_pal = col_pal, col_alpha = col_alpha, minQ = minQ, legend = TRUE)
  })

  # Add coordinates, if provided
  if(!is.null(coords)) points(coords, pch = 3)
}


#' Plot method: Q max
#'
#' @inheritParams tess_plot
#' @param K K value
#' @param poly whether to plot as polygon instead of continuous Q values
#'
#' @export
#' @family TESS functions
#' @noRd
tess_plot_max <- function(krig_admix, K, coords = NULL, poly = FALSE, col_pal = algatr_col_default("base"), col_breaks = 20, legend = TRUE){

  # Make and summarize dataframe by only retaining highest Q values for each point
  pop_df <-  krig_admix %>%
    terra::as.data.frame(x, xy = TRUE, na.rm = FALSE) %>%
    tidyr::as_tibble() %>%
    tidyr::gather("K", "Q", -c(x, y)) %>%
    dplyr::mutate(K = as.numeric(gsub("K", "", K))) %>%
    dplyr::group_by(x, y) %>%
    dplyr::top_n(1, Q)

  # Get extent of raster (to set up plot)
  ext <- extent(krig_admix)

  # Set up base plot (! important ! Don't remove or things will get wonky)
  plot(1, legend = FALSE, axes = FALSE, box = FALSE, type = "n",
       xlim = c(ext[1], ext[2]), ylim = c(ext[3], ext[4]),
       xlab = "", ylab = "")

  # Plot each kriged admixture layer one by one on top of each other
  purrr::walk(1:K, max_plot_helper, pop_df, poly = poly, col = col_pal(K), col_breaks = col_breaks)

  # Add legend
  if(legend) legend("topright", pch = 15, legend = paste0("K = ", 1:K), col = col_pal(K), bty = "n")

  # Add coordinates
  if(!is.null(coords)) points(coords, pch = 3)

}

#' K max plotting helper function
#'
#' @inheritParams tess_plot
#' @param K K value
#' @param pop_df SpatialPointsDataFrame with K and Q-values
#' @param poly whether to plot as polygon instead of continuous Q values
#' @param col single color code
#'
#' @export
#' @noRd
#' @family TESS functions
#'
max_plot_helper <- function(K, pop_df, poly, col, col_breaks = 20, zlim = NULL){

  # Subset by K
  pop_spdf <- pop_df[pop_df$K == K, ]

  # Skip to next iteration if there is no more than one value for that K
  if(nrow(pop_spdf) < 2){
    warning(paste0("less than two values found for K = ", K,", skipping..."))
    next
  }

  # Make into spdf and convert to raster
  sp::coordinates(pop_spdf) <- ~x+y
  rl <- terra::rast(pop_spdf[, "Q"], type="xyz")

  # If not poly plot, set zlim to range of all Q values
  if(!poly) zlim <- range(pop_df$Q)

  # Plot raster
  terra::plot(rl,
               add = TRUE,
               legend = FALSE,
               col = make_plot_col(K, col, col_breaks, poly),
               zlim = zlim)
}

#' Plot method: Q all
#'
#' @inheritParams tess_plot
#' @param K K value
#' @param poly whether to plot as polygon instead of continous Q values
#'
#' @export
#' @noRd
#' @family TESS functions
#'
tess_plot_all <- function(krig_admix, K = K, coords = NULL, poly = FALSE, col_pal = algatr_col_default("base"), col_breaks = 20, col_alpha = 0.50, minQ = 0.10, legend = TRUE){

  # Get max raster value for plotting (minQ defines minimum)
  maxr <- max(maxValue(krig_admix))

  # Get extent of raster (to set up plot)
  ext <- extent(krig_admix)

  # Set up base plot (! important ! Don't remove or things will get wonky)
  plot(1, legend = FALSE, axes = FALSE, box = FALSE, type = "n",
       xlim = c(ext[1], ext[2]), ylim = c(ext[3], ext[4]),
       xlab = "", ylab = "")

  # Plot kriged admixture layers while masking values < minQ
  purrr::walk(1:K, all_plot_helper, krig_admix, poly = poly, col = col_pal(K, alpha = col_alpha), col_breaks = col_breaks, zlim = c(minQ, maxr))

  # Add legend
  if(legend) legend("topright", pch = 15, legend = paste0("K = ", 1:K), col = col_pal(K), bty = "n")

  # Add coordinates
  if(!is.null(coords)) points(coords, pch = 3)

}

#' Helper function for max plotting
#'
#' @inheritParams tess_plot
#' @param K K value
#' @param poly whether to plot as polygon instead of continous Q values
#' @param col single color code
#'
#' @export
#' @noRd
#' @family TESS functions
#'
all_plot_helper <- function(K, krig_admix, poly, col, col_breaks = 20, zlim = NULL){

  # Plot raster
  terra::plot(krig_admix[[K]],
               col = make_plot_col(K, col, col_breaks, poly),
               add = TRUE,
               legend = FALSE,
               zlim = zlim)

}

#' Plot all kriged Q values for each K
#'
#' @param krig_admix RasterStack returned by \link[algatr]{tess_krig}
#' @param coords dataframe with x and y coordinates for plotting (optional)
#' @param ... Graphical parameters. Any argument that can be passed to image.plot and to base plot
#' @inheritParams tess_do_everything
#'
#' @return
#' @export
#' @family TESS functions
#'
#' @examples
tess_plot_allK <- function(krig_admix, coords = NULL, col_pal = algatr_col_default("base"), col_breaks = 20, ...){

  # Get K
  K <- terra::nlyr(krig_admix)

  # Plot kriged admixture maps while masking small values (e.g. < minQ)
  purrr::walk(1:K, allK_plot_helper, krig_admix, coords = coords,  col = col_pal(K), col_breaks = col_breaks, ...)

  }

#' Helper function for all K plotting
#'
#' @inheritParams tess_plot_allK
#' @param K K value
#' @param col single color code
#'
#' @export
#' @family TESS functions
#'
allK_plot_helper <- function(K, krig_admix, coords = NULL, col, col_breaks, ...){

  # Suppress irrelevant plot warnings
  suppressWarnings({terra::plot(krig_admix[[K]],
               col = make_plot_col(K, col, col_breaks, alpha = 1, start_col = rgb(0.94, 0.94, 0.95, 1)),
               zlim = c(0, max(maxValue(krig_admix))),
               axes = FALSE,
               box = FALSE,
               main = paste0("K = ",K),
               ...)})

  # Add coordinates, if provided
  if(!is.null(coords)) points(coords, pch = 3)
}

#' Make color vector for plotting
#'
#' @param K K value (used to index col)
#' @param col vector of colors
#' @param alpha transparency to start color scale at
#' @inheritParams tess_plot
#'
#' @export
#' @noRd
#' @family TESS functions
#'
make_plot_col <- function(K, col, col_breaks, poly = FALSE, alpha = 0, start_col = rgb(1, 1, 1, alpha)){

  if(poly){
    # Make color palette using only solid color
    plot_col <- col[K]
  } else {
    # Make color palette gradient from transparent to color defined by K
    kpal <- colorRampPalette(c(start_col, col[K]), interpolate = "linear", alpha = TRUE)
    plot_col <- kpal(col_breaks)
  }

  return(plot_col)
}

#' Create TESS Barplot
#'
#' Based on code from: https://github.com/bcm-uga/TESS3_encho_sen/blob/master/R/plotQ.R
#'
#' @param qmat Q matrix
#' @param sort_by_Q whether to sort bars by Q value (equivalent to \link[tess3r]{barplot} sort.by.Q)
#' @param legend whether to display legend (defaults to TRUE)
#' @param legend_position the x and y coordinates or keyword to determine legend position (defaults to bottom right)
#' @inheritParams tess_do_everything
#' @inheritParams graphics::barplot
#' @param ... other parameters of the function \code{\link{barplot.default}}.
#'
#' @return
#' @export
#'
#' @family TESS functions
#' @examples
tess_barplot <- function(qmat, col_pal = algatr_col_default("base"), sort_by_Q = TRUE, legend = TRUE, legend_position = "bottomright", border = NA, space = 0, ...){
  # CODE ADAPTED FROM: https://github.com/bcm-uga/TESS3_encho_sen/blob/master/R/plotQ.R

  # Get K
  K <- ncol(qmat)

  if (sort_by_Q) {
    gr = apply(qmat, MARGIN = 1, which.max)
    gm = max(gr)
    gr.o = order(sapply(1:gm, FUN = function(g) mean(qmat[,g])))
    gr = sapply(gr, FUN = function(i) gr.o[i])
    or = order(gr)
    Qm = t(qmat[or,])
    class(Qm) = "matrix"
    graphics::barplot(Qm, col =  col_pal(K), border = border, space = space, ...)
    legend("bottomright", pch = 15, legend = paste0("K = ", 1:K), col = col_pal(K))
    return(list(order = or))
  }
  else {
    Qm = t(qmat)
    class(Qm) = "matrix"
    graphics::barplot(Qm, col =  col_pal(ncol(qmat)), border = border, space = space, ...)
    legend(legend_position, pch = 15, legend = paste0("K = ", 1:K), col = col_pal(K))
    return(list(order = 1:nrow(qmat)))
  }

}


#' Best K Selection based on cross entropy
#'
#' @param tess3_obj list produced by \code{\link{tess3}}
#' @param Kvals vector of K values for testing
#'
#' @note (source: https://chazhyseni.github.io/NALgen/post/determining_bestk/)
#' @return
#' @export
#' @family TESS functions
#'
#' @examples
bestK <- function(tess3_obj, Kvals){
  ce <- list()
  for(k in Kvals) ce[[k]] <- tess3_obj[[k]]$crossentropy
  ce.K <- c()
  for(k in Kvals) ce.K[k] <- min(ce[[k]])
  diff <- ce.K[-1] - ce.K[-max(Kvals)]
  slope <- exp(-diff) - 1
  # K is selected based on the smallest slope value in the upper quartile
  K <- min(which(slope <= quantile(slope)[4]))
  return(K)
}


#' Create default TESS color palette
#' @param n number of colors to generate (must be less than 9)
#'
#' @export
#' @family TESS functions
tess_col_default <- function(n, alpha = 1){
  if (n > 9) stop("The default color palette expects less than 9 values")
  tessCP <- CreatePalette()
  tesscol <- sapply(1:n, function(x, tessCP){tessCP[[x]][9]}, tessCP)
  tesspal <- colorRampPalette(tesscol, interpolate = "linear", alpha = TRUE)
  return(tesspal(n))
}

#' Create default algatr color palette for TESS
#'
#' @param x whether to return ggplot or base color scale function
#'
#' @export
#' @noRd
#' @family TESS functions
algatr_col_default <- function(x){
  if(x == "ggplot") col <- ggplot2::scale_fill_viridis_d(option = "turbo", begin = 0.1, end = 0.9)
  if(x == "base") col <- function (n, alpha = 1, begin = 0, end = 1, direction = 1) viridis::viridis(n, alpha, begin = 0.1, end = 0.9, direction, option = "turbo")
  return(col)
}
