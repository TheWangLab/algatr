
#' TESS function to do everything
#'
#' @param gen genotype dosage matrix (rows = individuals & columns = snps) or `vcfR` object
#' @param coords dataframe with x and y coordinates
#' @param grid SpatRaster for kriging
#' @param Kvals vector of K values to test
#' @param K_selection how to perform K selection ("manual" to enter into console (default) or "auto" for automatic selection based on \link[algatr]{bestK})
#' @param plot_method method for making rainbow map of kriged layers (options: "maxQ" to only plot the max Q value for each cell (default), "allQ" to plot all Q values greater than \code{minQ}, "maxQ_poly" or "allQ_poly" to create the plots as previously described, but as polygons for each K instead of continuous Q values)
#' @param col_breaks number of breaks for plotting (defaults to 20)
#' @param minQ threshold for minimum Q value for rainbow plotting if \code{method = "all"} is used (defaults to 0.10)
#' @param tess_method the type of TESS method to be run ("projected.ls" for projected least squares algorithm (default) or "qp" for quadratic programming algorithm)
#' @param lambda numeric value for the spatial regularization parameter. The default value lambda = 1 attributes equal weights to the loss function and to the penalty function.
#' @param ploidy ploidy of data (defaults to 2)
#' @param correct_kriged_Q whether to correct kriged Q values so values greater than 1 are set to 1 and values less than 0 are set to 0 (defaults to TRUE)
#' @param quiet whether to print output tables and figures (defaults to FALSE)
#'
#' @family TESS functions
#'
#' @details
#' TESS is run using the tess3r package: Caye, K., François, O. (2016). tess3r: Inference of Spatial Population Genetic Structure. R package version 1.1.0.
#' See also: Caye, K., Deist, T.M., Martins, H., Michel, O., François, O. (2016). TESS3: fast inference of spatial population structure and genome scans for selection. Mol. Ecol. Res. 16(2):540-548. https://doi.org/10.1111/1755-0998.12471
#'
#' @return list with all TESS results, final K value, and final kriged raster
#' @export
tess_do_everything <- function(gen, coords, grid = NULL, Kvals = 1:10, K_selection = "manual",
                               plot_method = "maxQ", col_breaks = 20, minQ = 0.10,
                               tess_method = "projected.ls", lambda = 1, ploidy = 2, correct_kriged_Q = TRUE,
                               quiet = FALSE) {
  message("Please be aware: the do_everything functions are meant to be exploratory. We do not recommend their use for final analyses unless certain they are properly parameterized.")

  # RUN TESS ---------------------------------------------------------------------------------------------------

  # Convert vcf to dosage
  if (inherits(gen, "vcfR")) gen <- vcf_to_dosage(gen)

  # Convert sf coords to matrix
  if (inherits(coords, "sf")) {
    coords <- sf::st_coordinates(coords)
  }

  # Convert coords to matrix
  coords <- as.matrix(coords)

  # Test different k values, if more than one provided
  if (length(Kvals) > 1) {
    # Run TESS K test
    tess_results <- tess_ktest(gen, coords, Kvals = Kvals, tess_method = tess_method, lambda = lambda, K_selection = K_selection, ploidy = ploidy, quiet = quiet)

    # Get K
    K <- tess_results[["K"]]

    # Get tessobj
    tess3_obj <- tess_results[["tess3_obj"]]

    # Get population assignments
    pops <- tess_results[["pops"]]
  }

  # If only one K value is provided, just use that
  if (length(Kvals) == 1) {
    # K is just Kvals if there is only one value
    K <- Kvals

    # Run tess for given K value
    tess3_obj <- tess3r::tess3(X = gen, coord = coords, K = Kvals, method = tess_method, lambda = lambda, ploidy = ploidy)

    # Get population assignments
    pops <- pops_helper(gen = gen, tess3_obj = tess3_obj, K = Kvals)
  }

  # KRIGE QMATRIX  -----------------------------------------------------------------------------------------------

  # Get Qmatrix
  qmat <- tess3r::qmatrix(tess3_obj, K = K)

  # Give warning if K = 1
  if (K == 1) warning("K = 1, skipping kriging and plotting")

  # Give warning if grid is not provided
  if (is.null(grid)) warning("Grid not provided, skipping kriging")

  # Krige Qmatrix
  if (K != 1 & !is.null(grid)) krig_admix <- tess_krig(qmat = qmat, coords = coords, grid = grid, correct_kriged_Q = correct_kriged_Q) else krig_admix <- NULL

  # PLOTS --------------------------------------------------------------------------------------------------------

  # Plot Q-values
  if (!quiet) {
    # Make map
    if (K != 1 & !is.null(grid)) print(tess_ggplot(krig_admix, coords, plot_method = plot_method, ggplot_fill = algatr_col_default("ggplot")))

    # Make barplot
    if (K != 1) print(tess_ggbarplot(qmat = qmat, ggplot_fill = algatr_col_default("ggplot")))
  }

  # OUTPUTS ------------------------------------------------------------------------------------------------------

  # Create list with all outputs
  tess_results <- list(
    K = K,
    Qmatrix = qmat,
    krig_admix = krig_admix,
    tess_results = tess3_obj,
    coords = coords,
    Kvals = Kvals,
    grid = grid,
    pops = pops
  )

  return(tess_results)
}

#' Test multiple K values
#'
#' @inheritParams tess_do_everything
#' @return list with results from testing different K-values
#' @export
#'
#' @family TESS functions
tess_ktest <- function(gen, coords, Kvals = 1:10, grid = NULL, tess_method = "projected.ls", lambda = 1, K_selection = "manual", ploidy = 2, quiet = FALSE) {
  # Format coordinates
  coords <- as.matrix(coords)

  # Run tess for all K values
  tess3_obj <- tess3r::tess3(X = gen, coord = coords, K = Kvals, method = tess_method, lambda = lambda, ploidy = ploidy)

  # Plot CV results
  if (!quiet) {
    plot(tess3_obj,
      pch = 19, col = "blue",
      xlab = "Number of ancestral populations",
      ylab = "Cross-validation score",
      xaxt = "n"
    )
    axis(side = 1, at = Kvals)
  }

  # Get best K value
  if (K_selection == "auto") {
    K <- bestK(tess3_obj, Kvals)
  }
  if (K_selection == "manual") {
    K <- as.numeric(readline(prompt = "Enter K Value: "))
  }

  # Mark the K-value selected
  if (!quiet) abline(v = K, col = "red", lty = "dashed")

  # Get population assignments
  pops <- pops_helper(gen = gen, tess3_obj = tess3_obj, K = K)

  # Create list with tess3 object and K value
  tess_results <- list(
    K = K,
    tess3_obj = tess3_obj,
    coords = coords,
    Kvals = Kvals,
    grid = grid,
    pops = pops
  )

  return(tess_results)
}

#' Krige admixture values
#'
#' @param qmat qmatrix
#' @inheritParams tess_do_everything
#'
#' @return Raster\* type object of kriged Q values
#' @export
#'
#' @family TESS functions
tess_krig <- function(qmat, coords, grid = NULL, correct_kriged_Q = TRUE) {
  # Check CRS
  crs_check(coords, grid)

  # Get K
  K <- ncol(qmat)

  # Make grid for kriging
  if (!inherits(grid, "SpatRaster")) grid <- terra::rast(grid)
  krig_grid <- raster_to_grid(grid)

  # Convert coords
  krig_df <- coords_to_sp(coords)

  # Krige each K value
  krig_admix <-
    purrr::map(1:K, krig_K, qmat, krig_grid, krig_df) %>%
    terra::rast()

  # mask with original raster layer because the grid fills in all NAs
  #( note: we don't remove NAs because it can change the extent)
  grid <- terra::resample(grid, krig_admix[[1]])
  krig_admix <- terra::mask(krig_admix, grid)

  # Convert all values in raster greater than 1 to 1 and all values less than 0 to 0
  if (correct_kriged_Q) {
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
krig_K <- function(K, qmat, krig_grid, krig_df) {
  # Add Q values to spatial dataframe
  krig_df$Q <- qmat[, K]

  # Skip if all of the Q values are identical (kriging not possible)
  if (length(unique(krig_df$Q)) == 1) {
    warning(paste0("Only one unique Q value for K = ", K, ", returning NULL (note: may want to consider different K value)"))
    return(NULL)
  }

  # Krige (capture output so it is not printed automatically)
  co <- capture.output(krig_res <- automap::autoKrige(Q ~ 1, krig_df, new_data = krig_grid))

  # Get Krige output
  krig_spdf <- krig_res$krige_output

  # turn spdf into raster
  krig_r <- terra::rast(krig_spdf, type = "xyz", crs = terra::crs(krig_grid))
  # return just the prediction (may want to provide var/stdev in the future)
  krig_r <- krig_r[[1]]

  return(krig_r)
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
#' @param plot_method method for making rainbow map of kriged layers (options: "maxQ" to only plot the max Q value for each cell (default), "allQ" to plot all Q values greater than \code{minQ}, "maxQ_poly" or "allQ_poly" to create the plots as previously described, but as polygons for each K instead of continuous Q values)
#' @param ggplot_fill any ggplot2 scale fill discrete function (default: \link[algatr]{scale_fill_viridis_d}, \code{option = "turbo"})
#' @param minQ threshold for minimum Q-value for rainbow plotting if \code{plot_method = "allQ"} or \code{plot_method = "allQ_poly"} is used (defaults to 0.10)
#' @param plot_axes whether to plot axes or not (defaults to FALSE)
#' @param rel_widths if \code{plot_method = "maxQ"} or \code{plot_method = "allQ"} is used, sets relative widths of kriged TESS map and legend (defaults to 3:1), from \link[cowplot]{plot_grid}
#'
#' @family TESS functions
#'
#' @return ggplot object of TESS results
#' @export
tess_ggplot <- function(krig_admix, coords = NULL, plot_method = "maxQ", ggplot_fill = algatr_col_default("ggplot"), minQ = 0.10, plot_axes = FALSE, rel_widths = c(3, 1)) {
  # Set up ggplot df
  gg_df <- krig_admix %>%
    terra::as.data.frame(x, xy = TRUE, na.rm = FALSE) %>%
    tidyr::as_tibble() %>%
    tidyr::pivot_longer(names_to = "K", values_to = "Q", -c(x, y)) %>%
    dplyr::mutate(K = as.factor(gsub("K", "", K))) %>%
    dplyr::group_by(x, y)

  # Use max or all Q
  if (plot_method == "maxQ" | plot_method == "maxQ_poly") gg_df <- gg_df %>% dplyr::top_n(1, Q)
  if (plot_method == "allQ" | plot_method == "allQ_poly") gg_df <- gg_df %>% dplyr::filter(Q > minQ)

  # Set up base plot
  plt <- ggplot2::ggplot()

  # Plot as polygon or continuous Q
  if (plot_method == "maxQ_poly" | plot_method == "allQ_poly") {
    plt <- plt + ggplot2::geom_tile(data = gg_df, ggplot2::aes(x = x, y = y, fill = K), alpha = 0.5)
  } else {
    plt <- plt +
      ggplot2::geom_tile(data = gg_df, ggplot2::aes(x = x, y = y, fill = K, alpha = Q)) +
      ggplot2::scale_alpha_binned(breaks = round(seq(0, 1, by = 0.10), 1))
  }

  # Add color
  plt <- plt + ggplot_fill

  # Add themes and coord controls
  plt <- plt + ggplot2::coord_equal() + ggplot2::theme_bw()

  # Add axes
  if (plot_axes) {
    plt <- plt + ggplot2::theme(
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      aspect.ratio = 1
    )
  }

  if (!plot_axes) {
    plt <- plt + ggplot2::theme(
      panel.grid.minor.y = ggplot2::element_blank(),
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
      aspect.ratio = 1
    )
  }

  # Add coords
  if (!is.null(coords)) plt <- plt + ggplot2::geom_point(data = data.frame(coords), ggplot2::aes(x = x, y = y))

  # Produce plot with krig_legend for "allQ" or "maxQ"
  if (plot_method == "allQ" | plot_method == "maxQ") {
    # Remove existing legend
    plt <- plt +
      ggplot2::theme(legend.position = "none")

    # Add secondary plot (which will become the legend) with combined K and Q values using helper function
    plt_leg <- krig_legend(gg_df = gg_df, plot_method = plot_method, ggplot_fill = ggplot_fill, minQ = minQ)

    plt <- cowplot::plot_grid(plt, plt_leg, rel_widths = rel_widths)
  }

  return(plt)
}

#' Helper function to make a custom legend for TESS maps
#'
#' @param gg_df dataframe in tidy format of Q values from \link[algatr]{tess_ggplot}
#' @inheritParams tess_ggplot
#'
#' @family TESS functions
#'
#' @return legend for kriged map from TESS
#' @export
#' @family TESS functions
#' @importFrom ggplot2 '%+replace%'
krig_legend <- function(gg_df, plot_method, ggplot_fill, minQ){
  if (plot_method == "maxQ") vals <- seq(0, 1, by = 0.10)
  if (plot_method == "allQ") vals <- seq(minQ, 1, by = 0.10)
  kvals <- 1:length(unique(gg_df$K))
  dat <- as.data.frame(tidyr::expand_grid(vals, kvals))
  dat$kvals <- as.character(dat$kvals)
  dat$vals <- as.character(dat$vals)

  plt_leg <-
    dat %>%
    ggplot2::ggplot(ggplot2::aes(x = kvals, y = vals, fill = kvals, alpha = vals, group = kvals)) +
    ggplot2::geom_raster() +
    ggplot2::scale_y_discrete(expand = c(0, 0), name = "Q", breaks = vals) +
    ggplot2::scale_x_discrete(expand = c(0, 0), name = "K", breaks = kvals, labels = kvals) +
    ggplot2::coord_fixed(ratio = 1) +
    cowplot::theme_cowplot() %+replace% ggplot2::theme(legend.position = "none",
                                                       axis.ticks = ggplot2::element_blank(),
                                                       axis.line = ggplot2::element_blank()) +
    ggplot_fill

  return(plt_leg)
}

#' Plot all kriged Q values for each K
#'
#' @param krig_admix RasterStack returned by \link[algatr]{tess_krig}
#' @param coords dataframe with x and y coordinates for plotting (optional)
#' @param ... Graphical parameters. Any argument that can be passed to image.plot and to base plot
#' @inheritParams tess_do_everything
#'
#' @export
#' @family TESS functions
tess_plot_allK <- function(krig_admix, coords = NULL, col_pal = algatr_col_default("base"), col_breaks = 20, ...) {
  # Get K
  K <- terra::nlyr(krig_admix)

  # Plot kriged admixture maps while masking small values (e.g. < minQ)
  purrr::walk(1:K, allK_plot_helper, krig_admix, coords = coords, col = col_pal(K), col_breaks = col_breaks, ...)
}

#' Helper function for all K plotting
#'
#' @inheritParams tess_plot_allK
#' @param K K value
#' @param col single color code
#'
#' @export
#' @family TESS functions
allK_plot_helper <- function(K, krig_admix, coords = NULL, col, col_breaks, ...) {
  # Suppress irrelevant plot warnings
  suppressWarnings({
    terra::plot(krig_admix[[K]],
      col = make_plot_col(K, col, col_breaks, alpha = 1, start_col = rgb(0.94, 0.94, 0.95, 1)),
      range = c(0, max(terra::minmax(krig_admix)["max", ])),
      axes = FALSE,
      box = FALSE,
      main = paste0("K = ", K),
      ...
    )
  })

  # Add coordinates, if provided
  if (!is.null(coords)) points(coords, pch = 3)
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
make_plot_col <- function(K, col, col_breaks, poly = FALSE, alpha = 0, start_col = rgb(1, 1, 1, alpha)) {
  if (poly) {
    # Make color palette using only solid color
    plot_col <- col[K]
  } else {
    # Make color palette gradient from transparent to color defined by K
    kpal <- colorRampPalette(c(start_col, col[K]), interpolate = "linear", alpha = TRUE)
    plot_col <- kpal(col_breaks)
  }

  return(plot_col)
}

#' Create TESS barplot
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
#' @return STRUCTURE-style bar plot of TESS results
#' @export
#'
#' @family TESS functions
tess_barplot <- function(qmat, col_pal = algatr_col_default("base"), sort_by_Q = TRUE, legend = TRUE, legend_position = "bottomright", border = NA, space = 0, ...) {
  # CODE ADAPTED FROM: https://github.com/bcm-uga/TESS3_encho_sen/blob/master/R/plotQ.R

  # Get K
  K <- ncol(qmat)

  if (sort_by_Q) {
    gr <- apply(qmat, MARGIN = 1, which.max)
    gm <- max(gr)
    gr.o <- order(sapply(1:gm, FUN = function(g) mean(qmat[, g])))
    gr <- sapply(gr, FUN = function(i) gr.o[i])
    or <- order(gr)
    Qm <- t(qmat[or, ])
    class(Qm) <- "matrix"
    graphics::barplot(Qm, col = col_pal(K), border = border, space = space, ...)
    legend("bottomright", pch = 15, legend = paste0("K = ", 1:K), col = col_pal(K))
    return(list(order = or))
  } else {
    Qm <- t(qmat)
    class(Qm) <- "matrix"
    graphics::barplot(Qm, col = col_pal(ncol(qmat)), border = border, space = space, ...)
    legend(legend_position, pch = 15, legend = paste0("K = ", 1:K), col = col_pal(K))
    return(list(order = 1:nrow(qmat)))
  }
}

#' Create TESS barplot using ggplot2
#'
#' @param qmat Q matrix
#' @param ggplot_fill any ggplot2 scale fill discrete function (default: \link[algatr]{scale_fill_viridis_d}, \code{option = "turbo"})
#' @param sort_by_Q whether to sort bars by Q value (equivalent to \link[tess3r]{barplot} sort.by.Q)
#' @param legend whether to display legend (defaults to TRUE)
#'
#' @return ggplot object of TESS results as a barplot
#'
#' @family TESS functions
#' @export
tess_ggbarplot <- function(qmat, ggplot_fill = algatr_col_default("ggplot"), sort_by_Q = TRUE, legend = TRUE) {
  # Get K
  K <- ncol(qmat)

  dat <- as.data.frame(qmat)
  dat <- dat %>%
    tibble::rownames_to_column(var = "order")

  if (sort_by_Q) {
    gr <- apply(qmat, MARGIN = 1, which.max)
    gm <- max(gr)
    gr.o <- order(sapply(1:gm, FUN = function(g) mean(qmat[, g])))
    gr <- sapply(gr, FUN = function(i) gr.o[i])
    or <- order(gr)

    dat <- dat %>%
      dplyr::arrange(factor(order, levels = or))
    dat$order <- factor(dat$order, levels = dat$order)
  }

  # Make into tidy df
  gg_df <-
    dat %>%
    tidyr::pivot_longer(names_to = "K_value", values_to = "Q_value",
                        -c(order))

  # Build plot using helper function
  plt <- ggbarplot_helper(gg_df) + ggplot_fill

  # Remove legend
  if (!legend) plt <- plt + ggplot2::theme(legend.position = "none")

  return(plt)
}

#' Helper function for TESS barplots using ggplot
#'
#' @param dat Q matrix
#'
#' @return barplot with Q-values and individuals, colorized by K-value
#'
#' @family TESS functions
#' @export
ggbarplot_helper <- function(dat) {
  dat %>%
    ggplot2::ggplot(ggplot2::aes(x = order, y = Q_value, fill = K_value)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_y_continuous(expand = c(0,0)) +
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA, colour = "black", linetype = "solid", linewidth = 1.5),
                   strip.text.y = ggplot2::element_text(size = 30, face = "bold"),
                   strip.background = ggplot2::element_rect(colour = "white", fill = "white"),
                   panel.spacing = ggplot2::unit(-0.1, "lines"))
}

#' Best K Selection based on cross entropy
#'
#' @param tess3_obj list produced by \code{\link{tess3}}
#' @param Kvals vector of K values for testing
#'
#' @note (source: https://chazhyseni.github.io/NALgen/post/determining_bestk/)
#' @export
#' @family TESS functions
bestK <- function(tess3_obj, Kvals) {
  ce <- list()
  for (k in Kvals) ce[[k]] <- tess3_obj[[k]]$crossentropy
  ce.K <- c()
  for (k in Kvals) ce.K[k] <- min(ce[[k]])
  diff <- ce.K[-1] - ce.K[-max(Kvals)]
  slope <- exp(-diff) - 1
  # K is selected based on the smallest slope value in the upper quartile
  K <- min(which(slope <= quantile(slope)[4]))
  return(K)
}

#' Helper function to get population assignments for best K
#'
#' @param gen genotype dosage matrix (rows = individuals & columns = SNPs) or `vcfR` object
#' @param tess3_obj list produced by \code{\link{tess3}}
#' @param K K value
#'
#' @return population assignments for each individual based on max Q values
#'
#' @export
#' @family TESS functions
pops_helper <- function(gen, tess3_obj, K) {
  # Get individual names for population assignments
  if (inherits(gen, "vcfR")) names <- colnames(gen@gt[,-1])
  if (inherits(gen, "matrix")) names <- rownames(gen)

  # Get qmatrix and make into df
  qmat <- tess3r::qmatrix(tess3_obj, K = K)
  qmat <- as.data.frame(qmat)
  # Replace Vs with Ks for clarity
  colnames(qmat) <- stringr::str_replace_all(colnames(qmat), "V", "K")

  pops <- dplyr::bind_cols(names, qmat) %>%
    dplyr::rename(individual = `...1`)

  # Get population assignment based on max Q value
  pops %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pop_assignment = which.max(dplyr::c_across(-individual)))

  return(pops)
}

#' Create default TESS color palette
#' @param n number of colors to generate (must be less than 9)
#'
#' @export
#' @family TESS functions
tess_col_default <- function(n, alpha = 1) {
  if (n > 9) stop("The default color palette expects less than 9 values")
  tessCP <- CreatePalette()
  tesscol <- sapply(1:n, function(x, tessCP) {
    tessCP[[x]][9]
  }, tessCP)
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
algatr_col_default <- function(x) {
  if (x == "ggplot") col <- ggplot2::scale_fill_viridis_d(option = "turbo", begin = 0.1, end = 0.9)
  if (x == "base") col <- function(n, alpha = 1, begin = 0, end = 1, direction = 1) viridis::viridis(n, alpha, begin = 0.1, end = 0.9, direction, option = "turbo")
  return(col)
}
