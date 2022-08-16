
# TODO:
# - FIX PROJECTION WEIRDNESS -> have inputs be projected coords with same projection as SPDF or conversion?
# - FIX ARGUMENTS/WEIRD DF
# - CHANGE FUNCTION NAME TO tess_do_everything
# - change k_selection to K_selection

#' TESS function to do everything
#'
#' @param gen genotype matrix
#' @param coords dataframe with x and y coordinates
#' @param Kvals vector of K values to test
#' @param spdf Spatial Points Dataframe to create grid for kriging
#' @param ploidy ploidy of organism (defaults to 2)
#' @param n_cell number of grid cells to use when kriging
#' @param k_selection how to perform k selection (options: "auto" for automatic selection based on \code{\link{bestK}} or "manual" to enter into console)
#' @param plot_method method for making rainbow map of kriged layers (options: "max" to only plot the max Q value for each cell or "all" to plot all Qvalues greater than \code{minQ})
#' @param minQ threshold for minimum Q-value for rainbow plotting if \code{method = "all"} is used
#'
#' @family TESS functions
#'
#' @return list with all TESS results, final K value, and final kriged raster
#' @export
#'
#' @examples
tess_full <- function(gen, coords, spdf, Kvals = 1:10,
                      ploidy = 2, n_cell = 1000, k_selection = "auto",
                      plot_method = "max", minQ = 0,
                      tess_method = "projected.ls"){

  # RUN TESS ---------------------------------------------------------------------------------------------------
  # test different k values if more than one provided
  if(length(Kvals)>1){
    tess_results <- tess_ktest(gen, coords, Kvals = Kvals, tess_method = tess_method, ploidy = ploidy, k_selection = k_selection)

    # get K
    K <- tess_results[["K"]]

    # get tessobj
    tess3_obj <- tess_results[["tess3_obj"]]

    }

  # if only one value provided just use that
  # TODO: FIX THIS TO NOT MAKE WEIRD DF (MAYBE SWITCH BACK FUNCTIONS TO TAKE MORE ARGUMENTS)
  if(length(Kvals) == 1){

    # K is just Kvals if there is only one value
    K = Kvals

    # run tess for given K value
    tess3_obj <- tess3(X = gen, coord = coords, K = Kvals, method = tess_method, ploidy = ploidy)

  }

  # KRIGE QMATRIX  -----------------------------------------------------------------------------------------------

  # get Qmatrix
  qmat <- qmatrix(tess3_obj, K = K)

  # krige Qmatrix
  krig_admix <- tess_krig(qmat = qmat, coords = coords, spdf = spdf, n_cell = n_cell)


  # PLOTS --------------------------------------------------------------------------------------------------------

  # set col pal
  col_pal = turbo(K)

  # plot all layers
  # TODO: FIX THIS SO COORDS ARE PROJECTED Consistently
  coords_spdf <- coord_proj(coords, spdf, crop_to_spdf = TRUE)

  # plot all kriged maps for each K
  par(pty = "s", mar = rep(1,4), mfrow = c(1,2))
  tess_plot_allK(krig_admix, coords_spdf = coords_spdf, spdf = spdf, col_pal = col_pal)

  # make rainbow plots
  par(pty = "s", mar = rep(0,4), mfrow = c(1,2))
  tess_rainbow_plot(krig_admix, coords = coords_spdf, spdf = spdf, plot_method = "max", minQ = minQ, col_pal = col_pal)
  tess_rainbow_plot(krig_admix, coords = coords_spdf, spdf = spdf, plot_method = "all", minQ = minQ, col_pal = col_pal)

  # make barplot
  par(pty = "m", mfrow = c(1,1))
  tess_barplot(qmat = qmat, col_pal = col_pal)

  # make maxQ plot
  par(pty = "m", mfrow = c(1,1), mar = rep(0,4))
  maxQ <- map_maxQ(krig_admix, plot_map = TRUE)

  # OUTPUTS ------------------------------------------------------------------------------------------------------

  # create list with all outputs
  tess_results <- list(K = K,
                       Qmatrix = qmat,
                       krig_admix = krig_admix,
                       tess_results = tess3_obj,
                       coords = coords,
                       Kvals = Kvals,
                       spdf = spdf,
                       maxQ = maxQ)

  return(tess_results)
}


#' Test multiple K values
#'
#' @param gen genotype matrix
#' @param coords dataframe with x and y coordinates
#' @param Kvals vector of K values to test
#' @param tess_method method provided to tess
#' @param ploidy ploidy of organism (defaults to 2)
#' @param k_selection how to perform k selection (options: "auto" for automatic selection based on \code{\link{bestK}} or "manual" to enter into console)
#'
#' @return
#' @export
#'
#' @examples
tess_ktest <- function(gen, coords, Kvals = 1:10, spdf = NULL, tess_method = "projected.ls", ploidy = 2, k_selection = "auto"){

  # run tess for all K values
  tess3_obj <- tess3(X = gen, coord = coords, K = Kvals, method = tess_method, ploidy = ploidy)

  # plot CV results and mark the K-value automatically selected
  plot(tess3_obj, pch = 19, col = "blue",
       xlab = "Number of ancestral populations",
       ylab = "Cross-validation score")

  # get best K value
  if(k_selection == "auto"){K <-  bestK(tess3_obj, Kvals)}
  if(k_selection == "manual"){K <- as.numeric(readline(prompt = "Enter K Value: "))}

  abline(v = K, col = "red", lty = "dashed")

  # create list with tess3 object and K value
  tess_results <- list(K = K,
                    tess3_obj = tess3_obj,
                    coords = coords,
                    Kvals = Kvals,
                    spdf = spdf)

  return(tess_results)
}


#' Plot results from K value testing
#'
#' @param tess_results list produced by \code{\link{tess_ktest}} with Kvals, coords, and spdf
#' @param n_cell number of grid cells to use when kriging
#' @param plot_method method for making rainbow map of kriged layers (options: "max" to only plot the max Q value for each cell or "all" to plot all Qvalues greater than \code{minQ})
#'
#' @return
#' @export
#'
#' @examples
tess_ktest_plot <- function(tess_results, n_cell = 1000, plot_method = "max"){

  # get objects from tess results list
  Kvals <- tess_results[["Kvals"]]
  coords <- tess_results[["coords"]]
  spdf <- tess_results[["spdf"]]

  # krige and plot results
  par(pty = "s", mfrow = c(3,3), mar = rep(0,4), oma = rep(0,4))
  for(k in Kvals[-1]){
    # get Qmatrix
    qmat <-  qmatrix(tess_results[["tess3_obj"]], k)

    # krige qmatrix
    krig_admix <- tess_krig(qmat = qmat, coords = coords, spdf = spdf, n_cell = n_cell)

    # plot all layers
    raster::plot(krig_admix, col = plasma(100), zlim = c(0, max(maxValue(krig_admix))), axes = FALSE, box = FALSE)

    # make rainbow plots
    tess_rainbow_plot(krig_admix, spdf = spdf, plot_method = plot_method, minQ = minQ)

    # make barplot
    pal <- CreatePalette(turbo(K))
    barplot(qmat, sort.by.Q = TRUE,
            border = NA, space = 0,
            col.palette = pal,
            xlab = "Individuals", ylab = "Ancestry coefficients")
  }
}


#' Krige admixture values
#'
#' @param qmat qmatrix
#' @param coords dataframe with x and y coordinates
#' @param spdf Spatial Points Dataframe to create grid for kriging
#' @param n_cell number of grid cells to use when kriging
#'
#' @return
#' @export
#'
#' @examples
tess_krig <- function(qmat, coords, spdf, n_cell = 10000){

  # define K
  K <- ncol(qmat)

  # make grid for kriging
  krig_grid <- spdf_to_grid(spdf, n_cell = n_cell)

  krig_admix_r <- stack()
  krig_admix_df <- data.frame()
  for(k in 1:K){
    # make qmat into spdf
    # TODO: FIX PROJECTION STUFF
    krig_df <- coord_proj(coords, spdf)
    krig_df$prop <- qmat[,k]

    # Skip if props are identical (kriging not possible)
    if(unique(krig_df$prop) == 1){
      # TODO: COME BACK AND FIX THIS SO THAT A BLANK RASTER IS ADDED
      # COME BACK AND FIX THIS SO THAT A BLANK RASTER IS ADDED
      warning(paste0("Only one unique Q value for K = ", k, ", skipping (note: may want to consider different K value)"))
      next
    }

    # Krige
    krig_res <- autoKrige(prop ~ 1, krig_df, krig_grid)
    krig_spdf <- krig_res$krige_output

    # turn SPDF into raster and stack
    krig_raster <- raster::rasterFromXYZ(krig_spdf, crs = raster::crs(krig_grid))
    krig_admix_r <- raster::stack(krig_admix_r, krig_raster)

    # save SPDF
    krig_df <- data.frame(krig_spdf)
    # add K value
    krig_df$K <- k
    # add to df
    krig_admix_df <- bind_rows(krig_admix_df, krig_df)

  }

  # TODO: DECIDE IF YOU WANT THIS: convert all values greater than 1 to 1 and all values less than 0 to 0
  # DECIDE IF YOU WANT THIS: convert all values greater than 1 to 1 and all values less than 0 to 0

  # rename layers
  names(krig_admix_r) <- paste0("K",1:K)

  # convert K to factor
  krig_admix_df$K <- factor(krig_admix_df$K)
  # add column names
  colnames(krig_admix_df) <- c("x", "y", "Q", "var", "stdev", "optional", "K")

  # return stack and dataframe
  krig_admix <- list(raster = krig_admix_r, dataframe = krig_admix_df)

  return(krig_admix)
}

#' Make rainbow TESS plot from kriged admixture plots
#'
#' @param krig_admix list returned by  \code{\link{tess_krig}}
#' @param coords dataframe with x and y coordinates for plotting (not required)
#' @param spdf Spatial Points Dataframe for plotting (not required)
#' @param plot_method method for making rainbow map of kriged layers (options: "max" to only plot the max Q value for each cell or "all" to plot all Qvalues greater than \code{minQ})
#' @param minQ threshold for minimum Q-value for rainbow plotting if \code{method = "all"} is used
#' @param alpha transparency for plotting
#'
#' @return
#' @export
#'
#' @examples
tess_rainbow_plot <- function(krig_admix, coords = NULL, spdf = NULL, plot_method = "max", minQ = 0.10, alpha = 1, col_pal = "Default", reclassify = FALSE, poly = FALSE){

  # get raster stack
  r <- krig_admix[["raster"]]

  # determine K based on_cell number of layers
  K <- nlayers(r)

  # set col pal
  if(col_pal == "Default"){
    col_pal <- turbo(K, alpha)
  }

  # get extent of raster (to set up plot)
  ext <- extent(r)
  # set up base plot (! important ! Don't remove or things will get wonky)
  plot(1, legend = FALSE, axes = FALSE, box = FALSE, type = "n",
       xlim = c(ext[1], ext[2]), ylim = c(ext[3], ext[4]),
       xlab = "", ylab = "")

  # Select method and options
  if(plot_method == "max"){tess_plot_max(krig_admix, K = K, reclassify = reclassify, poly = FALSE, col_pal = col_pal)}
  if(plot_method == "all"){tess_plot_all(krig_admix, K = K, reclassify = reclassify, poly = FALSE, minQ = minQ, col_pal = col_pal)}
  if(plot_method == "max_poly"){tess_plot_max(krig_admix, K = K, reclassify = reclassify, poly = TRUE, col_pal = col_pal)}
  if(plot_method == "all_poly"){tess_plot_all(krig_admix, K = K, reclassify = reclassify, poly = TRUE, minQ = minQ, col_pal = col_pal)}


  # add coordinates if given
  if(!is.null(coords)){points(coords, pch = 16)}
  # add spdf boundary if
  if(!is.null(spdf)){raster::plot(spdf, add = TRUE)}
}


#' Plot method: Q max
#'
#' @describeIn tess_rainbow_plot plot based on max of Qvalues
#'
#' @param krig_admix
#' @param K
#' @param reclassify
#' @param poly
#' @param col_pal
#'
#' @return
#' @export
#'
#' @examples
tess_plot_max <- function(krig_admix, K = K, reclassify = FALSE, poly = FALSE, col_pal){

  # max and min raster values for plotting
  # TODO: FIX THIS
  maxr <- max(maxValue(krig_admix[["raster"]]))
  minr <- min(minValue(krig_admix[["raster"]]))

  # summarize dataframe by only retaining highest Q values for each point
  pop_df <-  krig_admix[["dataframe"]] %>%
    group_by(x, y) %>%
    top_n(1, Q)

  # plot each kriged admixture map one by one on top of each other
  # TODO: convert to purrr::map
  for(i in 1:K){
    pop_spdf <- pop_df[pop_df$K == i, ]

    # skip to next iteration if there is no more than one value for that K
    # TODO: add warning about skipping
    if(nrow(pop_spdf) < 2) next

    # make into spdf and convert to raster
    coordinates(pop_spdf) <- ~x+y
    rl <- rasterFromXYZ(pop_spdf)

    if(reclassify){
      rl <- tess_reclass(rl)

      # reset max and min raster values for plotting
      maxr <- max(maxValue(rl))
      minr <- min(minValue(rl))
    }

    # make color palette (if poly use solid color, if not use gradient)
    if(poly){
      cols <- col_pal[i]

    } else {
      colors <- c(rgb(1,1,1,alpha = 0), col_pal[i])
      kpal <- colorRampPalette(colors, interpolate="linear", alpha = TRUE)
      cols <- kpal(100)
      }

    # plot raster
    raster::plot(rl,
         add = TRUE,
         legend = FALSE,
         col = cols,
         main = plot_method,
         zlim = c(minr, maxr))

  }

}

#' Plot method: Q all
#'
#' @describeIn tess_rainbow_plot plot based on all Qvalues
#'
#' @param krig_admix
#' @param K
#' @param reclassify
#' @param poly
#' @param minQ
#' @param col_pal
#'
#' @return
#' @export
#'
#' @examples
tess_plot_all <- function(krig_admix, K = K, reclassify = FALSE, poly = FALSE, minQ = 0.1, col_pal){
tess_plot_all <- function(krig_admix, K = K, reclassify = FALSE, poly = FALSE, minQ = 0.10, col_pal){

  rl <- krig_admix[["raster"]]

  if(reclassify) rl <- tess_reclass(rl)

  # max of raster values for plotting
  maxr <- max(maxValue(rl))

  if(reclassify){
    rl <- tess_reclass(rl)

    # max and raster values for plotting
    maxr <- max(maxValue(rl))
  }


  # plot kriged admixture maps while masking small values (e.g. < minQ)
  for(i in 1:K){

    # make color palette (if poly use solid color, if not use gradient)
    if(poly){
      cols <- col_pal[i]
    } else {
      colors <- c(rgb(1,1,1,alpha = 0), col_pal[i])
      kpal <- colorRampPalette(colors, interpolate="linear", alpha = TRUE)
      cols <- kpal(100)
    }

    raster::plot(rl[[i]],
         col = cols,
         main = plot_method,
         add = TRUE,
         legend = FALSE,
         zlim = c(minQ, maxr))
         zlim = c(minr, maxr))

  }
}

tess_plot_allK <- function(krig_admix, K = K, coords_spdf = NULL, spdf = NULL, col_pal = turbo(K)){

  r <- krig_admix[["raster"]]
  K <- nlayers(r)

  for(i in 1:K){
    cols <- c(rgb(1,1,1,alpha = 0), col_pal[i])

    kpal <- colorRampPalette(cols, interpolate="linear", alpha = TRUE)

    raster::plot(krig_admix[["raster"]][[i]],
         col = kpal(100),
         zlim = c(0, max(maxValue(r))),
         axes = FALSE,
         box = FALSE,
         main = paste0("K=",i))

    # add coordinates if given
    if(!is.null(coords_spdf)){points(coords_spdf, pch = 16)}
    # add spdf boundary if
    if(!is.null(spdf)){raster::plot(spdf, add = TRUE)}
  }
}


#' Creat TESS barplot
#'
#' @inheritParams tess_doEverything
#'
#' @return
#' @export
#'
#' @examples
tess_barplot <- function(qmat, col_pal = turbo(ncol(qmat))){
  pal <- CreatePalette(col_pal)
  barplot(qmat,
          sort.by.Q = TRUE,
          border = NA,
          space = 0,
          col.palette = pal,
          xlab = "Individuals", ylab = "Ancestry coefficients")
}


# TODO: CHECK THIS
#' Best K Selection based on cross entropy
#'
#' @param tess3_obj list produced by \code{\link{tess3}}
#' @param Kvals vector of K values for testing
#'
#' @note  (source: https://chazhyseni.github.io/NALgen/post/determining_bestk/)
#' @return
#' @export
#'
#' @examples
bestK <- function(tess3_obj, Kvals){
  ce <- list()
  for(k in Kvals) ce[[k]] <- tess3_obj[[k]]$crossentropy
  ce.K <- c()
  for(k in Kvals) ce.K[k] <- min(ce[[k]])
  diff <- ce.K[-1] - ce.K[-max(Kvals)]
  slope <- exp(-diff) - 1
  #K is selected based on the smallest slope value in the upper quartile
  K <- min(which(slope <= quantile(slope)[4]))
  return(K)
}

#' Reclassify Kriged Qmatrix Raster
#'
#' @param r kriged Q raster
#' @param inc increment for reclassification bins
#' @param rec reclassify matrix (if not provided, will be set to default)
#'
#' @return
#' @export
#'
#' @examples
tess_reclass <- function(r, inc = 0.05, rec = "Default"){

  if(rec == "Default"){
    v1 <- seq(0,0.9,inc)
    v2 <- seq(0.1,1,inc)
    rec <- cbind(v1, v2, v1)
  }

  r <- reclassify(r, rec)
}



# FUNCTIONS IN TESTING ---------------------------------------------------------
stack_maxQ <- function(krig_list){
  K <- length(krig_list)
  rs <- stack()
  for(i in 2:K){
    krig_admix <- krig_list[[i]]
    r <- map_maxQ(krig_admix)
    rs <- stack(rs, r)
  }
  return(rs)
}

#' Make raster with max Q values
#'
#' @param krig_admix
#' @param alpha
#' @param zlim
#' @param add
#' @param cols
#' @param plot_map
#'
#' @return
#' @export
#'
#' @examples
map_maxQ <- function(krig_admix, alpha = alpha, add = FALSE, cols = "Default", plot_map = FALSE){
  # summarize dataframe by only retaining highest Q values for each point
  pop_spdf <-  krig_admix[["dataframe"]] %>%
    group_by(x, y) %>%
    top_n(1, Q)

  # make into spdf and convert to raster
  coordinates(pop_spdf) <- ~x+y
  rl <- rasterFromXYZ(pop_spdf)

  if(plot_map ){
    if(cols == "Default"){
      cols <- magma(100)
    }

    raster::plot(rl, col = cols, add = add, legend = FALSE, axes = FALSE, box = FALSE)
  }


  return(rl)
}

#' Make raster of population boundaries
#'
#' @param krig_admix
#' @param K
#' @param m
#' @param minQ
#' @param boundaries
#'
#' @return
#' @export
#'
#' @examples
make_boundary_raster <- function(krig_admix, K, m = matrix(rep(1,9), nrow = 3), minQ = NULL, boundaries = FALSE){
  r <-  krig_admix[["raster"]][[1]]
  r[] <- NA

  pop_df <-  krig_admix[["dataframe"]] %>%
    group_by(x, y) %>%
    top_n(1, Q)

  rs <- r
  # plot each kriged admixture map one by one on top of each other
  for(i in 1:K){
    pop_spdf <- pop_df[pop_df$K == i, ]

    # skip to next iteration if there are no values for that K
    if(nrow(pop_spdf) == 0) next
    if(length(unique(pop_spdf$x)) == 1) next
    if(length(unique(pop_spdf$y)) == 1) next

    # make into spdf and convert to raster
    coordinates(pop_spdf) <- ~x+y
    rl <- rasterFromXYZ(pop_spdf)

    # mask areas where Q is high (pop assignment certain)
    if(is.numeric(minQ)){rl[rl > minQ] <- minQ}

    # only plot boundaries beween pops
    if(boundaries){rl <- rl*0+i}

    # stitch rasters together
    rs <- mosaic(rs, rl, fun = min)

  }

  rv <- focal(rs, w = m, fun = var)

  # only plot boundaries beween pops
  if(boundaries){rv[rv != 0] <- 1}

  return(rv)

}

#' Make boundary rasters for a set of K values
#'
#' @param krig_list
#' @param minQ
#' @param boundaries
#'
#' @return
#' @export
#'
#' @examples
make_boundaries_allK <- function(krig_list, minQ = NULL, boundaries = FALSE){
  K = length(krig_list)

  rv_stack <- stack()
  for(k in 2:K){
    krig_admix <- krig_list[[k]]
    rv <- make_boundary_raster(krig_admix, K = k, minQ = minQ, boundaries = boundaries)
    rv_stack <- stack(rv_stack, rv)
  }
  names(rv_stack) <- paste0("K", 2:K)

  return(rv_stack)
}

#' Plot boundaries
#'
#' @param rv_stack
#' @param option
#' @param log_transform
#' @param col
#' @param weights
#' @param tess3_obj
#' @param main
#'
#' @return
#' @export
#'
#' @examples
plot_boundaries <- function(rv_stack, option = "all", log_transform = FALSE,
                            col = magma(100),
                            weights = NULL, tess3_obj = NULL,
                            main = ""){
  if(option == "all"){
    rv_summary <- rv_stack
  }

  if(option == "mean"){
    rv_summary <-  raster::mean(rv_stack)
  }

  if(option == "weighted.mean"){
    if(!is.null(tess3_obj)){
      ce_vec <- unlist(lapply(tess3_obj, function(x){x <- x$crossentropy}))
      weights <- 1/ce_vec[-1]
    }

    rv_summary <-  weighted.mean(rv_stack, w = weights[1:nlayers(rv_stack)])
  }

  if(log_transform){
    rv_summary <- log(rv_summary)
  }

  raster::plot(rv_summary, col = magma(100), axes = FALSE, box = FALSE, main = main)
}

