#' MMRR function to do everything ()
#' TODO: gendist matrix must range between 0 and 1 as with GDM? I don't think so - APB
#'
#' @param gendist matrix of genetic distances
#' @param coords dataframe with x and y coordinates
#' @param env dataframe with environmental values for each coordinate; if not provided it will be calculated based on coords/envlayers
#' @param envlayers rasters for for extracting environmental values using coordinates if `env` isn't provided
#' @param model whether to fit the model with all variables (`"full"`) or to perform variable selection to determine the best set of variables (`"best"`); default = "best"
#' @param nperm number of permutations to use to calculate variable importance; only used if `model = "best"` (default = 999)
#' @param stdz if TRUE then matrices will be standardized. Default = TRUE.
#'
#' @return
#' @export
#'
#' @examples
mmrr_do_everything <- function(gendist, coords, env = NULL, envlayers = NULL, model = "best", nperm = 999, stdz = TRUE){

  # If not provided, make env data frame from layers and coords
  if(is.null(env)){env <- raster::extract(envlayers, coords)}

  # Convert to tibble for map function
  env <- dplyr::as_tibble(env)

  # Make env dist matrix
  Xdist <- purrr::map(env, env_dist)
  Xdist[["geodist"]] <- geo_dist(coords)

  # Make geodist mat
  Ydist <- as.matrix(gendist)

  # Run MMRR with variable selection
  if(model == "best"){
    mod <- mmrr_var_sel(Ydist, Xdist, nperm = nperm, stdz = stdz)

    # If NULL, exit with NULL
    if(is.null(mod)) return(NULL)

    # subset Xdist with significant variables
    # note: the list() and names() part are necessary to ensure a named list is produced if there is only one significant variable
    Xdist_best <- list(Xdist[[names(mod$coefficients)[-1]]])
    names(Xdist_best) <- names(mod$coefficients)[-1]

    # plot results
    plot_mmrr(Y = Ydist, X = Xdist_best, mod = mod, stdz = stdz)

  }

  if(model == "full"){
    mod <- MMRR(Ydist, Xdist, nperm = nperm)

    # If NULL, exit with NULL
    if(is.null(mod)) return(NULL)

    # plot results
    plot_mmrr(Y = Ydist, X = Xdist, mod = mod, stdz = stdz)
    }

  # Make nice data frame
  coeff_df <- mmrr_df(mod)

  # Make results list
  results <- list(coeff_df = coeff_df,
                  mod = mod,
                  Ydist = Ydist,
                  Xdist = Xdist)

  if(model == "best") results[["Xdist_best"]] <- Xdist_best

  return(results)

}

#' mmrr_var_sel performs MMRR with backward elimination variable selection
#' @param Y is a dependent distance matrix
#' @param X is a list of independent distance matrices (with optional names)
#' @param nperm number of permutations to use to calculate variable importance; only used if `model = "best"` (default = 999)
#' @param stdz if TRUE then matrices will be standardized. Default = TRUE.
mmrr_var_sel <- function(Y, X, nperm = 999, stdz = TRUE){
  # Fit full model
  mmrr.model <- MMRR(Y, X, nperm = nperm, scale = stdz)
  pvals <- mmrr.model$tpvalue[-1] # Remove intercept p-value

  # Eliminate variable with highest p-value, re-fit, and continue until only significant variables remain or no variables remain
  while((max(pvals) > 0.05) & (length(pvals) > 1)){
    print(pvals)
    rem.var <- which(pvals == max(pvals))
    X <- X[-rem.var]
    if(length(X) == 0) break
    mmrr.model <- MMRR(Y, X, nperm = nperm, scale = stdz)
    pvals <- mmrr.model$tpvalue[-1]
  }

  # Repetitive test, but just to be sure
  if(length(X) == 0 | (length(pvals) == 1) & all(pvals > 0.05)){warning("No significant variable combo found, returning NULL object"); mmrr.model <- NULL}

  return(mmrr.model)
}

#' MMRR performs Multiple Matrix Regression with Randomization analysis
#' @param Y is a dependent distance matrix
#' @param X is a list of independent distance matrices (with optional names)
#' @param nperm is the number of permutations to be used in significance tests. Default = 999.
#' @param scale if TRUE then matrices will be standardized. Default = TRUE.
#' @details
#' When using MMRR, please cite the original citation:
#' Wang I.J. (2013) Examining the full effects of landscape heterogeneity on spatial genetic variation: a multiple matrix regression approach for quantifying geographic and ecological isolation. Evolution, 67: 3403-3411.
#' @export
MMRR <- function(Y, X, nperm = 999, scale = TRUE){
  # Compute regression coefficients and test statistics
  nrowsY <- nrow(Y)
  y <- unfold(Y, scale)
  if(is.null(names(X)))names(X) <- paste("X", 1:length(X), sep="")
  Xmats <- sapply(X, unfold, scale = scale)
  fit <- stats::lm(y ~ Xmats)
  coeffs <- fit$coefficients
  summ <- summary(fit)
  r.squared <- summ$r.squared
  tstat <- summ$coefficients[, "t value"]
  Fstat <- summ$fstatistic[1]
  tprob <- rep(1,length(tstat))
  Fprob <- 1

  # Get confidence interval
  conf_df <- stats::confint(fit, names(fit$coefficients), level = 0.90)
  rownames(conf_df) <-  c("Intercept", names(X))

  # Perform permutations
  for(i in 1:nperm){
    rand <- sample(1:nrowsY)
    Yperm <- Y[rand, rand]
    yperm <- unfold(Yperm, scale)
    fit <- stats::lm(yperm ~ Xmats)
    summ <- summary(fit)
    Fprob <- Fprob + as.numeric(summ$fstatistic[1] >= Fstat)
    tprob <- tprob + as.numeric(abs(summ$coefficients[, "t value"]) >= abs(tstat))
  }

  # Return values
  tp <- tprob / (nperm + 1)
  Fp <- Fprob / (nperm + 1)
  names(r.squared) <- "r.squared"
  names(coeffs) <- c("Intercept", names(X))
  names(tstat) <- paste(c("Intercept", names(X)), "(t)", sep="")
  names(tp) <- paste(c("Intercept", names(X)), "(p)", sep="")
  names(Fstat) <- "F-statistic"
  names(Fp) <- "F p-value"

  return(list(r.squared = r.squared,
              coefficients = coeffs,
              tstatistic = tstat,
              tpvalue = tp,
              Fstatistic = Fstat,
              Fpvalue = Fp,
              conf_df = conf_df))
}

#' unfold converts the lower diagonal elements of a matrix into a vector
#' @param X is a distance matrix.
#' @param scale if TRUE then matrices will be standardized. Default = TRUE.

unfold <- function(X, scale = TRUE){
  x <- vector()
  for(i in 2:nrow(X)) x <- c(x, X[i, 1:i-1])
  if(scale == TRUE) x <- raster::scale(x, center = TRUE, scale = TRUE)
  return(x)
}


#' Make nice dataframe from MMRR results
#'
#' @param mod The fitted MMRR model
#'
#' @return
#' @export
#'
#' @examples
mmrr_df <- function(mod){
  coeff_df <- data.frame(coeff = mod$coefficients, p = mod$tpvalue)
  coeff_df$var <- rownames(coeff_df)
  ci_df <- data.frame(mod$conf_df)
  ci_df$var <- rownames(ci_df)
  coeff_df <- merge(coeff_df, ci_df, by = "var")
  rownames(coeff_df) <- NULL
  return(coeff_df)
}

#' Plot MMRR results
#'
#' @param Y The dependent variable in the form of a distance matrix
#' @param X A list of independent variables in the form of distance matrices (required if `plot_type = "fitted"`, `vars` or `"all"`)
#' @param mod The fitted MMRR model (required if `plot_type = "fitted"` or `"all"`)
#' @param plot_type which plots to produce (options: (1) "vars" to plot single variable relationships, (2) "fitted" to plot the fitted relationship, (3) "cov" to plot covariances between the predictor variables, (4) "all" to produce all plots (default))
#' @param scale If TRUE, all variables are scaled
#' @param varNames A vector of names for the variables in the model (optional)
#'
#' @return
#' @export
#'
#' @examples
plot_mmrr <- function(Y, X = NULL, mod = NULL, plot_type = "all", stdz = TRUE, var_names = NULL){

  # Plot single variable relationships
  if("all" %in% plot_type | "vars" %in% plot_type) print(mmrr_plot_vars(Y, X, stdz = TRUE))

  # Plot fitted relationship
  if("all" %in% plot_type | "fitted" %in% plot_type) print(mmrr_plot_fitted(mod, Y, X, stdz = TRUE))

  # Plot fitted relationship
  if("all" %in% plot_type | "cov" %in% plot_type) print(mmrr_plot_cov(X, stdz = TRUE))

  return()
}

#' Plot single variable relationships
#'
#' @inheritParams plot_mmrr
#'
#' @export
#' @noRd
mmrr_plot_vars <- function(Y, X, stdz = TRUE){
  # Unfold X and Y
  y <- unfold(Y, scale = stdz)
  dfX <- purrr::map_dfc(X, unfold, scale = stdz) %>% purrr::map_dfc(as.numeric)

  # Make single variable dataframe
  df <- dfX %>%
    dplyr::mutate(Y = y) %>%
    tidyr::gather("var", "X", -Y)

  # Plot single variable relationships
  plt_lm <- ggplot2::ggplot(df, ggplot2::aes(X, Y)) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
    ggplot2::facet_wrap(~var, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

  return(plt_lm)
}

#' Plot fitted relationship
#'
#' @inheritParams plot_mmrr
#'
#' @export
#' @noRd
mmrr_plot_fitted <- function(mod, Y, X, stdz = TRUE){

  # Make model dataframe
  coeff_df <- mmrr_df(mod)

  # Make fitted dataframe
  df_fitted <- purrr::map_dfc(X, unfold, scale = stdz) %>%
    purrr::map_dfc(as.numeric) %>%
    dplyr::mutate(Y = unfold(Y, scale = stdz)) %>%
    tidyr::gather("var", "X", -Y) %>%
    dplyr::left_join(coeff_df, by = "var") %>%
    dplyr::mutate(coeffX = coeff*X) %>%
    dplyr::select(Y, coeffX) %>%
    dplyr::group_by(Y) %>%
    dplyr::summarise(Yfitted = sum(coeffX))

  # Plot fitted relationship
  plt_fitted <- ggplot2::ggplot(data = df_fitted, ggplot2::aes(x = Yfitted, y = Y)) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Observed Genetic Distance") +
    ggplot2::xlab("Predicted Genetic Distance") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

  return(plt_fitted)
}

#' Plot covariances
#'
#' @inheritParams plot_mmrr
#'
#' @export
#' @noRd
mmrr_plot_cov <- function(X, stdz = TRUE){

  # Unfold X
  dfX <- purrr::map_dfc(X, unfold, scale = stdz) %>% purrr::map_dfc(as.numeric)

  # Plot covariances
  plt_cor <- GGally::ggpairs(dfX, progress = FALSE,
                             lower = list(continuous = GGally::wrap("points", col = "#6464c8", alpha = 0.1, cex = 0.9)),
                             diag = list(continuous = GGally::wrap("densityDiag",  fill = "blue", alpha = 0.1)))
  plt_cor <- plt_cor + ggplot2::theme_bw() + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank())

  return(plt_cor)
}
