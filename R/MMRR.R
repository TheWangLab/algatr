#' MMRR function to do everything ()
#'
#' @param gendist matrix of genetic distances
#' @param coords dataframe with x and y coordinates
#' @param env dataframe with environmental values for each coordinate; if not provided it will be calculated based on coords/envlayers
#' @param envlayers rasters for for extracting environmental values using coordinates if `env` isn't provided
#' @param model whether to fit the model with all variables (`"full"`) or to perform variable selection to determine the best set of variables (`"best"`); default = "best"
#' @param nperm number of permutations to use to calculate variable importance; only used if `model = "best"` (default = 999)
#' @param stdz if TRUE then matrices will be standardized (defaults to TRUE)
#' @param plot whether to plot results (defaults to TRUE)
#' @param plot_type which plots to produce (options: (1) "vars" to plot single variable relationships, (2) "fitted" to plot the fitted relationship, (3) "cov" to plot covariances between the predictor variables, (4) "all" to produce all plots (default))
#'
#' @return
#' @export
#'
#' @examples
mmrr_do_everything <- function(gendist, coords, env = NULL, envlayers = NULL, model = "best", nperm = 999, stdz = TRUE, plot = TRUE, plot_type = "all"){

  # If not provided, make env data frame from layers and coords
  if(is.null(env)){env <- raster::extract(envlayers, coords)}

  # Make env dist matrix
  X <- env_dist(env)
  X[["geodist"]] <- geo_dist(coords)

  # Make geodist mat
  Y <- as.matrix(gendist)

  # Run MMRR
  if(model == "best") results <- mmrr_best(Y, X, nperm = nperm, stdz = stdz, plot = plot, plot_type = plot_type)

  if(model == "full") results <- mmrr_full(Y, X, nperm = nperm, stdz = stdz, plot = plot, plot_type = plot_type)

  return(results)

}

#' Run MMRR with variable selection
#'
#' @param Y dependent distance matrix
#' @param X list of independent distance matrices (with optional names)
#' @inheritParams mmrr_do_everything
#'
#' @return
#' @export
#'
#' @examples
mmrr_best <- function(Y, X, nperm = 999, stdz = TRUE, plot = TRUE, plot_type = "all"){

  # fil model with variable selection
  mod <- mmrr_var_sel(Y, X, nperm = nperm, scale = stdz)

  # If NULL, exit with NULL
  if(is.null(mod)) return(NULL)

  # subset X with significant variables
  X_best <- X[names(mod$coefficients)[-1]]
  # plot results
  if (plot) mmrr_plot(Y = Y, X = X_best, mod = mod, plot_type = plot_type, stdz = stdz)

  # Make nice data frame
  coeff_df <- mmrr_df(mod)

  # Make results list
  results <- list(coeff_df = coeff_df,
                  mod = mod,
                  Y = Y,
                  X = X,
                  X_best = X_best)

  return(results)
}

#' Run MMRR with all variables
#'
#' @param Y dependent distance matrix
#' @param X list of independent distance matrices (with optional names)
#' @inheritParams mmrr_do_everything
#'
#' @return
#' @export
#'
#' @examples
mmrr_full <- function(Y, X, nperm = nperm, stdz = TRUE, plot = TRUE, plot_type = "all"){

  # Run full model
  mod <- MMRR(Y, X, nperm = nperm, scale = stdz)

  # If NULL, exit with NULL
  if (is.null(mod)) return(NULL)

  # plot results
  if (plot) mmrr_plot(Y = Y, X = X, mod = mod, plot_type = plot_type, stdz = stdz)

  # Make nice data frame
  coeff_df <- mmrr_df(mod)

  # Make results list
  results <- list(coeff_df = coeff_df,
                  mod = mod,
                  Y = Y,
                  X = X)

  return(results)
}

#' mmrr_var_sel performs MMRR with backward elimination variable selection
#' @param Y is a dependent distance matrix
#' @param X is a list of independent distance matrices (with optional names)
#' @inheritParams mmrr_do_everything
#'
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

  # Repetitive test for no significant variables, but just to be sure
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
mmrr_plot <- function(Y = NULL, X, mod = NULL, plot_type = "all", stdz = TRUE, var_names = NULL){

  # Plot single variable relationships
  if("all" %in% plot_type | "vars" %in% plot_type) print(mmrr_plot_vars(Y, X, stdz = TRUE))

  # Plot fitted relationship
  if("all" %in% plot_type | "fitted" %in% plot_type) print(mmrr_plot_fitted(mod, Y, X, stdz = TRUE))

  # Plot fitted relationship
  if("all" %in% plot_type | "cov" %in% plot_type) print(mmrr_plot_cov(X, stdz = TRUE))
}

#' Plot single variable relationships
#'
#' @inheritParams mmrr_plot
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
#' @inheritParams mmrr_plot
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
#' @inheritParams mmrr_plot
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
