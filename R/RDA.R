#' RDA function to do everything
#'
#' @param gen genotype dosage matrix (rows = individuals & columns = SNPs) or `vcfR` object
#' @param env dataframe with environmental data or a Raster* type object from which environmental values for the coordinates can be extracted
#' @param coords dataframe with coordinates (only needed if correctGEO = TRUE) or if env is a Raster* from which values should be extracted
#' @param model whether to fit the model with all variables ("full") or to perform variable selection to determine the best set of variables ("best"); defaults to "best"
#' @param correctGEO whether to condition on geographic coordinates
#' @param correctPC whether to condition on PCs from PCA of genotypes
#' @param outlier_method method to determine outliers. Can either be "p" to use the p-value method from [here](https://github.com/Capblancq/RDA-landscape-genomics) or "z" to use the z-score based method from [here](https://popgen.nescent.org/2018-03-27_RDA_GEA.html)
#' @param sig if `outlier_method = "p"`, the significance level to use to identify SNPs (defaults to 0.05)
#' @param p_adj if `outlier_method = "p"`, method to use for p-value correction (defaults to "fdr"); other options can be found in [p.adjust()]
#' @param z if `outlier_method = "z"`, the number of standard deviations to use to identify SNPs (defaults to 3)
#' @param cortest whether to create table of correlations for SNPs and environmental variable (defaults to TRUE)
#' @param nPC number of PCs to use if correctPC = TRUE (defaults to 3); if set to "manual" a selection option with a terminal prompt will be provided
#' @param varpart whether to perform variance partitioning (defaults to FALSE)
#' @param naxes number of RDA axes to use (defaults to "all" to use all axes), if set to "manual" a selection option with a terminal prompt will be given, otherwise can be any integer that is less than or equal to the total number of axes
#' @param Pin if `model = "best"`, limits of permutation P-values for adding (`Pin`) a term to the model, or dropping (`Pout`) from the model. Term is added if` P <= Pin`, and removed if `P > Pout` (see \link[vegan]{ordiR2step})
#' @param R2permutations if `model = "best"`, number of permutations used in the estimation of adjusted R2 for cca using RsquareAdj (see \link[vegan]{ordiR2step})
#' @param R2scope if `model = "best"`, use adjusted R2 as the stopping criterion: only models with lower adjusted R2 than scope are accepted (see \link[vegan]{ordiR2step})
#' @param stdz whether to center and scale environmental data (defaults to TRUE)
#' @param quiet whether to print output tables and figures (defaults to FALSE)
#'
#' @inheritParams vegan::ordiR2step
#'
#' @importFrom vegan rda
#'
#' @return list containing (1) outlier SNPs, (2) dataframe with correlation test results, if `cortest = TRUE`, (3) the RDA model, (4) results from outlier analysis (output from \link[algatr]{rda_getoutliers}), (5) RDA R-Squared, (6) RDA ANOVA, (7) p-values if `outlier_method = "p"`, and (8) results from variance partitioning analysis, if `varpart = TRUE`
#' @export
#' @details
#' Much of algatr's code is adapted from Capblancq T., Forester B.R. 2021. Redundancy analysis: A swiss army knife for landscape genomics. Methods Ecol. Evol. 12:2298-2309. doi: https://doi.org/10.1111/2041-210X.13722.
#'
#' @family RDA functions
#'
#' @examples
rda_do_everything <- function(gen, env, coords = NULL, model = "best", correctGEO = FALSE, correctPC = FALSE,
                              outlier_method = "p", sig = 0.05, z = 3,
                              p_adj = "fdr", cortest = TRUE, nPC = 3, varpart = FALSE, naxes = "all",
                              Pin = 0.05, R2permutations = 1000, R2scope = T, stdz = TRUE, quiet = FALSE) {
  # Modify environmental data --------------------------------------------------------------------------------------------------

  # Extract environmental data if env is a raster
  if (inherits(env, "Raster")) env <- terra::rast(env)
  if (inherits(env, "SpatRaster")) crs_check(coords = coords, lyr = env)
  if (inherits(env, "SpatRaster")) env <- terra::extract(env, coords_to_sf(coords), ID = FALSE)

  # Standardize environmental variables
  if (stdz) env <- terra::scale(env, center = TRUE, scale = TRUE)
  env <- data.frame(env)

  # Format coords
  if (!is.null(coords)) coords <- coords_to_df(coords)

  # Modify genetic data -----------------------------------------------------

  # Convert vcf to dosage
  if (inherits(gen, "vcfR")) gen <- wingen::vcf_to_dosage(gen)

  # Perform imputation with warning
  if (any(is.na(gen))) {
    gen <- simple_impute(gen, median)
    warning("NAs found in genetic data, imputing to the median (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
  }

  # Check for NAs
  if (any(is.na(gen))) {
    stop("NA values found in gen data")
  }

  if (any(is.na(env))) {
    warning("NA values found in env data, removing rows with NAs for RDA")
    gen <- gen[complete.cases(env), ]
    coords <- coords[complete.cases(env), ]
    # NOTE: this must be last
    env <- env[complete.cases(env), ]
  }

  # Running RDA ----------------------------------------------------------------------------------------------------------------

  # Run model
  mod <- rda_run(gen, env, coords,
    model = model,
    correctGEO = correctGEO,
    correctPC = correctPC,
    nPC = nPC,
    Pin = Pin,
    R2permutations = R2permutations,
    R2scope = R2scope
  )

  # If NULL, exit
  if (is.null(mod)) {
    warning("Model is NULL, returning NULL object")
    return(NULL)
  }

  # get R-squared and run ANOVA
  mod_rsq <- vegan::RsquareAdj(mod)
  mod_aov <- stats::anova(mod)

  # Variance partitioning ---------------------------------------------------

  if (varpart) {
    varpart_df <- rda_varpart(gen, env, coords, Pin = Pin, R2permutations = R2permutations, R2scope = R2scope, nPC = nPC)
    if (!quiet) print(rda_varpart_table(varpart_df))
  } else {
    varpart_df <- NULL
  }

  # Identify candidate SNPs ----------------------------------------------------------------------------------------------------

  # Running with all axes
  rda_sig <- rda_getoutliers(mod, naxes = naxes, outlier_method = outlier_method, p_adj = p_adj, sig = sig)

  # Get SNPs
  rda_snps <- rda_sig$rda_snps

  # Summarize results ---------------------------------------------------------------------------------------------------------------

  # Plot all axes
  if (any("pvalues" %in% names(rda_sig))) pvalues <- rda_sig[["pvalues"]] else pvalues <- NULL
  if (!quiet) rda_plot(mod, rda_snps = rda_snps, pvalues = pvalues, axes = "all", biplot_axes = NULL, sig = sig, manhattan = TRUE, rdaplot = TRUE)

  # Get correlations -----------------------------------------------------------------------------------------------------------
  rda_gen <- gen[, rda_snps]
  if (cortest) {
    cor_df <- rda_cor(rda_gen, env)
    if (!quiet) print(rda_table(cor_df, top = TRUE, order = TRUE, nrow = 10))
  } else {
    cor_df <- NULL
  }

  # Compile results ------------------------------------------------------------------------------------------------------------

  results <- list(
    rda_snps = rda_snps,
    cor_df = cor_df,
    rda_mod = mod,
    rda_outlier_test = rda_sig,
    rsq = mod_rsq,
    anova = mod_aov,
    pvalues = pvalues,
    varpart = varpart_df
  )

  return(results)
}


#' Run RDA
#'
#' @inheritParams rda_doEverything
#'
#' @return RDA model
#' @export
#'
#' @family RDA functions
#'
rda_run <- function(gen, env, coords = NULL, model = "full",
                    correctGEO = FALSE, correctPC = FALSE, nPC = 3,
                    Pin = 0.05, R2permutations = 1000, R2scope = T) {
  if (!correctPC & !correctGEO) {
    moddf <- data.frame(env)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+")))
  }

  if (correctPC & !correctGEO) {
    pcres <- stats::prcomp(gen)
    stats::screeplot(pcres, type = "barplot", npcs = length(pcres$sdev), main = "PCA Eigenvalues")
    if (nPC == "manual") nPC <- readline("Number of PC axes to retain:")
    pc <- pcres$x[, 1:nPC]
    moddf <- data.frame(env, pc)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+"), "+ Condition(", paste(colnames(pc), collapse = "+"), ")"))
  }

  if (!correctPC & correctGEO) {
    if (is.null(coords)) stop("Coordinates must be provided if correctGEO is TRUE")
    moddf <- data.frame(env, coords)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+"), "+ Condition(x + y)"))
  }

  if (correctPC & correctGEO) {
    if (is.null(coords)) stop("Coordinates must be provided if correctGEO is TRUE")
    pcres <- stats::prcomp(gen)
    stats::screeplot(pcres, type = "barplot", npcs = length(pcres$sdev), main = "PCA Eigenvalues")
    if (nPC == "manual") nPC <- readline("Number of PC axes to retain:")
    pc <- pcres$x[, 1:nPC]
    moddf <- data.frame(env, coords, pc)
    f <- as.formula(paste0("gen ~ ", paste(colnames(env), collapse = "+"), "+ Condition(", paste(colnames(pc), collapse = "+"), "+ x + y)"))
  }

  if (model == "best") {
    mod_full <- vegan::rda(f, data = moddf)
    mod_null <- vegan::rda(gen ~ 1, data = moddf)
    mod <- vegan::ordiR2step(mod_null, mod_full, Pin = Pin, R2permutations = R2permutations, R2scope = R2scope)
    if (mod$call == mod_null$call) {
      mod <- NULL
      warning("Best model is NULL model, returning NULL")
    }
  } else {
    mod <- vegan::rda(f, data = moddf)
  }

  return(mod)
}



#' Get significant outliers from RDA model
#'
#' @param plot whether to produce scree plot of RDA axes (defaults to TRUE)
#' @inheritParams rda_do_everything
#'
#' @return results from outlier tests. If `outlier_method = "p"`, a list of outlier SNPs, p-values, and results from rdadapt (see [Capblancq et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12906)). If `outlier_method = "z"`, a dataframe with outlier SNP Z-scores for each axis
#' @export
#'
#' @family RDA functions
#'
rda_getoutliers <- function(mod, naxes = "all", outlier_method = "p", p_adj = "fdr", sig = 0.05, z = 3, plot = TRUE) {
  # Running the function with all axes
  if (plot) stats::screeplot(mod, main = "Eigenvalues of constrained axes")
  if (naxes == "manual") naxes <- readline("Number of RDA axes to retain:")
  if (naxes == "all") naxes <- ncol(mod$CCA$v)

  if (outlier_method == "p" & naxes == 1) warning("Cannot compute p-values (outlier_method = \"p\") when the number of RDA axes is less than two, using the standard deviation based method (outlier_method = \"z\") instead")
  if (outlier_method == "p" & naxes != 1) results <- p_outlier_method(mod, naxes, sig, p_adj)
  if (outlier_method == "z" | naxes == 1) results <- z_outlier_method(mod, naxes, z)

  return(results)
}



#' Determine RDA outliers based on p-values
#'
#' @inheritParams rda_do_everything
#' @export
#' @noRd
#'
#' @family RDA functions
#'
p_outlier_method <- function(mod, naxes, sig = 0.05, p_adj = "fdr") {
  rdadapt_env <- rdadapt(mod, naxes)

  # p-value threshold after p-value adjustment (different from Capblancq & Forester 2021)
  pvalues <- p.adjust(rdadapt_env$p.values, method = p_adj)

  # Get SNP names
  snp_names <- rownames(vegan::scores(mod, choices = naxes, display = "species"))

  # Restore SNP names
  names(pvalues) <- snp_names

  # Identifying the SNPs that are below the p-value threshold
  rda_snps <- snp_names[which(pvalues < sig)]
  if (length(rda_snps) == 0) {
    warning("No significant SNPs found, returning NULL object")
    return(NULL)
  }

  results <- list(
    rda_snps = rda_snps,
    pvalues = pvalues,
    rdadapt = rdadapt_env
  )

  return(results)
}

#' Determine RDA outliers based on Z-scores
#'
#' @inheritParams rda_do_everything
#'
#' @export
#' @noRd
#'
#' @family RDA functions
#'
z_outlier_method <- function(mod, naxes, z = 3) {
  load.rda <- vegan::scores(mod, choices = naxes, display = "species")

  results <- purrr::map_dfr(data.frame(1:ncol(load.rda)), z_outlier_helper, load.rda, z)

  return(results)
}

#' z_outlier_method helper function
#'
#' @export
#' @noRd
#'
#' @family RDA functions
#'
z_outlier_helper <- function(axis, load.rda, z) {
  x <- load.rda[, axis]
  out <- outliers(x, z)
  cand <- cbind.data.frame(names(out), rep(axis, times = length(out)), unname(out))
  colnames(cand) <- c("rda_snps", "axis", "loading")
  cand$rda_snps <- as.character(cand$rda_snps)
  return(cand)
}

#' Z-scores outlier finder
#'
#' @details code adapted from [Forester et al. 2018](https://popgen.nescent.org/2018-03-27_RDA_GEA.html)
#'
#' @export
#' @noRd
#'
#' @family RDA functions
#'
outliers <- function(x, z) {
  lims <- mean(x) + c(-1, 1) * z * sd(x) # find loadings +/-z sd from mean loading
  x[x < lims[1] | x > lims[2]] # SNP names in these tails
}

#' Function to conduct a RDA-based genome scan
#'
#' @param mod model object of class `rda`
#' @param K number of RDA axes to retain when detecting outliers
#'
#' @details Method developed by [Capblancq et al. 2018](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12906)
#' Code provided in [Capblancq & Forester 2021](https://github.com/Capblancq/RDA-landscape-genomics/blob/main/RDA_landscape_genomics.Rmd)
#'
#' @export
#' @noRd
#' @family RDA functions
rdadapt <- function(mod, K) {
  # Extract scores based on number of specified RDA axes
  zscores <- mod$CCA$v[, 1:as.numeric(K)]
  # Standardize by scaling
  resscale <- apply(zscores, 2, scale)
  # Calculate squared Mahalanobis distances for each locus
  resmaha <- robust::covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  # Distribute Mahalanobis distances as chi-sq dist'n with K DF; calculate genomic inflation factor (lambda)
  lambda <- median(resmaha) / qchisq(0.5, df = K)
  # Rescale distances according to lambda (genomic inflation factor); resulting values are p-values
  reschi2test <- pchisq(resmaha / lambda, K, lower.tail = FALSE)
  # Obtain q-values
  qval <- qvalue::qvalue(reschi2test)
  q.values_rdadapt <- qval$qvalues

  return(data.frame(p.values = reschi2test, q.values = q.values_rdadapt))
}

#' Genotype-environment correlation test
#'
#' @param gen dosage matrix
#' @param var dataframe with predictor variables
#'
#' @return dataframe with r and p-values from correlation test
#' @export
#' @family RDA functions
rda_cor <- function(gen, var) {
  cor_df <- purrr::map_dfr(colnames(gen), rda_cor_env_helper, gen, var)
  rownames(cor_df) <- NULL
  colnames(cor_df) <- c("r", "p", "snp", "var")
  return(cor_df)
}

#' Helper function for rda_cor_test
#'
#' @export
#' @noRd
#' @family RDA functions
rda_cor_env_helper <- function(snp_name, snp_df, env) {
  cor_df <- data.frame(t(apply(env, 2, rda_cor_helper, snp_df[, snp_name])))
  cor_df$snp <- snp_name
  cor_df$env <- colnames(env)
  return(cor_df)
}

#' Helper function for rda_cor_test
#'
#' @export
#' @noRd
#' @family RDA functions
rda_cor_helper <- function(envvar, snp) {
  if (sum(!is.na(envvar)) < 3 | sum(!is.na(snp)) < 3) {
    return(c(r = NA, p = NA))
  }
  # kendall is used instead of pearson because it is non-parameteric and doesn't require vars to be continuous
  mod <- stats::cor.test(envvar, snp, alternative = "two.sided", method = "kendall", na.action = "na.omit")
  pvalue <- mod$p.value
  r <- mod$estimate
  results <- c(r, pvalue)
  names(results) <- c("r", "p")
  return(results)
}

#' Plot RDA results
#'
#' @param mod model object of class `rda`; if this is all that's provided, histograms with loadings will be generated
#' @param rda_snps vector of outlier SNPs (defaults to NULL)
#' @param pvalues if creating a Manhattan plot (i.e., `manhattan = TRUE`), a matrix of p-values (defaults to NULL)
#' @param axes which RDA axes to include while plotting (defaults to `"all"`)
#' @param biplot_axes if creating an RDA biplot (i.e., `rdaplot = TRUE`), which pairs of axes to plot. Defaults to plotting all pairs of axes possible, otherwise can be set to a single pair of axes (e.g., c(1,2)) or a list of axes pairs (e.g., list(c(1,2), c(2,3))))
#' @param manhattan whether to produce Manhattan plot (defaults to `TRUE`)
#' @param rdaplot whether to produce an RDA biplot (defaults to `TRUE`). If only one axis is provided, instead of a biplot, a histogram will be created
#' @param sig if creating a Manhattan plot, significance threshold for y axis (defaults to 0.05)
#' @param binwidth width of bins for histograms (defaults to NULL)
#'
#' @export
#'
#' @family RDA functions
#'
rda_plot <- function(mod, rda_snps = NULL, pvalues = NULL, axes = "all", biplot_axes = NULL, sig = 0.05, manhattan = NULL, rdaplot = NULL, binwidth = NULL) {
  # Get axes
  if (axes == "all") axes <- 1:ncol(mod$CCA$v)

  # Histograms with loadings ------------------------------------------------

  if (is.null(rda_snps)) {
    # Extract loadings from model (RDA places SNP names within "species")
    loadings <- vegan::scores(mod, choices = axes, display = "species")

    # Tidy data
    loadings <- loadings %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "SNP") %>%
      tidyr::pivot_longer(!SNP, names_to = "axis", values_to = "loading")

    # Generate plot, faceting on RDA axis
    print(rda_hist(loadings, binwidth = binwidth))
  }


  # If outliers found -------------------------------------------------------

  if (!is.null(rda_snps)) {
    # Make and get tidy dataframes for plotting
    tidy_list <- rda_ggtidy(mod, rda_snps, axes = axes)
    TAB_snps <- tidy_list[["TAB_snps"]]
    TAB_var <- tidy_list[["TAB_var"]]

    # Make RDA plots
    if (rdaplot) {
      if (length(axes) == 1) {
        print(rda_hist(TAB_snps, binwidth = binwidth))
      } else if (!is.null(biplot_axes)) {
        if (is.vector(biplot_axes)) print(rda_biplot(TAB_snps, TAB_var, biplot_axes = biplot_axes))
        if (is.list(biplot_axes)) {
          lapply(biplot_axes, function(x) {
            print(rda_biplot(TAB_snps, TAB_var, biplot_axes = x))
          })
        }
      } else {
        cb <- combn(length(axes), 2)
        if (!is.null(dim(cb))) {
          apply(cb, 2, function(x) {
            print(rda_biplot(TAB_snps, TAB_var, biplot_axes = x))
          })
        } else {
          print(rda_biplot(TAB_snps, TAB_var, biplot_axes = cb))
        }
      }
    }

    # Make Manhattan plot
    if (manhattan & !is.null(pvalues)) print(rda_manhattan(TAB_snps, rda_snps, pvalues, sig = sig))
  }
}


#' Make dataframe for ggplot from RDA results
#'
#' @export
#' @noRd
#' @family RDA functions
rda_ggtidy <- function(mod, rda_snps, axes) {
  snp_scores <- vegan::scores(mod, choices = axes, display = "species", scaling = "none") # vegan references "species", here these are the snps
  TAB_snps <- data.frame(names = row.names(snp_scores), snp_scores)

  TAB_snps$type <- "Neutral"
  TAB_snps$type[TAB_snps$names %in% rda_snps] <- "Outliers"
  TAB_snps$type <- factor(TAB_snps$type, levels = c("Neutral", "Outliers"))
  TAB_var <- as.data.frame(vegan::scores(mod, choices = axes, display = "bp")) # pull the biplot scores

  tidy_list <- list(TAB_snps = TAB_snps, TAB_var = TAB_var)
  return(tidy_list)
}


#' Helper function to plot RDA biplot
#'
#' @export
#' @noRd
#' @family RDA functions
rda_biplot <- function(TAB_snps, TAB_var, biplot_axes = c(1, 2)) {
  # Select axes for plotting
  xax <- paste0("RDA", biplot_axes[1])
  yax <- paste0("RDA", biplot_axes[2])
  TAB_snps_sub <- TAB_snps[, c(xax, yax, "type")]
  colnames(TAB_snps_sub) <- c("x", "y", "type")
  TAB_var_sub <- TAB_var[, c(xax, yax)]
  colnames(TAB_var_sub) <- c("x", "y")

  # Scale the variable loadings for the arrows
  TAB_var_sub$x <- TAB_var_sub$x * max(TAB_snps_sub$x) / stats::quantile(TAB_var_sub$x)[4]
  TAB_var_sub$y <- TAB_var_sub$y * max(TAB_snps_sub$y) / stats::quantile(TAB_var_sub$y)[4]

  ## Biplot of RDA SNPs and scores for variables
  ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = gray(.80), size = 0.6) +
    ggplot2::geom_point(data = TAB_snps_sub, ggplot2::aes(x = x, y = y, colour = type), size = 1.4) +
    ggplot2::scale_color_manual(values = c(rgb(0.7, 0.7, 0.7, 0.1), "#F9A242FF")) +
    ggplot2::geom_segment(data = TAB_var_sub, ggplot2::aes(xend = x, yend = y, x = 0, y = 0), colour = "black", size = 0.15, linetype = 1, arrow = ggplot2::arrow(length = ggplot2::unit(0.02, "npc"))) +
    ggrepel::geom_text_repel(data = TAB_var_sub, ggplot2::aes(x = x, y = y, label = row.names(TAB_var_sub)), size = 4) +
    ggplot2::xlab(xax) +
    ggplot2::ylab(yax) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "snp type")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = ggplot2::rel(.8)),
      strip.text = ggplot2::element_text(size = 11)
    )
}

#' Helper function to plot RDA manhattan plot
#'
#' @export
#' @noRd
#' @family RDA functions
rda_manhattan <- function(TAB_snps, rda_snps, pvalues, sig = 0.05) {
  TAB_manhattan <- data.frame(
    pos = 1:nrow(TAB_snps),
    pvalues = pvalues,
    type = factor(TAB_snps$type, levels = c("Neutral", "Outliers"))
  )

  TAB_manhattan <- TAB_manhattan[order(TAB_manhattan$pos), ]

  ggplot2::ggplot(data = TAB_manhattan) +
    ggplot2::geom_point(ggplot2::aes(x = pos, y = -log10(pvalues), col = type), size = 1.4) +
    ggplot2::scale_color_manual(values = c(rgb(0.7, 0.7, 0.7, 0.5), "#F9A242FF", "#6B4596FF")) +
    ggplot2::xlab("position") +
    ggplot2::ylab("-log10(p)") +
    ggplot2::geom_hline(yintercept = -log10(sig), linetype = "dashed", color = "black", size = 0.6) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "SNP type")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      legend.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = ggplot2::rel(.8)),
      strip.text = ggplot2::element_text(size = 11)
    )
}

#' Helper function to plot RDA histogram
#'
#' @param data TAB_snps if `rda_getoutliers()` has been run, otherwise can be loadings
#' @param binwidth width of bins for histogram
#'
#' @export
#' @noRd
#' @family RDA functions
rda_hist <- function(data, binwidth = NULL) {
  if ("type" %in% names(data)) {
    ggplot2::ggplot() +
      ggplot2::geom_histogram(data = data, ggplot2::aes(fill = type, x = get(colnames(data)[2])), binwidth = binwidth) +
      ggplot2::scale_fill_manual(values = c(rgb(0.7, 0.7, 0.7, 0.5), "#F9A242FF")) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "SNP type")) +
      ggplot2::xlab(colnames(data)[2]) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = "right",
        legend.background = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        legend.box.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        legend.text = ggplot2::element_text(size = ggplot2::rel(.8)),
        strip.text = ggplot2::element_text(size = 11)
      )
  } else {
    ggplot2::ggplot() +
      ggplot2::geom_histogram(data = data, ggplot2::aes(x = loading), bins = binwidth) +
      ggplot2::facet_wrap(~axis) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 11)
      )
  }
}

#' Create `gt` table of RDA results
#'
#' @param cor_df dataframe of correlation results output from \link[algatr]{rda_cor}
#' @param sig_only whether to only include loci with p-values less than `sig` (defaults to TRUE)
#' @param top whether to only include only keep the top variable for each snp in the table by the strength of the correlation (defaults to FALSE)
#' @param order whether to order by the magnitude of the correlation (defaults to FALSE)
#' @param var which variables to include (defaults to including all variables)
#' @param nrow number of rows to display (defaults to displaying all rows)
#' @param digits number of digits to include (defaults to 2)
#' @inheritParams rda_do_everything
#'
#' @return An object of class `gt_tbl`
#' @export
#' @family RDA functions
rda_table <- function(cor_df, sig = 0.05, sig_only = TRUE, top = FALSE, order = FALSE, var = NULL, nrow = NULL, digits = 2) {
  if (!is.null(var)) cor_df <- cor_df[cor_df$var %in% var, ]
  if (sig_only) cor_df <- cor_df[cor_df$p < sig, ]

  if (nrow(cor_df) == 0) {
    warning("no significant variants found, returning NULL object")
    return(NULL)
  }

  if (order) cor_df <- cor_df[order(abs(cor_df$r), decreasing = TRUE), ]
  if (top) {
    cor_df <- cor_df %>%
      dplyr::group_by(snp) %>%
      dplyr::filter(abs(r) == max(abs(r)))
  }
  if (!is.null(nrow)) {
    if (nrow > nrow(cor_df)) nrow <- nrow(cor_df)
    cor_df <- cor_df[1:nrow, ]
  }

  cor_df <- cor_df %>% dplyr::as_tibble()
  if (!is.null(digits)) cor_df <- cor_df %>% dplyr::mutate(dplyr::across(-c(var, snp), round, digits))

  d <- max(abs(min(cor_df$r)), abs(max(cor_df$r)))

  suppressWarnings(
    tbl <- cor_df %>%
      gt::gt() %>%
      gtExtras::gt_hulk_col_numeric(r, trim = TRUE, domain = c(-d, d))
  )

  tbl
}

#' Partial RDA variance partitioning
#'
#' @inheritParams rda_doEverything
#'
#' @return df with relevant statistics from variance partitioning analysis
#' @export
#'
#' @family RDA functions
#'
#' @examples
rda_varpart <- function(gen, env, coords, Pin, R2permutations, R2scope, nPC) {
  moddf <- data.frame(env)

  # Run best ----------------------------------------------------------------

  f <- as.formula(paste0("gen ~ ", paste(colnames(moddf), collapse = "+")))
  mod_best <- rda_run(gen, env,
    model = "best",
    Pin = Pin,
    R2permutations = R2permutations,
    R2scope = R2scope
  )

  # Extract sig enviro vars
  sig_vars <- as.character(mod_best$terms)[3]
  if (length(sig_vars) == 0) {
    warning("No significant terms found in best model, returning NULL object")
    return(NULL)
  }

  # Run PCA for pop structure -----------------------------------------------

  pcres <- stats::prcomp(gen)
  stats::screeplot(pcres, type = "barplot", npcs = length(pcres$sdev), main = "PCA Eigenvalues")
  if (nPC == "manual") nPC <- readline("Number of PC axes to retain:")
  pc <- pcres$x[, 1:nPC]

  moddf_covar <- data.frame(env, coords, pc)

  # Run RDAs ----------------------------------------------------------------

  # Full model with covariables as full expl vars; only sig enviro vars
  f <- as.formula(paste0("gen ~ ", paste(sig_vars), " + ", paste(colnames(pc), collapse = "+"), "+ x + y"))
  full <- vegan::rda(f, data = moddf_covar)

  # Pure env
  f <- as.formula(paste0("gen ~ ", paste(sig_vars), " + Condition(", paste(colnames(pc), collapse = "+"), "+ x + y)"))
  pure_env <- vegan::rda(f, data = moddf_covar)

  # Pure structure
  f <- as.formula(paste0("gen ~ ", paste(colnames(pc), collapse = "+"), "+ Condition(", paste(sig_vars), "+ x + y)"))
  pure_str <- vegan::rda(f, data = moddf_covar)

  # Pure geo
  f <- as.formula(paste0("gen ~ x + y + Condition(", paste(colnames(pc), collapse = " + "), " + ", paste(sig_vars), ")"))
  pure_geo <- vegan::rda(f, data = moddf_covar)

  # Run helper function on models -------------------------------------------

  df <- rbind(
    rda_varpart_helper(full),
    rda_varpart_helper(pure_env),
    rda_varpart_helper(pure_str),
    rda_varpart_helper(pure_geo)
  )

  # Calculate relevant stats ------------------------------------------------

  total_inertia <- mod_best$tot.chi
  full_inertia <- df$inertia[1]
  confounded <- as.numeric(
    full_inertia - (df %>%
      dplyr::filter(rownames(df) %in% c("pure_env", "pure_str", "pure_geo")) %>%
      dplyr::summarise(sum(inertia)))
  )
  total_unexpl <- total_inertia - full_inertia
  results <- data.frame(total_inertia, full_inertia, confounded, total_unexpl)

  # Compile df --------------------------------------------------------------

  df <- df %>%
    dplyr::mutate(
      prop_expl_var = inertia / max(df$inertia),
      prop_total_var = inertia / total_inertia
    )

  # Add additional rows
  df <-
    df %>%
    tibble::add_row(
      inertia = results$confounded,
      prop_expl_var = results$confounded / results$full_inertia,
      prop_total_var = results$confounded / results$total_inertia
    ) %>%
    tibble::add_row(
      inertia = results$total_unexpl,
      prop_total_var = results$total_unexpl / results$total_inertia
    ) %>%
    tibble::add_row(
      inertia = results$total_inertia,
      prop_total_var = 1
    )

  rownames(df) <- c("full", "pure_env", "pure_str", "pure_geo", "confounded", "total_unexplained", "total")

  return(df)
}



#' Helper function for `rda_varpart()`
#'
#' Extracts relevant statistics from variance partitioning analysis
#'
#' @param mod RDA model results
#'
#' @return df with relevant statistics (call, R2, adjusted R2, inertia)
#' @export
#'
#' @noRd
#' @family RDA functions
#' @examples
rda_varpart_helper <- function(mod) {
  call <- paste(mod$call)[2]
  R2adj <- vegan::RsquareAdj(mod)
  results <- anova(mod)
  inertia <- results$Variance[1]
  p <- results$`Pr(>F)`[1]
  df <- data.frame(call, inertia, R2adj, p)
  rownames(df) <- deparse(substitute(mod))

  return(df)
}

#' Create `gt` table with RDA variance partitioning results
#'
#' @param df dataframe of variance partitioning results output from \link[algatr]{rda_varpart}
#' @param digits number of digits to include (defaults to 2)
#' @param call_col whether to include column with RDA call (defaults to FALSE)
#'
#' @return
#' @export
#'
#' @family RDA functions
#'
#' @examples
rda_varpart_table <- function(df, digits = 2, call_col = FALSE) {
  # Replace row and column names
  rownames(df) <- c("Full model", "Pure enviro. model", "Pure pop. structure model", "Pure geography model", "Confounded variance", "Total unexplained variance", "Total inertia")
  colnames(df) <- c(
    "pRDA model call", "Inertia", "R2", "Adjusted R2", "p (>F)",
    "Prop. of explainable variance", "Prop. of total variance"
  )

  df <- df %>%
    tibble::rownames_to_column(var = "Model")

  if (!is.null(digits)) df <- df %>% dplyr::mutate(dplyr::across(-c("Model", "pRDA model call"), ~ round(.x, digits)))

  d <- max(abs(min(df$Inertia, na.rm = TRUE)), abs(max(df$Inertia, na.rm = TRUE)))

  suppressWarnings({
    tbl <- df %>%
      dplyr::select(-"pRDA model call") %>%
      gt::gt() %>%
      gtExtras::gt_hulk_col_numeric("Inertia", trim = TRUE, domain = c(-d, d)) %>%
      gt::sub_missing(missing_text = "") %>%
      gt::tab_header(title = "Variance partitioning")

    # Add column with call
    if (call_col) {
      tbl <- df %>%
        gt::gt() %>%
        gtExtras::gt_hulk_col_numeric("Inertia", trim = TRUE, domain = c(-d, d)) %>%
        gt::sub_missing(missing_text = "") %>%
        gt::tab_header(title = "Variance partitioning")
    }
  })

  tbl
}
