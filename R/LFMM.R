
#' Run LFMM
#'
#' @param gen genotype dosage matrix (rows = individuals & columns = snps) or `vcfR` object
#' @param env dataframe with environmental data or a Raster* type object from which environmental values for the coordinates can be extracted
#' @param coords dataframe with coordinates (only needed if K selection is performed with TESS or if environmental values aren't provided)
#' @param K number of latent factors (if left as NULL (default), K value selection will be conducted)
#' @param lfmm_method lfmm method (either \code{"ridge"} (default) or \code{"lasso"})
#' @param K_selection method for performing k selection (can either by "tracy_widom" (default), "quick_elbow", "tess", or "find_clusters")
#' @param sig alpha level for determining candidate SNPs (defaults to 0.05)
#' @param p_adj method to use for p-value correction (defaults to "fdr")
#' @inheritParams lfmm::lfmm_test
#' @inheritParams select_K
#'
#' @return
#' @export
#'
#' @examples
lfmm_do_everything <- function(gen, env, coords = NULL, K = NULL, lfmm_method = "ridge",
                     K_selection = "tracy_widom", Kvals = 1:10, sig = 0.05,
                     p_adj = "fdr", calibrate = "gif", criticalpoint = 2.0234,
                     low = 0.08, max.pc = 0.9, pca.select = "percVar", perc.pca = 90,
                     choose.n.clust = FALSE, criterion = "diffNgroup", max.n.clust = 10){

  # Get environmental data
  if(inherits(env, "Raster")) env <- raster::extract(env, coords)

  # Convert vcf to dosage matrix
  if(inherits(gen, "vcfR")) gen <- vcf_to_dosage(gen)

  # Perform imputation with warning
  if(any(is.na(gen))){
    gen <- simple_impute(gen, median)
    warning("NAs found in genetic data, imputing to the median (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
  }

  # PCA to determine number of latent factors
  # If K is not specified, it is calculated based on given K selection method
  if(is.null(K)){
    K <- select_K(gen, K_selection = K_selection, coords = coords,
               Kvals = Kvals, criticalpoint = criticalpoint, low = low,
               max.pc = max.pc, pca.select = pca.select,  perc.pca = perc.pca,
               choose.n.clust = choose.n.clust, criterion = criterion, max.n.clust = max.n.clust)
  }

  # Run LFMM
  results <- lfmm_run(gen, env, K = K, lfmm_method = lfmm_method, p_adj = p_adj, sig = sig, calibrate = calibrate)

  # Check qqplots
  print(lfmm_qqplot(results$df))

  # Make Manhattan plots
  print(lfmm_manhattanplot(results$df, sig))

  # Make table
  print(lfmm_table(results$df, top = TRUE, order = TRUE, nrow = 10))

  return(results)
}


#' Run LFMM
#'
#' @inheritParams lfmm_do_everything
#'
#' @export
#'
lfmm_run <- function(gen, env, K, lfmm_method = "ridge", p_adj = "fdr", sig = 0.05, calibrate = "gif"){
  # gen matrix
  genmat <- as.matrix(gen)
  # env matrix
  envmat <- as.matrix(env)

  # Remove NAs
  if(any(is.na(envmat))){
    warning("Missing values found in environmental data, removing rows with NAs")
    genmat <- genmat[complete.cases(envmat),]
    envmat <- envmat[complete.cases(envmat),]
  }

  # Run model
  if(lfmm_method == "ridge"){lfmm_mod <- lfmm::lfmm_ridge(genmat, envmat, K = K)}
  if(lfmm_method == "lasso"){lfmm_mod <- lfmm::lfmm_lasso(genmat, envmat, K = K)}

  # Perform association testing using the fitted model:
  lfmm_test_result <- lfmm::lfmm_test(Y = genmat,
                  X = envmat,
                  lfmm = lfmm_mod,
                  calibrate = calibrate)

  # If p_adj method is specified, perform p-value correction by column (by env variable)
  lfmm_test_result$adjusted.pvalue <- apply(dplyr::as_tibble(lfmm_test_result$calibrated.pvalue), 2, p.adjust, method = p_adj)

  # Stop if all p-values are NA
  if(all(is.na(lfmm_test_result$adjusted.pvalue))) stop("All p-values are NA")

  # Transfer column names
  colnames(lfmm_test_result$adjusted.pvalue) <- colnames(envmat)
  # Transfer rownames
  rownames(lfmm_test_result$adjusted.pvalue) <- colnames(genmat)

  # Make tidy dataframe of results
  result_df <- lfmm_df(lfmm_test_result)

  # Subset out candidate SNPs
  cand_snps <- result_df %>% dplyr::filter(adjusted.pvalue < 0.05)

  return(list(cand_snps = cand_snps, df = result_df, model = lfmm_mod, lfmm_test_result = lfmm_test_result, K = K))
}


#' Convert LFMM results into a tidy dataframe for downstream processing
#'
#' @param x df with LFMM test result (`lfmm_test_result` element)
#'
#' @return tidy dataframe with LFMM results with each SNP, its p-value, association with env var and other relevant statistics
#' @export
#'
#' @examples
lfmm_df <- function(x){
  # Extract names of elements from lfmm_test_result
  df_names <- names(x)[purrr::map_lgl(x, function(x) !is.null(rownames(x)))]

  # Apply helper function lfmm_test_tidy to extract stats for each SNP
  df <- purrr::map(df_names, lfmm_test_tidy, x) %>%
    purrr::reduce(dplyr::left_join, by = c("snp", "var"))

  return(df)
}

#' Helper function to tidy LFMM test results
#'
#' @param colname names of elements within lfmm_test_result
#' @param lfmm_test_result LFMM test results
#'
#' @return
#' @export
#'
#' @examples
lfmm_test_tidy <- function(colname, lfmm_test_result){
  x <- lfmm_test_result[[colname]]
  if(is.null(rownames(x))) rownames(x) <- paste0("snp", 1:nrow(x))
  df <- x %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(snp = rownames(x)) %>%
    tidyr::gather(var, stat, -snp)

  colnames(df) <- c("snp", "var", colname)

  return(df)
}

#' Create `gt` table of LFMM results
#'
#' @param df df element from \code{\link{lfmm_run}} results
#' @param sig alpha level for determining candidate snps (defaults to 0.5)
#' @param sig_only only include SNPs that exceeded the significance threshold in the table
#' @param top if there are SNPs that are significantly associated with multiple environmental variables, only display the top association (i.e., variable with the maximum B value; defaults to FALSE)
#' @param order if TRUE, will order rows by decreasing B value (defaults to FALSE and orders rows based on variable)
#' @param var display significant SNPs associated with particular environmental variable (defaults to NULL)
#' @param rows number of rows to include in table (defaults to NULL; will only include significant SNPs)
#' @param digits number of decimal points to include (defaults to 2)
#'
#' @return table of LFMM results
#' @export
#'
#' @examples
lfmm_table <- function(df, sig = 0.05, sig_only = TRUE, top = FALSE, order = FALSE, var = NULL, nrow = NULL, digits = 2, footnotes = TRUE){

  if(!is.null(var)) df <- df[df$var %in% var, ]
  if(sig_only) df <- df[df$adjusted.pvalue < sig, ]
  if(order) df <- df[order(abs(df$B), decreasing = TRUE),]
  if(top) df <- df %>%
      dplyr::group_by(snp) %>%
      dplyr::filter(abs(B) == max(abs(B)))
  if(!is.null(nrow)) {
    if(nrow > nrow(df)) nrow <- nrow(df)
    df <- df[1:nrow, ]
  }

  df <- df %>% dplyr::as_tibble() %>% dplyr::filter(dplyr::if_any(dplyr::everything(), ~ !is.na(.)))
  if(!is.null(digits)) df <- df %>% dplyr::mutate(dplyr::across(-c(var, snp), round, digits))

  if(nrow(df) == 0) {
    warning("No significant variants found, returning NULL object")
    return(NULL)
  }

  d <- max(abs(min(df$B, na.rm = TRUE)), abs(max(df$B, na.rm = TRUE)))

  # TODO: come back and remove adjusted p-value if no correction is applied
  colnames(df) <- c("snp", "variable", "B", "z-score", "p-value", "calibrated z-score", "calibrated p-value", "adjusted p-value")

  suppressWarnings(
    tbl <- df  %>%
      gt::gt() %>%
      gtExtras::gt_hulk_col_numeric("B", trim = TRUE, domain = c(-d,d))

  )

  if(footnotes) tbl <- tbl %>% gt::tab_footnote(footnote = "LFMM effect size", locations = gt::cells_column_labels(columns = B))

  tbl
}


#' K selection
#'
#' @param gen genotype matrix
#' @param K_selection method for performing K selection (can either by "tracy_widom" (default), "quick_elbow", or "tess")
#' @param coords if "tess" method is used, coordinates for TESS based K selection (defaults to NULL)
#' @param Kvals values of K to test if using "tess" method of K selection (defaults to 1:10)
#' @param criticalpoint if "tracy_widom" method is used, a numeric value corresponding to the significance level. If the significance level is 0.05, 0.01, 0.005, or 0.001, the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, respectively (defaults to 2.0234)
#' @param low if "quick_elbow" method is used, numeric, between zero and one, the threshold that defines whether a principal component explains 'much' of the variance (defaults to 0.08).
#' @param max.pc if "quick_elbow" method is used, maximum percentage of the variance to capture before the elbow (cumulative sum to PC 'n'; defaults to 0.90).
#' @param pca.select if "find_clusters" method is used, a character indicating the mode of selection of PCA axes, matching either "nbEig" or "percVar" (default). For "nbEig", the user has to specify the number of axes retained (interactively, or via n.pca). For "percVar", the user has to specify the minimum amount of the total variance to be preserved by the retained axes, expressed as a percentage (interactively, or via perc.pca).
#' @param perc.pca if "find_clusters" method is used, a numeric value between 0 and 100 indicating the minimal percentage of the total variance of the data to be expressed by the retained axes of PCA (defaults to 90).
#' @param choose.n.clust if "find_clusters" method is used, a logical indicating whether the number of clusters should be chosen by the user (defaults to FALSE), or automatically, based on a given criterion (argument criterion). It is HIGHLY RECOMMENDED to choose the number of clusters INTERACTIVELY, since i) the decrease of the summary statistics (BIC by default) is informative, and ii) no criteria for automatic selection is appropriate to all cases (see details in \code{find.cluster} documentation).
#' @param criterion if "find_clusters" method is used, a logical indicating whether the number of clusters should be chosen by the user (defaults to FALSE), or automatically, based on a given criterion (argument criterion). It is HIGHLY RECOMMENDED to choose the number of clusters INTERACTIVELY, since i) the decrease of the summary statistics (BIC by default) is informative, and ii) no criteria for automatic selection is appropriate to all cases (see details in \code{find.cluster} documentation).
#' @param max.n.clust if "find_clusters" method is used, an integer indicating the maximum number of clusters to be tried. Values of 'k' will be picked up between 1 and max.n.clust (defaults to 10)
#'
#' @return prints the best K value given the specified K selection procedure
#' @export
#'
#' @examples
select_K <- function(gen, K_selection = "tracy_widom", coords = NULL, Kvals = 1:10, criticalpoint = 2.023,
                  low = 0.08, max.pc = 0.9, pca.select = "percVar", perc.pca = 90, choose.n.clust = FALSE,
                  criterion = "diffNgroup", max.n.clust = 10){

  if(K_selection == "tracy_widom") K <- select_K_tw(gen, criticalpoint)

  if(K_selection == "quick_elbow") K <- select_K_elbow(gen, low, max.pc)

  if(K_selection == "tess") K <- select_K_tess(gen, coords, Kvals)

  if(K_selection == "find_clusters") K <- select_K_fc(gen,
                                                   pca.select = pca.select,
                                                   perc.pca = perc.pca,
                                                   choose.n.clust = choose.n.clust,
                                                   criterion = criterion,
                                                   max.n.clust = max.n.clust)

  return(K)
}


#' @describeIn select_K select K using Tracy-Widom Test
#' @param gen genotype matrix
#' @inheritParams AssocTests::tw
#'
#' @return
#' @export
#'
#' @examples
select_K_tw <- function(gen, criticalpoint = 2.0234){
  # Turn gen into df
  df <- data.frame(gen)

  # Run PCA
  pc <- prcomp(~., df, na.action = na.omit)

  # Get eig
  eig <- pc$sdev^2

  # Run Tracy-Widom test
  tw_result <- AssocTests::tw(eig, eigenL = length(eig), criticalpoint = criticalpoint)

  # Get K based on number of significant eigenvalues
  K <- tw_result$SigntEigenL

  # If K is zero, return 1
  if(K == 0) K <- 1

  return(K)
}

#' @describeIn select_K select K using PCA and \code{quick_elbow} method
#' @param gen genotype matrix
#' @inheritParams quick_elbow
#'
#' @return
#' @export
#'
#' @examples
select_K_elbow <- function(gen, low = 0.08, max.pc = 0.9){
  # Run PCA
  pc <- prcomp(gen)

  # Get eig
  eig <- pc$sdev^2

  # Estimate number of latent factors using quick_elbow (see general functions for description of how this function works)
  # This is a crude way to determine the number of latent factors that is based on an arbitrary "low" value
  K <- quick_elbow(eig, low = low, max.pc = max.pc)

  par(pty = "s",mfrow = c(1,1))
  plot(eig, xlab = 'PC', ylab = "Variance explained")
  abline(v = K, col = "red", lty = "dashed")

  return(K)
}


#' @describeIn select_K select K using TESS and \code{bestK} method
#' @param gen genotype matrix
#' @param coords coordinates for "tess"
#' @param Kvals values of K to test for "tess"
#' @param tess_method method to use for "tess"
#' @param ploidy ploidy for "tess"
#'
#' @return
#' @export
#'
#' @examples
select_K_tess <- function(gen, coords, Kvals = 1:10, tess_method = "projected.ls", ploidy = 2){
  # Run TESS for all K values
  tess3_obj <- tess3r::tess3(X = gen, coord = as.matrix(coords), K = Kvals, method = tess_method, ploidy = ploidy)

  # Plot x-validation results and indicate K-value that is automatically selected
  plot(tess3_obj, pch = 19, col = "blue",
       xlab = "Number of ancestral populations",
       ylab = "Cross-validation score")

  # Get best K value
  K <-  bestK(tess3_obj, Kvals)

  return(K)
}



#' @describeIn select_K select K using find.clusters method
#' @param gen a genotype matrix
#' @inheritParams adegenet::find.clusters
#' @return
#' @export
#'
#' @examples
select_K_fc <- function(gen, pca.select = "percVar", perc.pca = 90, choose.n.clust = FALSE,
                        criterion = "diffNgroup", max.n.clust = 10){

  fc <- adegenet::find.clusters(gen,
                                pca.select = pca.select,
                                perc.pca = perc.pca,
                                choose.n.clust = choose.n.clust,
                                criterion = criterion,
                                max.n.clust = max.n.clust)

  K <- max(as.numeric(fc$grp))
  return(K)
}


#' Quickly choose an elbow for a PC.
#' At variance below 5% per component, choose the largest % drop
#' Designed for variance percentages, but will also work given a full set of Evalues
#' Quickly estimate the 'elbow' of a scree plot (PCA)
#'
#' This function uses a rough algorithm to estimate a sensible 'elbow' to
#' choose for a PCA screeplot of eigenvalues. The function looks at an initial arbitrarily 'low'
#' level of variance and looks for the first eigenvalue lower than this. If the very first eigenvalue
#' is actually lower than this (i.e, when the PCs are not very explanatory) then this 'low' value is
#' iteratively halved until this is no longer the case. After starting below this arbitrary threshold
#' the drop in variance explained by each pair of consecutive PCs is standardized by dividing over the
#' larger of the pair. The largest percentage drop in the series below 'low' % is selected as the 'elbow'.
#'
#' @param varpc numeric, vector of eigenvalues, or 'percentage of variance' explained by datapoints for each principal component. If only using a partial set of components, should first pass to `estimate.eig.vpcs()` to estimate any missing eigenvalues
#' @param low numeric (between zero and one); the threshold that defines whether a principal component explains 'much' of the variance
#' @param max.pc maximum percentage of the variance to capture before the elbow (cumulative sum to PC 'n')
#' @return the number of principal components to keep, prior to the determined elbow cutoff
#' @export
#' @seealso \code{\link{estimate.eig.vpcs}}
#' @author Nicholas Cooper
#' @examples
#' # correlated data
#' mat <- sim.cor(100,50)
#' result <- princomp(mat)
#' eig <- result$sdev^2
#' elb.a <- quick_elbow(eig)
#' pca.scree.plot(eig,elbow=elb.a,M=mat)
#' elb.b <- quick_elbow(eig,low=.05) # decrease 'low' to select more components
#' pca.scree.plot(eig,elbow=elb.b,M=mat)
#' # random (largely independent) data, usually higher elbow #
#' mat2 <- generate.test.matrix(5,3)
#' result2 <- princomp(mat2)
#' eig2 <- result2$sdev^2
#' elb2 <- quick_elbow(result2$sdev^2)
#' pca.scree.plot(eig2,elbow=elb2,M=mat2)
quick_elbow <- function(varpc, low = 0.08, max.pc = 0.9) {
  # Ensure below sums to 1
  ee <- varpc/sum(varpc)

  # When no big components, then adjust 'low'
  while(low>=max(ee)) { low <- low/2 }
  lowie <- (ee<low) ; highie <- ee>low/8
  low.ones <- which(lowie & highie)
  others <- length(which(!lowie))
  if(length(low.ones)>0) {
    if(length(low.ones)==1) {
      elbow <- low.ones
    } else {
      set <- ee[low.ones]
      pc.drops <- abs(diff(set))/(set[1:(length(set)-1)])
      infz <- is.infinite(pc.drops)
      elbow <- which(pc.drops==max(pc.drops[!infz],na.rm=T))[1]+others
    }
  } else {
    # If somehow there are no small eigenvalues, just choose the elbow as the second last
    cat("No eigenvalues were significantly smaller than the previous\n")
    elbow <- length(ee)
  }
  if(tail(cumsum(ee[1:elbow]),1)>max.pc) {
    elbow <- which(cumsum(ee)>max.pc)[1]-1
  }
  if(elbow<1) {
    warning("Elbow calculation failed, return zero")
    return(0)
  }
  names(elbow) <- NULL
  return(elbow)
}

#' LFMM QQplot
#'
#' @param df dataframe of LFMM test results produced by \code{lfmm_df}
#'
#' @return a QQ-plot
#' @export
#'
#' @examples
lfmm_qqplot <- function(df){

  plt <- ggplot2::ggplot(df, ggplot2::aes(sample = -log10(adjusted.pvalue))) +
    ggplot2::stat_qq() +
    ggplot2::geom_abline(intercept = 0, slope = 1) +
    ggplot2::facet_wrap( ~ var, nrow = 1) +
    ggplot2::labs(x = NULL,
         y = "-log10(p)") +
    ggplot2::theme_bw() +
    ggplot2::theme(line = ggplot2::element_blank())

  return(plt)
}

#' LFMM Manhattan Plot
#'
#' @param df dataframe of lfmm test results produced by \code{lfmm_df}
#' @param sig significance cutoff
#'
#' @return
#' @export
#'
#' @examples
lfmm_manhattanplot <- function(df, sig, group = NULL, var = NULL){

  # subset variables
  if(!is.null(var)) df <- df[df$var %in% var, ]

  # Convert to df to not get tidy warnings about uninitialized columns
  df <- data.frame(df)
  df$type[df$adjusted.pvalue < sig] <- "Outlier"
  df$type[!(df$adjusted.pvalue < sig)] <- "Neutral"
  df$index <- 1:length(unique(df$snp))

  plt <-
    ggplot2::ggplot(df, ggplot2::aes(x = index, y = -log10(adjusted.pvalue))) +
    ggplot2::geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") +
    ggplot2::geom_point(alpha = 0.75, pch = 16, ggplot2::aes(col = type)) +
    ggplot2::scale_color_manual(values = c(rgb(0.7,0.7,0.7,0.5), "#F9A242FF"), na.translate = F) +
    ggplot2::xlab("snps") + ggplot2::ylab("-log10(p)") +
    ggplot2::geom_hline(yintercept=-log10(sig), linetype="dashed", color = "black", size=0.6) +
    ggplot2::guides(color=ggplot2::guide_legend(title="snp type")) +
    ggplot2::facet_wrap( ~ var, nrow = length(unique(df$var))) +
    ggplot2::xlab("position") + ggplot2::ylab("-log10(p)") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position="right",
                   legend.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   legend.box.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size=ggplot2::rel(.8)),
                   strip.text = ggplot2::element_text(size=11))
  return(plt)
}

