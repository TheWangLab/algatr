#' Calculate genetic distances
#'
#' @param vcf path to vcf file or a `vcfR` type object
#' @param dist_type the type of genetic distance to calculate (options: "euclidean" (default), "bray_curtis", "dps" for proportion of shared alleles, "plink", or "pc" for PC-based)
#' @param plink_file if "plink" dist_type is used, path to plink distance file (typically ".dist"; required only for calculating plink distance)
#' @param plink_id_file if "plink" dist_type is used, path to plink id file (typically ".dist.id"; required only for calculating plink distance)
#' @param npc_selection if "pc" dist_type is used, how to perform K selection (options: "auto" for automatic selection based on significant eigenvalues from Tracy-Widom test (default), or "manual" to examine PC screeplot and enter no. PCs into console)
#' @param criticalpoint if "pc" dist_type is used with "auto" npc_selection, the critical point for the significance threshold for the Tracy-Widom test within the PCA (defaults to 2.0234 which corresponds to an alpha of 0.01)
#'
#' @details
#' Euclidean and Bray-Curtis distances calculated using the ecodist package: Goslee, S.C. and Urban, D.L. 2007. The ecodist package for dissimilarity-based analysis of ecological data. Journal of Statistical Software 22(7):1-19. DOI:10.18637/jss.v022.i07.
#' Proportions of shared alleles calculated using the adegenet package: Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics. doi:10.1093/bioinformatics/btr521.
#' For calculating proportions of shared alleles, missing values are ignored (i.e., prop shared alleles calculated from present values; no scaling performed)
#'
#' @return pairwise distance matrix for given distance metric
#' @export
gen_dist <- function(vcf = NULL, dist_type = "euclidean", plink_file = NULL, plink_id_file = NULL, npc_selection = "auto", criticalpoint = 2.0234) {
  # Import vcf if provided --------------------------------------------------

  if (!is.null(vcf)) if (!inherits(vcf, "vcfR")) vcf <- vcfR::read.vcfR(vcf)

  # Calculate Euclidean distances -------------------------------------------

  if (dist_type == "euclidean") {
    # Convert to genlight and matrix
    gl <- vcfR::vcfR2genlight(vcf)
    mat <- as.matrix(gl)

    # Perform imputation with warning
    if (any(is.na(mat))) {
      mat <- simple_impute(mat, median)
      warning("NAs found in genetic data, imputing to the median (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
    }

    dists <- ecodist::distance(mat, method = "euclidean")
    dists <- as.matrix(dists)
    return(as.data.frame(dists))
  }

  # Calculate Bray-Curtis distances -----------------------------------------

  if (dist_type == "bray_curtis") {
    # Convert to genlight and matrix
    gl <- vcfR::vcfR2genlight(vcf)
    mat <- as.matrix(gl)

    # Perform imputation with warning
    if (any(is.na(mat))) {
      mat <- simple_impute(mat, median)
      warning("NAs found in genetic data, imputing to the median (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
    }

    # Check for NAs
    if (any(is.na(mat))) {
      stop("NA values found in genetic data")
    }

    dists <- ecodist::distance(mat, method = "bray-curtis")
    dists <- as.matrix(dists)
    return(as.data.frame(dists))
  }

  # Calculate proportion of shared alleles ----------------------------------

  if (dist_type == "dps") {
    # Convert to genind
    genind <- vcfR::vcfR2genind(vcf)
    dists <- adegenet::propShared(genind)
    return(as.data.frame(dists))
  }

  # Process Plink distance output files -------------------------------------

  if (dist_type == "plink") {
    dists <- as.data.frame(readr::read_tsv(plink_file, col_names = FALSE))
    plink_names <- readr::read_tsv(plink_id_file, col_names = FALSE) %>%
      dplyr::select(-`X1`) %>%
      as.matrix()
    # Assign row and col names according to sampleID
    rownames(dists) <- plink_names
    colnames(dists) <- plink_names
    return(dists)
  }

  # PC-based dist -----------------------------------------------------------

  if (dist_type == "pc") {
    # Convert to genlight
    gl <- vcfR::vcfR2genlight(vcf)
    mat <- as.matrix(gl) # to check for NAs

    # Perform imputation with warning
    if (any(is.na(mat))) {
      mat <- simple_impute(mat, median)
      gl <- adegenet::as.genlight(mat)
      warning("NAs found in genetic data, imputing to median (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
    }

    # Check for NAs
    if (any(is.na(mat))) {
      stop("NA values found in genetic data")
    }

    # Perform PCA
    pc <- stats::prcomp(gl)

    # Get eig
    eig <- pc$sdev^2

    # Automatic npc selection based on number of significant eigenvalues
    if (npc_selection == "auto") {
      # Run Tracy-Widom test
      # NOTE: critical point corresponds to significance level.
      # If the significance level is 0.05, 0.01, 0.005, or 0.001,
      # the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, respectively.
      # The default is 2.0234.
      tw_result <- AssocTests::tw(eig, eigenL = length(eig), criticalpoint = criticalpoint)
      npc <- tw_result$SigntEigenL
    }

    # Manual npc selection: screeplot printout and selecting no. PCs to retain
    if (npc_selection == "manual") {
      stats::screeplot(pc, type = "barplot", npcs = 10, main = "PCA Eigenvalues")
      npc <- as.numeric(readline(prompt = "Number of PC axes to retain:"))
    }

    # Calculate PC-based distance
    dists <- as.matrix(dist(pc$x[, 1:npc], diag = TRUE, upper = TRUE))

    return(as.data.frame(dists))
  }
}

#' Plot the relationship between two distance metrics
#'
#' @param dist_x df containing square distance matrix for x axis
#' @param dist_y df containing square distance matrix for y axis
#' @param metric_name_x name of distance metric for x axis; if DPS used, must be `"dps"`
#' @param metric_name_y name of distance metric for y axis; if DPS used, must be `"dps"`
#'
#' @export
gen_dist_corr <- function(dist_x, dist_y, metric_name_x, metric_name_y) {
  # Check to ensure sample IDs match ----------------------------------------

  if (all(rownames(dist_x) == rownames(dist_y)) == FALSE) {
    stop("Sample IDs do not match")
  }

  # Melt data from square to long -------------------------------------------

  # Assign NAs to upper triangle of square matrix
  dist_x[upper.tri(dist_x, diag = FALSE)] <- NA

  melt_x <- dist_x %>%
    tibble::rownames_to_column(var = "comparison") %>%
    tidyr::pivot_longer(cols = -(comparison)) %>%
    na.omit() %>%
    dplyr::filter(comparison != name) %>%
    dplyr::rename(!!metric_name_x := value)

  melt_y <- dist_y %>%
    tibble::rownames_to_column(var = "comparison") %>%
    tidyr::pivot_longer(cols = -(comparison)) %>%
    na.omit() %>%
    dplyr::filter(comparison != name) %>%
    dplyr::rename(!!metric_name_y := value)


  # Build plots -------------------------------------------------------------

  if (metric_name_x == "dps" || metric_name_y == "dps") {
    joined <- dplyr::full_join(melt_x, melt_y) %>%
      dplyr::mutate(rev_dps = (1 - dps))
    joined %>%
      ggplot2::ggplot(ggplot2::aes_string(x = metric_name_x, y = metric_name_y)) +
      ggplot2::geom_abline(ggplot2::aes(intercept = 0.0, slope = 1), color = "gray") +
      ggplot2::geom_point(color = "black", size = .2, alpha = .5)
  } else {
    joined <- dplyr::full_join(melt_x, melt_y)
    joined %>%
      ggplot2::ggplot(ggplot2::aes_string(x = metric_name_x, y = metric_name_y)) +
      ggplot2::geom_abline(ggplot2::aes(intercept = 0.0, slope = 1), color = "gray") +
      ggplot2::geom_point(color = "black", size = .2, alpha = .5)
  }
}

#' Make heatmap of genetic distances
#'
#' @param dist Matrix of genetic distances
#'
#' @return heatmap of genetic distances
#' @export
gen_dist_hm <- function(dist) {
  if (!is.null(dist)) if (!inherits(dist, "data.frame")) dist <- as.data.frame(dist)

  dist %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::gather("sample_comp", "dist", -"sample") %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = sample_comp, fill = dist)) +
    ggplot2::geom_tile() +
    ggplot2::coord_equal() +
    viridis::scale_fill_viridis(option = "inferno") +
    ggplot2::xlab("Sample") +
    ggplot2::ylab("Sample") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))
}
