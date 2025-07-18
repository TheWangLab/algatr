#' Calculate genetic distances
#'
#' @param gen path to vcf file, a `vcfR` type object, or a dosage matrix
#' @param dist_type the type of genetic distance to calculate (options: `"euclidean"` (default), `"bray_curtis"`, `"dps"` for proportion of shared alleles (requires vcf), `"plink"`, or `"pc"` for PC-based)
#' @param plink_file if `"plink"` dist_type is used, path to plink distance file (typically ".dist"; required only for calculating plink distance). File must be a **square** distance matrix.
#' @param plink_id_file if `"plink"` dist_type is used, path to plink id file (typically ".dist.id"; required only for calculating plink distance)
#' @param npc_selection if `dist_type = "pc"`, how to perform K selection (options: `"auto"` for automatic selection based on significant eigenvalues from Tracy-Widom test (default), or `"manual"` to examine PC screeplot and enter no. PCs into console)
#' @param criticalpoint if `dist_type = "pc"` used with `npc_selection = "auto"`, the critical point for the significance threshold for the Tracy-Widom test within the PCA (defaults to 2.0234 which corresponds to an alpha of 0.01)
#'
#' @details
#' Euclidean and Bray-Curtis distances calculated using the ecodist package: Goslee, S.C. and Urban, D.L. 2007. The ecodist package for dissimilarity-based analysis of ecological data. Journal of Statistical Software 22(7):1-19. DOI:10.18637/jss.v022.i07.
#' Proportions of shared alleles calculated using the adegenet package: Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics. doi:10.1093/bioinformatics/btr521.
#' For calculating proportions of shared alleles, missing values are ignored (i.e., prop shared alleles calculated from present values; no scaling performed)
#'
#' @return pairwise distance matrix for given distance metric
#' @export
gen_dist <- function(gen = NULL, dist_type = "euclidean", plink_file = NULL, plink_id_file = NULL,
                     npc_selection = "auto", criticalpoint = 2.0234) {

  # Process input data ------------------------------------------------------
  # Read in vcf if path provided
  if (is.character(gen)) gen <- vcfR::read.vcfR(gen)
  # Convert vcf to dosage matrix
  if (inherits(gen, "vcfR") & (dist_type == "euclidean" | dist_type == "bray_curtis" | dist_type == "pc")) {
    gen <- vcf_to_dosage(gen)
    # Perform imputation with warning
    if (any(is.na(gen))) {
      gen <- simple_impute(gen, median)
      warning("NAs found in genetic data, imputing to the median (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
    }
  }

  # Calculate Euclidean distances -------------------------------------------
  if (dist_type == "euclidean") {
    # Check for NAs
    if (any(is.na(gen))) {
      stop("NA values found in genetic data")
    }
    dists <- ecodist::distance(gen, method = "euclidean")
    dists <- as.matrix(dists)
    return(as.data.frame(dists))
  }

  # Calculate Bray-Curtis distances -----------------------------------------
  if (dist_type == "bray_curtis") {
    # Check for NAs
    if (any(is.na(gen))) {
      stop("NA values found in genetic data")
    }

    dists <- ecodist::distance(gen, method = "bray-curtis")
    dists <- as.matrix(dists)
    return(as.data.frame(dists))
  }

  # Calculate proportion of shared alleles ----------------------------------
  if (dist_type == "dps") {
    if (!inherits(gen, "vcfR")) stop("VCF file required for calculating DPS distances")
    # Convert to genind
    genind <- vcfR::vcfR2genind(gen)
    # Show DPS warning about previously incorrect calculation
    dps_warning()
    dists <- 1 - adegenet::propShared(genind)
    return(as.data.frame(dists))
  }

  # Process Plink distance output files -------------------------------------
  if (dist_type == "plink") {
    if (is.null(plink_file)) stop("No plink distance file provided")
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
    # Check for NAs
    if (any(is.na(gen))) {
      stop("NA values found in genetic data")
    }

    gl <- adegenet::as.genlight(gen)

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
#' @return scatterplot comparing two user-defined genetic distance metrics
#' @export
gen_dist_corr <- function(dist_x, dist_y, metric_name_x, metric_name_y) {
  if (!is.null(dist_x)) if (!inherits(dist_x, "data.frame")) dist_x <- as.data.frame(dist_x)
  if (!is.null(dist_y)) if (!inherits(dist_y, "data.frame")) dist_y <- as.data.frame(dist_y)

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
