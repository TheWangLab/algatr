
#' Calculate genetic distances
#'
#' @param vcf path to vcf file or a `vcfR` type object
#' @param plink_file path to plink distance file (typically ".dist"; required only for calculating plink distance)
#' @param plink_id_file path to plink id file (typically ".dist.id"; required only for calculating plink distance)
#' @param dist_type the type of genetic distance to calculate (TODO[EAC]: Switch the order so that dist_type comes before the plink files since dist_type is always required but plink files are not)
#' @param criticalpoint the critical point for the significance threshold for the Tracy Widom test within the PCA used to determine the number of PCs automatically (TODO[EAC]: add argument for automatic vs manual selection of PCs (e.g. have the function print a scree plot and then users can select a number of PCs and also have an argument for this function that is nPC ( i think I have this added in TESS)))
#'
#' @details
#' Euclidean and Bray-Curtis distances calculated using the ecodist package: Goslee, S.C. and Urban, D.L. 2007. The ecodist package for dissimilarity-based analysis of ecological data. Journal of Statistical Software 22(7):1-19. DOI:10.18637/jss.v022.i07.
#' Proportions of shared alleles calculated using the adegenet package: Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics. doi:10.1093/bioinformatics/btr521.
#' For PC-based distances, mean-based imputation calculated using the dartR package: Gruber, B, Unmack, PJ, Berry, OF, Georges, A. dartr: An r package to facilitate analysis of SNP data generated from reduced representation genome sequencing. Mol Ecol Resour. 2018; 18: 691-699. https://doi.org/10.1111/1755-0998.12745

#' @return pairwise distance matrix for given distance metric
#'
gen_dist <- function(vcf = NULL, plink_file = NULL, plink_id_file = NULL, dist_type, criticalpoint = 2.0234){

  # Import vcf if provided --------------------------------------------------

  if(!is.null(vcf)) if(!inherits(vcf, "vcfR")) vcf <- vcfR::read.vcfR(vcf)

  # Calculate Euclidean distances -------------------------------------------

  if (dist_type == "euclidean") {
    # Convert to genlight and matrix
    gl <- vcfR::vcfR2genlight(vcf)
    mat <- as.matrix(gl)

    # Perform imputation with warning
    if(any(is.na(mat))){
      mat <- simple_impute(mat, median)
      warning("NAs found in genetic data, imputing to the median (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
    }

    # Check for NAs
    # TODO[EAC]: this code is repetitive with the above expression
    if(any(is.na(mat))){
      stop("NA values found in genetic data")
    }

    dists <- ecodist::distance(mat, method = "euclidean")
    dists <- as.matrix(dists)
    return(as.data.frame(dists))
  }

  # Calculate Bray-Curtis distances -----------------------------------------

  if (dist_type == "bray-curtis") {
    # Convert to genlight and matrix
    gl <- vcfR::vcfR2genlight(vcf)
    mat <- as.matrix(gl)

    # Perform imputation with warning
    if(any(is.na(mat))){
      mat <- simple_impute(mat, median)
      warning("NAs found in genetic data, imputing to the median (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
    }

    # Check for NAs
    if(any(is.na(mat))){
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
    # TODO[EAC]: include in function description how adegent deals with NAs?
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

  if (dist_type == "pc"){
    # Convert to genlight
    gl <- vcfR::vcfR2genlight(vcf)
    mat <- as.matrix(gl) # to check for NAs

    # Perform imputation with warning
    if(any(is.na(mat))){
      length <- rep(1, length(gl$ind.names))
      adegenet::strata(gl) <- as.data.frame(length)
      adegenet::setPop(gl) <- ~length
      gl <- dartR::gl.impute(gl, method = "frequency")
      mat <- as.matrix(gl)
      # TODO[EAC]: why is imputation to the mean performed for PC and then imputation to the median performed for other measures? Probably should make consistent or give an explanation or make the choice of summarizing function an argument
      warning("NAs found in genetic data, imputing to mean (NOTE: this simplified imputation approach is strongly discouraged. Consider using another method of removing missing data)")
    }

    # Check for NAs
    if(any(is.na(mat))){
      stop("NA values found in genetic data")
    }

    # Perform PCA
    pc <- stats::prcomp(gl)

    # Get eig
    eig <- pc$sdev^2

    # Run Tracy-Widom test
    # NOTE: critical point corresponds to significance level.
    # If the significance level is 0.05, 0.01, 0.005, or 0.001,
    # the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, respectively.
    # The default is 2.0234.
    tw_result <- AssocTests::tw(eig, eigenL = length(eig), criticalpoint = criticalpoint)

    # TODO: see note above about adding automatic vs manual seleciton options as in TESS
    # Get K based on number of significant eigenvalues
    K <- tw_result$SigntEigenL

    # Calculate PC distance based on significant PCs
    dists <- as.matrix(dist(pc$x[,1:K], diag = TRUE, upper = TRUE))
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
gen_dist_corr <- function(dist_x, dist_y, metric_name_x, metric_name_y){

  # Check to ensure sample IDs match ----------------------------------------

  if (all(rownames(dist_x) == rownames(dist_y)) == FALSE) {
    stop("Sample IDs do not match")
  }

  # Melt data from square to long -------------------------------------------

  melt_x = harrietr::melt_dist(as.matrix(dist_x)) %>%
    dplyr::rename(!!metric_name_x := dist)
  melt_y = harrietr::melt_dist(as.matrix(dist_y)) %>%
    dplyr::rename(!!metric_name_y := dist)


  # Build plots -------------------------------------------------------------

  if (metric_name_x == "dps" || metric_name_y == "dps") {
    joined <- dplyr::full_join(melt_x, melt_y) %>%
      dplyr::mutate(rev_dps = (1-dps))
    joined %>%
      ggplot2::ggplot(ggplot2::aes_string(x=metric_name_x, y=metric_name_y)) +
      ggplot2::geom_abline(ggplot2::aes(intercept=0.0, slope=1), color="gray") +
      ggplot2::geom_point(color="black", size=.2, alpha = .5)
  } else {
    joined <- dplyr::full_join(melt_x, melt_y)
    joined %>%
      ggplot2::ggplot(ggplot2::aes_string(x=metric_name_x, y=metric_name_y)) +
      ggplot2::geom_abline(ggplot2::aes(intercept=0.0, slope=1), color="gray") +
      ggplot2::geom_point(color="black", size=.2, alpha = .5)
  }
}

#' Make heatmap of genetic distances
#'
#' @param dist Matrix of genetic distances
#'
#' @return
#' @export
#'
#' @examples
gen_dist_hm <- function(dist){

  if(!is.null(dist)) if(!inherits(dist, "data.frame")) dist <- as.data.frame(dist)

  dist %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::gather("sample_comp", "dist", -"sample") %>%
    ggplot2::ggplot(ggplot2::aes(x = sample, y = sample_comp, fill = dist)) +
    ggplot2::geom_tile() +
    ggplot2::coord_equal() +
    viridis::scale_fill_viridis(option = "inferno") +
    xlab("Sample") +
    ylab("Sample") +
    theme(axis.text.x = element_text(angle = 90))
}
