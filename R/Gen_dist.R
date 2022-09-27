
#' Calculate genetic distances
#'
#' @param dist_type is the type of genetic distance to calculate
#' @param criticalpoint is the critical point for the significance threshold for the TW test within the PCA
#' @param vcf_file path to vcf file
#' @param plink_file path to plink distance file (typically ".dist")
#' @param plink_id_file path to plink id file (typically ".dist.id")
#'
#' @details
#' Euclidean and Bray-Curtis distances calculated using the ecodist package: Goslee, S.C. and Urban, D.L. 2007. The ecodist package for dissimilarity-based analysis of ecological data. Journal of Statistical Software 22(7):1-19. DOI:10.18637/jss.v022.i07.
#' Proportions of shared alleles calculated using the adegenet package: Jombart T. and Ahmed I. (2011) adegenet 1.3-1: new tools for the analysis of genome-wide SNP data. Bioinformatics. doi:10.1093/bioinformatics/btr521.
#' For PC-based distances, mean-based imputation calculated using the dartR package: Gruber, B, Unmack, PJ, Berry, OF, Georges, A. dartr: An r package to facilitate analysis of SNP data generated from reduced representation genome sequencing. Mol Ecol Resour. 2018; 18: 691-699. https://doi.org/10.1111/1755-0998.12745

#'
#' @return pairwise distance matrix for given distance metric
#'
gen_dist_calc <- function(vcf_file, plink_file, plink_id_file, dist_type, criticalpoint = 2.0234){

  # Calculate Euclidean distances -------------------------------------------

  if (dist_type == "euclidean") {
    # Read in vcf file
    vcf <- vcfR::read.vcfR(vcf_file)
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

    dists <- ecodist::distance(mat, method = "euclidean")
    dists <- as.matrix(dists)
    return(as.data.frame(dists))
  }

  # Calculate Bray-Curtis distances -----------------------------------------

  if (dist_type == "bray-curtis") {
    # Read in vcf file
    vcf <- vcfR::read.vcfR(vcf_file)
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
    # Read in vcf file
    vcf <- vcfR::read.vcfR(vcf_file)
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

  if (dist_type == "pc"){
    # Read in vcf file
    vcf <- vcfR::read.vcfR(vcf_file)
    # Convert to genlight
    gl <- vcfR::vcfR2genlight(vcf)
    mat <- as.matrix(gl) # to check for NAs

    # Perform imputation with warning
    if(any(is.na(mat))){
      length <- rep(1, length(gl$ind.names))
      strata(gl) <- as.data.frame(length)
      setPop(gl) <- ~length
      gl <- dartR::gl.impute(gl, method = "frequency")
      mat <- as.matrix(gl)
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
#' @param metric_name_x name of distance metric for x axis
#' @param metric_name_y name of distance metric for y axis
#'
gen_dist_corr <- function(dist_x, dist_y, metric_name_x, metric_name_y){
  # Melt from square to long
  melt_x = harrietr::melt_dist(as.matrix(dist_x)) %>%
    dplyr::rename(!!metric_name_x := dist)
  melt_y = harrietr::melt_dist(as.matrix(dist_y)) %>%
    dplyr::rename(!!metric_name_y := dist)

  if (metric_name_x == "dps" || metric_name_y == "dps"){
    # TODO [EAC]: should check whether inds are the same across datasets, return something if not
    joined <- dplyr::full_join(melt_x, melt_y) %>%
      dplyr::mutate(rev_dps = (1-dps))
    joined %>%
      ggplot(aes_string(x=metric_name_x, y=metric_name_y)) +
      geom_abline(aes(intercept=0.0, slope=1), color="gray") +
      geom_point(color="black", size=.2, alpha = .5)
  } else {
    # TODO [EAC]: should check whether inds are the same across datasets, return something if not
    joined <- dplyr::full_join(melt_x, melt_y)
    joined %>%
      ggplot(aes_string(x=metric_name_x, y=metric_name_y)) +
      geom_abline(aes(intercept=0.0, slope=1), color="gray") +
      geom_point(color="black", size=.2, alpha = .5)
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
  as.data.frame(env_dist$CA_rPCA1) %>%
    tibble::rownames_to_column("sample") %>%
    tidyr::gather("sample_comp", "dist", -"sample") %>%
    ggplot2::ggplot(ggplot2::aes(x = as.numeric(sample), y = as.numeric(sample_comp), fill = dist)) +
    ggplot2::geom_tile() +
    ggplot2::coord_equal() +
    viridis::scale_fill_viridis() +
    xlab("Sample") +
    ylab("Sample")
}
