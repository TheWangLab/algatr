#' Calculate genetic distances
#'
#' @param vcf path to vcf file or a `vcfR` type object
#' @param plink_file path to plink distance file (typically ".dist")
#' @param plink_id_file path to plink id file (typically ".dist.id")
#' @param dist_type is the type of genetic distance to calculate
#' @param criticalpoint is the critical point for the significance threshold for the TW test within the PCA
#'
#' @return pairwise distance matrix for given distance metric
#'
gen_dist <- function(vcf = NULL, plink_file = NULL, plink_id_file = NULL, dist_type, criticalpoint = 2.0234){


  # Import vcf data if provided  --------------------------------------------

  if(!is.null(vcf)) if(!inherits(vcf, "vcfR")) vcf <- vcfR::read.vcfR(vcf)

  # Calculate Euclidean distances -------------------------------------------

  if (dist_type == "euclidean") {
    # Convert to genlight and matrix
    gl <- vcfR::vcfR2genlight(vcf)
    mat <- as.matrix(gl)
    dists <- ecodist::distance(mat, method="euclidean")
    dists <- as.matrix(dists)
    return(as.data.frame(dists))
  }

  # Calculate Bray-Curtis distances -----------------------------------------

  if (dist_type == "bray-curtis") {
    # Convert to genlight and matrix
    gl <- vcfR::vcfR2genlight(vcf)
    mat <- as.matrix(gl)
    dists <- ecodist::distance(mat, method="bray-curtis")
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

  if (dist_type == "pc"){
    # Convert to genlight
    gl <- vcfR::vcfR2genlight(vcf)
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
