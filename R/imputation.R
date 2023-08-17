#' Impute NA values
#' NOTE: use extreme caution when using this form of simplistic imputation. We mainly provide this code for creating test datasets and highly discourage its use in analyses.
#' @param x matrix
#' @param f function to use for imputation (defaults to median)
#'
#' @return matrix of values with missing values imputed
#' @export
#' @family Imputation functions
simple_impute <- function(x, FUN = median) {
  x_noNA <- apply(x, 2, impute_helper, FUN)
  return(x_noNA)
}

#' Helper function for imputation
#' @export
#' @noRd
#' @family Imputation functions
impute_helper <- function(i, FUN = median) {
  i[which(is.na(i))] <- FUN(i[-which(is.na(i))], na.rm = TRUE)
  return(i)
}

#' Imputation of missing values using population structure inferred with `LEA::snmf`
#'
#' @param gen a dosage matrix, an object of class 'vcfR', a path to a .vcf file, or an object of type snmfProject
#' @param quiet whether to print results of cross-entropy scores (defaults to TRUE; only does so if more than one K-value); only displays run with minimum cross-entropy
#' @param save_output name prefix for saved .geno file, SNMF project file, and SNMF output file results (defaults to FALSE, in which no files are saved)
#' if TRUE, saves SNP GDS and ped (plink) files with retained SNPs in new directory; if FALSE returns object (defaults to TRUE)
#' gen <- load.snmfProject("example.snmfProject")
#'
#' @inheritParams LEA::snmf
#'
#' @return dosage matrix with imputed missing values
#' @export
#' @family Imputation functions
str_impute <- function(gen, K, entropy = TRUE, repetitions = 10, project = "new", quiet = TRUE, save_output = FALSE, output_filename = NULL) {
  if (is.null(output_filename)) filename <- "tmp" else filename <- output_filename

  if (inherits(gen, "snmfProject")) snmf_proj <- gen

  # Convert gen to .geno type file unless snmfProject provided
  if (!inherits(gen, "snmfProject")) {
    geno <- gen_to_geno(gen)
    # SNMF requires an input file saved to file (cannot accept an R object)
    LEA::write.geno(geno, here::here(paste0(filename, ".geno")))

    # Run SNMF
    snmf_proj <- LEA::snmf(input.file = here::here(paste0(filename, ".geno")), K = K, entropy = entropy, repetitions = repetitions, project = project)
  }

  # Look through directories
  bestK <- snmf_bestK(snmf_proj, K = K, quiet = quiet)

  if (!quiet) {
    print(plot_crossent(bestK$ce_values))
  }

  # Impute missing values based on snmf groupings
  LEA::impute(object = snmf_proj, input.file = here(paste0(filename, ".geno")), method = "mode", K = as.integer(bestK$K_value), run = as.integer(bestK$bestrun))

  # Read .lfmm file in
  imputed <- LEA::read.lfmm(here(paste0(filename, ".lfmm_imputed.lfmm")))

  # Add individual and variant names back in
  rownames(imputed) <- rownames(gen)
  colnames(imputed) <- colnames(gen)

  imputed <- geno_to_dosage(imputed)

  # Remove created files
  if (!save_output) {
    # Removes snmf project and associated snmf files
    LEA::remove.snmfProject(here::here(paste0(filename, ".snmfProject")))

    # Delete other associated files
    unlink(here::here(paste0(filename, ".geno")))
    unlink(here::here(paste0(filename, ".lfmm")))
    unlink(here::here(paste0(filename, ".lfmm_imputed.lfmm")))
  }
  return(imputed)
}

#' Helper function to select "best" K based on minimizing cross-entropy criteria from SNMF results
#'
#' @param snmf_proj object of type snmfProject
#' @param K integer corresponding to K-value
#' @param quiet whether to print results of cross-entropy scores (defaults to TRUE; only does so if more than one K-value); only displays run with minimum cross-entropy
#'
#' @return list with best K-value and run number and all cross-entropy scores
#' @export
#' @family Imputation functions
snmf_bestK <- function(snmf_proj, K, quiet) {
  if (length(K) == 1) {
    bestrun <- which.min(LEA::cross.entropy(snmf_proj, K = K))
    results <- list(K_value = K, run = bestrun)
  }

  if (length(K) > 1) {
    ce_values <- as.data.frame(purrr::map(K, snmf_crossent_helper, snmf_proj = snmf_proj, select_min = FALSE))
    ce_values <-
      ce_values %>%
      tibble::rownames_to_column(var = "run") %>%
      tidyr::pivot_longer(names_to = "K_value", values_to = "cross_entropy", -run)
    ce_values$run <- stringr::str_replace_all(ce_values$run, "run ", "")
    ce_values$K_value <- stringr::str_replace_all(ce_values$K_value, "K...", "")

    best <- ce_values %>% dplyr::slice(which.min(cross_entropy))
    results <- list(K_value = best$K_value, bestrun = best$run, ce_values = ce_values)

    if (!quiet) {
      print(plot_crossent(ce_values))
    }
  }
  return(results)
}

#' Helper function to retrieve cross entropy scores from SNMF project
#'
#' @param snmf_proj object of type snmfProject
#' @param K K-value(s)
#' @param select_min whether to return minimum
#'
#' @return cross entropy scores for given K
#' @export
#' @family Imputation functions
snmf_crossent_helper <- function(snmf_proj, K, select_min = TRUE) {
  if (select_min) results <- which.min(LEA::cross.entropy(snmf_proj, K = K))
  if (!select_min) results <- LEA::cross.entropy(snmf_proj, K = K)
  return(results)
}

#' Helper function to plot cross entropy scores from SNMF
#'
#' @param ce_values df with run, K-value, and cross entropy created in \link[algatr]{snmf_bestK}
#'
#' @return ggplots of cross entropy values compared to K-values (and across runs)
#' @export
plot_crossent <- function(ce_values) {
  if (length(unique(ce_values$run)) == 1) {
    plt <-
      ce_values %>%
      dplyr::group_by(K_value) %>%
      dplyr::slice(which.min(cross_entropy)) %>%
      ggplot2::ggplot(ggplot2::aes(x = K_value, y = cross_entropy)) +
      ggplot2::geom_point(size = 3, color = "red") +
      ggplot2::theme_bw() +
      ggplot2::ylab("Cross entropy value") +
      ggplot2::xlab("K value") +
      ggplot2::geom_vline(xintercept = best$K_value, linetype = "dashed", color = "blue")
    print(plt)
  }

  if (length(unique(ce_values$run)) > 1) {
    plt <-
      ce_values %>%
      ggplot2::ggplot(ggplot2::aes(x = run, y = cross_entropy)) +
      ggplot2::geom_point(size = 3, color = "red") +
      ggplot2::theme_bw() +
      ggplot2::ylab("Cross entropy value") +
      ggplot2::xlab("Run") +
      ggplot2::facet_grid(~K_value)
    print(plt)
  }
}

#' Convert dosage matrix or vcf to geno type object (N.B.: this only works for diploids!)
#'
#' @inheritParams str_impute
#'
#' @return matrix encoded as geno type object
#' @export
#' @family Imputation functions
gen_to_geno <- function(gen) {
  if (!is.matrix(gen)) if (!inherits(gen, "vcfR")) {
    gen <- vcfR::read.vcfR(gen)
    gen <- vcf_to_dosage(gen)
  }
  if (inherits(gen, "vcfR")) gen <- vcf_to_dosage(gen)

  # Recode data for geno type object
  gen[gen == 2] <- 555 # tmp placeholder so 0s aren't overwritten
  gen[gen == 0] <- 2 # two ref alleles (vcf 0/0 or dosage 0)
  gen[gen == 555] <- 0 # no ref alleles (vcf 1/1 or dosage 2)
  gen[is.na(gen)] <- 9 # missing data encoded as 9

  return(gen)
}

#' Convert lfmm/geno matrix to dosage matrix (N.B.: this only works for diploids!)
#'
#' @param geno matrix of LEA geno or lfmm format (i.e., 0 corresponds to zero reference alleles)
#'
#' @return matrix encoded as dosage type object (0 corresponds to two reference alleles)
#' @export
#' @family Imputation functions
geno_to_dosage <- function(geno) {
  # Recode data for geno type object
  geno[geno == 2] <- 555 # tmp placeholder so 2s aren't overwritten
  geno[geno == 0] <- 2 # two ref alleles (vcf 0/0 or dosage 0)
  geno[geno == 555] <- 0 # no ref alleles (vcf 1/1 or dosage 2)
  geno[geno == 9] <- NA # missing data encoded as NA

  return(geno)
}
