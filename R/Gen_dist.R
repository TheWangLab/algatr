library("adegenet")

#' Function to measure genetic distances between samples, used as input for various downstream analyses
#' Genetic distances are measured by proportions of shared alleles (DPS)
#' 
#' @param genind genind object containing genotypes (can contain missing values)
#' @param species specifies species name (Genus_species) to save output file
#' 
#' @return square matrix with pairwise proportions of shared alleles (0 to 1), saved as csv

genDist <- function(genind, species){
  
  # calculate shared alleles from genind object
  DPSdists <- propShared(genind)
  
  # define path to save output files to
  path = "../data/GenDist"
  
  write_csv(as.data.frame(DPSdists), (paste(path, species, "_DPS_data.csv", sep = "")))
}