#' LDprune prunes SNPs based on linkage disequilibrium using snpggdsLDpruning() in SNPRelate package.
#'
#' @param vcf is the vcf file containing all SNPs
#' @param species is the species name; used for file naming
#' @param nodes is the number of parallel processors (numeric)
#' @param ld.threshold is the threshold for LD pruning (numeric; 0 - 1)
#' @param slide.max.n is the maximum number of SNPs in a sliding window (numeric)
#' @param maf is the minor allele frequency cutoff (numeric)
#' @param seed is the random starting seed
#' @param method is the LD threshold method; corr is r2 correlation coefficient
#' @param output_dir is the output file directory
#'
#' Saves SNP GDS and ped (plink) files with retained SNPs as "Species_LDpruned_params"

LDprune <- function(vcf, species, output_dir, nodes = 1, ld.threshold = 0.6, slide.max.n = 1000,
                    maf = 0.05, seed = 1234, method = c("corr")){

  # make outfile name containing params and species
  outfile_name <- paste(output_dir, species, "_LDpruned_r", ld.threshold, "_n", slide.max.n, sep="")

  # write log file to track output
  sink(file = paste(outfile_name, "_LOGFILE.txt", sep=""))

  # convert vcf to gds file
  snpgdsVCF2GDS(vcf, paste(output_dir, species, ".gds", sep=""))

  # open gds file
  genofile <- snpgdsOpen(paste(output_dir, species, ".gds", sep=""), allow.duplicate = FALSE)

  set.seed(seed)

  # define set of snps based on maf, ld threshold, and window size
  snpset <- snpgdsLDpruning(genofile,
                            ld.threshold = ld.threshold,
                            slide.max.n = slide.max.n,
                            maf = maf,
                            autosome.only = FALSE,
                            method = method)

  snp.id <- unlist(snpset)

  # prune gds based on the LD-pruned snp set
  snpgdsCreateGenoSet(src.fn = paste(output_dir, species, ".gds", sep=""), dest.fn = paste(outfile_name, ".gds", sep=""),
                      snp.id = snp.id)

  # save LD-pruned snp set as gds
  LDgenofile <- snpgdsOpen(paste(outfile_name, ".gds", sep=""),
                           allow.duplicate = FALSE)

  # convert LD-pruned gds to ped file
  snpgdsGDS2PED(LDgenofile, ped.fn = paste(outfile_name, ".gds", sep=""),
                sample.id = NULL, snp.id = NULL, use.snp.rsid = FALSE, verbose = TRUE)

  showfile.gds(closeall = TRUE)

  sink()
}

