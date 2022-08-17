
#' ld_prune prunes SNPs based on linkage disequilibrium using `snpggdsLDpruning()` in SNPRelate package.
#' TODO: if vcf is path convert, otherwise no need to convert (does functionality exist within snpgdsVCF2GDS?)
#' TODO: should intermediate gds files be deleted? set to temp and overwrite?
#' TODO: print number of SNPs to be retained after LD-pruning
#'
#' @param vcf is the path to the vcf file containing all SNPs
#' @param nodes is the number of parallel processors (numeric)
#' @param ld.threshold is the threshold for LD pruning (numeric; 0 - 1; defaults to 0.6)
#' @param slide.max.n is the maximum number of SNPs in a sliding window (numeric; defaults to 100)
#' @param maf is the minor allele frequency cutoff (numeric; defaults to 0.05)
#' @param seed is the random starting seed (defaults to 1234)
#' @param method is the LD threshold method; default to corr which is r2 correlation coefficient
#' @param output_dir is the output file directory (default is here)
#' @param outname prefix name of output files (will append with param settings)
#' @param out_format output file format ("plink" will produce ped and map files while "vcf" will produce a vcf and a GDS)
#' Saves SNP GDS and ped (plink) files with retained SNPs

ld_prune <- function(vcf, outname, out_format, nodes = 1, ld.threshold = 0.6, slide.max.n = 100,
                    maf = 0.05, seed = 1234, method = c("corr")){

  # Specify output file name
  outfile_name <- paste(outname, "_LDpruned_r", ld.threshold, "_n", slide.max.n, sep="")

  # Write log file to track output
  sink(file = paste(outfile_name, "_LOGFILE.txt", sep=""))

  # Convert vcf to GDS file
  SNPRelate::snpgdsVCF2GDS(vcf, paste(outname, ".gds", sep=""))

  # Open GDS file
  genofile <- SNPRelate::snpgdsOpen(paste(outname, ".gds", sep=""), allow.duplicate = FALSE)

  set.seed(seed)

  # Define set of SNPs based on MAF, LD threshold, and window size
  snpset <- SNPRelate::snpgdsLDpruning(genofile,
                            ld.threshold = ld.threshold,
                            slide.max.n = slide.max.n,
                            maf = maf,
                            autosome.only = FALSE,
                            method = method)

  snp.id <- unlist(snpset)

  # Prune GDS based on the LD-pruned SNP set
  SNPRelate::snpgdsCreateGenoSet(src.fn = paste(outname, ".gds", sep=""), dest.fn = paste(outfile_name, ".gds", sep=""),
                      snp.id = snp.id)

  if (out_format == "plink") {
    # Open LD-pruned SNP set
    LDgenofile <- SNPRelate::snpgdsOpen(paste(outfile_name, ".gds", sep=""),
                                        allow.duplicate = FALSE)
    # Convert LD-pruned GDS to ped file
    SNPRelate::snpgdsGDS2PED(LDgenofile, ped.fn = paste(outfile_name, ".gds", sep=""),
                             sample.id = NULL, snp.id = NULL, use.snp.rsid = FALSE, verbose = TRUE)
  }

  if (out_format == "vcf") {
    # Convert from SNP GDS to GDS
    SeqArray::seqSNP2GDS(paste(outfile_name, ".gds", sep=""), paste(outfile_name, "seqarray.gds", sep=""))
    # Convert GDS to vcf
    SeqArray::seqGDS2VCF(paste(outfile_name, "seqarray.gds", sep=""), paste(outfile_name, ".vcf", sep=""), info.var = NULL, fmt.var = NULL, chr_prefix = "",
                         use_Rsamtools = TRUE, verbose = TRUE)
  }

  gdsfmt::showfile.gds(closeall = TRUE)

  sink()
}
