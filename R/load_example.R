#' Load example data
#'
#' Example dataset from [Bouzid et al. 2022](https://doi.org/10.1111/mec.15836). The code used to create this dataset can be found under data-raw/DATASET
#'
#' @param quiet whether to load without messages
#'
#' @return
#' @export
#'
#' @examples
load_example <- function(quiet = FALSE){
  # load all data
  utils::data(list = c("liz_vcf", "liz_coords", "liz_gendist", "CA_env"))

  if (!quiet) {
    # Give message with information about objects
    return(message(cat(
      crayon::cyan(crayon::bold("\n---------------- example dataset ----------------\n")),
      crayon::blue("\nObjects loaded:"),
      crayon::yellow(crayon::bold("\n*liz_vcf*")),
      crayon::yellow(paste0("vcfR object (1000 loci x 53 samples)")),
      crayon::yellow(crayon::bold("\n*liz_gendist*")),
      crayon::yellow(paste0("genetic distance matrix (Plink Distance)")),
      crayon::green(crayon::bold("\n*liz_coords*")), crayon::green("dataframe with x and y coordinates"),
      crayon::magenta(crayon::bold("\n*CA_env*")), crayon::magenta("RasterStack with PC environmental layers"),
      crayon::cyan(crayon::bold("\n\n-------------------------------------------------\n"))
    )))
  }
}
