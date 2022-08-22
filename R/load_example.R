#' Load example data
#'
#' @param quiet whether to load quietly or not
#'
#' @return
#' @export
#'
#' @examples
load_example <- function(quiet = FALSE){
  # load all data
  utils::data(list = c("liz_vcf", "liz_coords", "CA_env", "gendist"))

  if (!quiet) {
    # give message with information about objects
    return(message(cat(
      crayon::cyan(crayon::bold("\n---------------- example dataset ----------------\n")),
      crayon::silver("\nObjects loaded:"),
      crayon::yellow(crayon::bold("\n*liz_vcf*")),
      crayon::yellow(paste0("vcfR object (1000 loci x 53 samples)")),
      crayon::green(crayon::bold("\n*liz_coords*")), crayon::green("dataframe with x and y coordinates"),
      crayon::blue(crayon::bold("\n*CA_env*")), crayon::blue("RasterStack with PC environmental layers"),
      crayon::magenta(crayon::bold("\n*gendist*")), crayon::magenta("Genetic distance matrix"),
      crayon::cyan(crayon::bold("\n\n-------------------------------------------------\n"))
    )))
  }
}
