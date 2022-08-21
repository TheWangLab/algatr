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
  utils::data(list = c("liz_vcf", "liz_coords", "CA_env"))

  # temp fix for rda issue with env data
  CA_env <- raster::stack(system.file("extdata", "CA_env.tif", package = "algatr"))
  assign("CA_env", CA_env, envir = .GlobalEnv)

  if (!quiet) {
    # give message with information about objects
    return(message(cat(
      crayon::cyan(crayon::bold("\n---------------- example dataset ----------------\n")),
      crayon::silver("\nObjects loaded:"),
      crayon::yellow(crayon::bold("\n*liz_vcf*")),
      crayon::yellow(paste0("vcfR object (1000 loci x 53 samples)")),
      crayon::green(crayon::bold("\n*liz_coords*")), crayon::green("dataframe with x and y coordinates"),
      crayon::blue(crayon::bold("\n*CA_env*")), crayon::blue("RasterStack with PC environmental layers"),
      crayon::cyan(crayon::bold("\n\n-------------------------------------------------\n"))
    )))
  }
}
