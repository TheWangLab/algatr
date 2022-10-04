#' Helper function to create statistic vector for gt() tables
#'
#' @param stat_name name of statistic
#' @param stat value of statistic
#' @param df dataframe you will eventually add statistic to
#'
#' @return vector with statistic value and name and NA fillers for remaining columns
#' @export
#' @noRd
#'
make_stat_vec <- function(stat_name, stat, df){
  stat_vec <- rep(NA, ncol(df))
  stat_vec[1] <- stat_name
  stat_vec[2] <- stat
  names(stat_vec) <- colnames(df)
  return(stat_vec)
}
