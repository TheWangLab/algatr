# Internal environment to track warnings (not exported)
.my_pkg_env <- new.env(parent = emptyenv())

# Function to emit a DPS warning once per session
dps_warning <- function() {
  if (isTRUE(getOption("wingen.quiet_dps_warning", FALSE))) {
    return(invisible(NULL))  # User opted out of warning
  }

  if (!isTRUE(.my_pkg_env$dps_warning)) {
    warning(
      paste0(
        "Prior to May 2025, this function incorrectly returned the proportion of shared alleles (PS) ",
        "instead of the genetic distance measure: DPS = 1 - PS. ",
        "Please review results from prior versions accordingly. ",
        "This warning will appear once per session. ",
        "To suppress, set options(wingen.quiet_dps_warning = TRUE)."
      ),
      call. = FALSE
    )
    .my_pkg_env$dps_warning <- TRUE
  }
}
