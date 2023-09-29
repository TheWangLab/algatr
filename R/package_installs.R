
#' Install devtools and BiocManager
#'
#' @noRd
dev_packages <- function(){
  if (!require("devtools", quietly = TRUE)) install.packages("devtools")
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
}

#' Install alazygatr packages
#'
#' Checks for the presence of packages required for running alazygatr.
#' If any of these packages are not already installed, it will automatically install them.
#'
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{alazygatr_packages()}
alazygatr_packages <- function(){
  # To install subsets of packages, run the following depending on what methods you want to use
  ## Genetic distance processing:
  genetic_distance_packages()
  ## Genetic data processing:
  data_processing_packages()
  ## Environmental and geographic data processing:
  envirodata_packages()
  ## LFMM:
  lfmm_packages()
  ## RDA:
  rda_packages()
  ## MMRR:
  mmrr_packages()
  ## GDM:
  gdm_packages()
  ## TESS:
  tess_packages()
  ## wingen:
  wingen_packages()
}

#' Install genetic distance packages
#'
#' Checks for the presence of packages required for genetic distance calculations.
#' If any of these packages are not already installed, it will automatically install them.
#'
#' @details
#' The following packages will be installed if not already present:
#' \itemize{
#'   \item "adegenet"
#'   \item "AssocTests"
#'   \item "readr"
#'   \item "tibble"
#'   \item "ecodist"
#'   \item "cowplot"
#' }
#'
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{genetic_distance_packages()}
genetic_distance_packages <- function(){
  if (!require("adegenet", quietly = TRUE)) install.packages("adegenet")
  if (!require("AssocTests", quietly = TRUE)) install.packages("AssocTests")
  if (!require("readr", quietly = TRUE)) install.packages("readr")
  if (!require("tibble", quietly = TRUE)) install.packages("tibble")
  if (!require("ecodist", quietly = TRUE)) install.packages("ecodist")
  if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")
}


#' Install data processing packages
#'
#' Checks for the presence of packages required for genetic data processing.
#' If any of these packages are not already installed, it will automatically install them.
#'
#' @details
#' The following packages will be installed if not already present:
#' \itemize{
#'   \item "gdsfmt" (from Bioconductor repository)
#'   \item "SeqArray" (from Bioconductor repository)
#'   \item "SNPRelate" (from Bioconductor repository)
#' }
#'
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{data_processing_packages()}
data_processing_packages <- function(){
  dev_packages()
  if (!require("gdsfmt", quietly = TRUE)) BiocManager::install("gdsfmt")
  if (!require("SeqArray", quietly = TRUE)) BiocManager::install("SeqArray")
  if (!require("SNPRelate", quietly = TRUE)) BiocManager::install("SNPRelate")
}


#' Install environmental and geographic data processing packages
#'
#' Checks for the presence of packages required for environmental and geographic data processing.
#' If any of these packages are not already installed, it will automatically install them.
#'
#' @details
#' The following packages will be installed if not already present:
#' \itemize{
#'   \item "RStoolbox" (from GitHub repository bleutner/RStoolbox)
#'   \item "geodata"
#'   \item "corrplot"
#'   \item "vegan"
#'   \item "gdistance"
#'   \item "topoDistance"
#'   \item "rmapshaper"
#'   \item "wingen"(from GitHub repository AnushaPB/wingen)
#' }
#'
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{envirodata_packages()}
envirodata_packages <- function(){
  dev_packages()
  if (!require("RStoolbox", quietly = TRUE))devtools::install_github("bleutner/RStoolbox")
  if (!require("geodata", quietly = TRUE)) install.packages("geodata")
  if (!require("corrplot", quietly = TRUE)) install.packages("corrplot")
  if (!require("vegan", quietly = TRUE)) install.packages("vegan")
  if (!require("gdistance", quietly = TRUE)) install.packages("gdistance")
  if (!require("topoDistance", quietly = TRUE)) install.packages("topoDistance")
  if (!require("rmapshaper", quietly = TRUE)) install.packages("rmapshaper")
  if (!require("wingen", quietly = TRUE)) devtools::install_github("wingen")
}


#' Install LFMM packages
#'
#' Checks for the presence of packages required for LFMM.
#' If any of these packages are not already installed, it will automatically install them.
#'
#' @details
#' The following packages will be installed if not already present:
#' \itemize{
#'   \item "adegenet"
#'   \item "AssocTests"
#'   \item "lfmm"
#'   \item "TESS3_encho_sen" (from GitHub repository bcm-uga/TESS3_encho_sen)
#'   \item "LEA" (from Bioconductor repository)
#' }
#'
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{lfmm_packages()}
lfmm_packages <- function(){
  dev_packages()
  if (!require("adegenet", quietly = TRUE)) install.packages("adegenet")
  if (!require("AssocTests", quietly = TRUE)) install.packages("AssocTests")
  if (!require("lfmm", quietly = TRUE)) install.packages("lfmm")
  if (!require("TESS3_encho_sen", quietly = TRUE)) devtools::install_github("bcm-uga/TESS3_encho_sen")
  if (!require("LEA", quietly = TRUE)) BiocManager::install("LEA")
}

#' Install RDA packages
#'
#' Checks for the presence of packages required for RDA.
#' If any of these packages are not already installed, it will automatically install them.
#'
#' @details
#' The following packages will be installed if not already present:
#' \itemize{
#'   \item "ggrepel"
#'   \item "qvalue" (from Bioconductor repository)
#'   \item "robust"
#'   \item "tibble"
#'   \item "vegan"
#' }
#'
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{rda_packages()}
rda_packages <- function(){
  dev_packages()
  if (!require("ggrepel", quietly = TRUE)) install.packages("ggrepel")
  if (!require("qvalue", quietly = TRUE)) BiocManager::install("qvalue")
  if (!require("robust", quietly = TRUE)) install.packages("robust")
  if (!require("tibble", quietly = TRUE)) install.packages("tibble")
  if (!require("vegan", quietly = TRUE)) install.packages("vegan")
}


#' Install GDM packages
#'
#' Checks for the presence of packages required for GDM.
#' If any of these packages are not already installed, it will automatically install them.
#'
#' @details
#' The following packages will be installed if not already present:
#' \itemize{
#'   \item "cowplot"
#'   \item "gdm"
#' }
#'
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{gdm_packages()}
gdm_packages <- function(){
  if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")
  if (!require("gdm", quietly = TRUE)) install.packages("gdm")
}


#' Install MMRR packages
#'
#' Checks for the presence of packages required for MMRR.
#' If the package is not already installed, it will automatically install it.
#'
#' @details
#' The following package will be installed if not already present:
#' - "GGally"
#'
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{mmrr_packages()}
mmrr_packages <- function(){
  if (!require("GGally", quietly = TRUE)) install.packages("GGally")
}


#' Install TESS packages
#'
#' Checks for the presence of packages required for TESS.
#' If any of these packages are not already installed, it will automatically install them.
#'
#' @details
#' The following packages will be installed if not already present:
#' \itemize{
#'   \item "automap"
#'   \item "graphics"
#'   \item "LEA" (from Bioconductor repository)
#'   \item "TESS3_encho_sen" (from GitHub repository bcm-uga/TESS3_encho_sen)
#'   \item "fields"
#'   \item "rworldmap"
#'   \item "cowplot"
#' }
#'
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{tess_packages()}
tess_packages <- function(){
  dev_packages()
  if (!require("automap", quietly = TRUE)) install.packages("automap")
  if (!require("graphics", quietly = TRUE)) install.packages("graphics")
  if (!require("LEA", quietly = TRUE)) BiocManager::install("LEA") # required by tess3r
  if (!require("TESS3_encho_sen", quietly = TRUE)) devtools::install_github("bcm-uga/TESS3_encho_sen")
  if (!require("fields", quietly = TRUE)) install.packages("fields")
  if (!require("rworldmap", quietly = TRUE)) install.packages("rworldmap")
  if (!require("cowplot", quietly = TRUE)) install.packages("cowplot")
}

#' Install wingen packages
#'
#' This function checks for the presence of packages required for wingen.
#' If the packages are not already installed, it will automatically install them.
#'
#' @details
#' The following package will be installed if not already present:
#' - "wingen" (from GitHub repository AnushaPB/wingen)
#'
#' @return None
#'
#' @export
#'
#' @examples
#' \dontrun{wingen_packages()}
wingen_packages <- function(){
  if (!require("wingen", quietly = TRUE)) devtools::install_github("AnushaPB/wingen")
  }
