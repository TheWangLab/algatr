
#' Function to Select K value
#'
#' @param gen genotype matrix
#' @param k_selection method for performing K selection: "tracy_widom", "quick_elbow", or "tess"
#' @param coords coordinates for TESS-based K selection
#'
#' @return
#' @export
#'
#' @examples
get_K <- function(gen, coords = NULL, k_selection = "tracy_widom", Kvals = Kvals, ...){

  if(k_selection == "tracy_widom"){K <- get_K_tw(gen)}

  if(k_selection == "quick_elbow"){K <- get_K_elbow(gen)}

  if(k_selection == "tess"){K <- get_K_tess(gen, coords, Kvals)}

  return(K)
}

#' @describeIn get_K Determine K using Tracy-Widom Test
#' @param gen genotype matrix
#'
#' @return
#' @export
#'
#' @examples
get_K_tw <- function(gen){
  # Run PCA
  pc <- prcomp(gen)

  # Get eig
  eig <- pc$sdev^2

  # Run Tracy-Widom test
  # NOTE:
  # The critical point is a numeric value corresponding to the significance level.
  # If the significance level is 0.05, 0.01, 0.005, or 0.001,
  # the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, accordingly.
  # The default is 2.0234.
  tw_result <- AssocTests::tw(eig, eigenL = length(eig), criticalpoint = 0.9793)

  # Get K based on the number of significant eigenvalues
  K <- tw_result$SigntEigenL

  return(K)
}

#' @describeIn get_K Determine K using PCA and \code{quick.elbow}
#'
#' @param gen genotype matrix
#'
#' @return
#' @export
#'
#' @examples
get_K_elbow <- function(gen){
  # Run PCA
  pc <- prcomp(gen)

  # Get eig
  eig <- pc$sdev^2
  # Estimate number of latent factors using quick.elbow (see general functions for description of how this function works)
  # This is a crude way to determine the number of latent factors that is based on an arbitrary "low" value
  K <- quick.elbow(eig, low = 0.08, max.pc = 0.9)

  par(pty = "s",mfrow = c(1,1))
  plot(eig, xlab = 'PC', ylab = "Variance explained")
  abline(v = K, col = "red", lty = "dashed")

  return(K)
}


#' @describeIn get_K Determine K using TESS and \code{bestK}
#'
#' @param gen genotype matrix
#' @param coords coordinates for "tess"
#' @param Kvals values of K to test for "tess"
#' @param tess_method method to use for "tess"
#' @param ploidy ploidy for "tess"
#'
#' @return
#' @export
#'
#' @examples
get_K_tess <- function(gen, coords, Kvals = 1:10, tess_method = "projected.ls", ploidy = 2){
  # Run tess for all K values
  tess3_obj <- tess3r::tess3(X = gen, coord = coords, K = Kvals, method = tess_method, ploidy = ploidy)

  # Plot CV results and mark the K-value that is automatically selected
  plot(tess3_obj, pch = 19, col = "blue",
       xlab = "Number of ancestral populations",
       ylab = "Cross-validation score")

  # Get best K value
  K <-  bestK(tess3_obj, Kvals)

  return(K)
}
