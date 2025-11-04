##' Find the significant eigenvalues of a matrix (from AssocTests).
##'
##' This function implements the Tracy–Widom test to determine how many leading
##' eigenvalues of a matrix are statistically significant. This function and its 
##' documentation are from the \code{tw} function in the archived CRAN package 
##' \emph{AssocTests} (version 1.0-1).
##'
##' @title Tracy–Widom test
##' @description
##' This function is adapted from the \code{tw} function in the archived CRAN
##' package \emph{AssocTests}. The implementation and original documentation
##' were authored by Lin Wang, Wei Zhang, and Qizhai Li. The AssocTests package
##' is archived on CRAN; see the CRAN archive for details.
##'
##' @details
##' The input \code{eigenvalues} should be sorted in descending order. The test
##' statistic is computed for each leading eigenvalue and compared to a supplied
##' Tracy–Widom critical value.
##'
##' @param eigenvalues A numeric vector of eigenvalues (sorted decreasing).
##' @param eigenL The number of eigenvalues.
##' @param criticalpoint A numeric value corresponding to the significance
##' level. If the significance level is 0.05, 0.01, 0.005, or 0.001, set
##' \code{criticalpoint} to \code{0.9793}, \code{2.0234}, \code{2.4224}, or
##' \code{3.2724}, respectively. Default: \code{2.0234}.
##'
##' @return A list with class \code{"htest"}:
##' \tabular{llll}{
##' \code{statistic} \tab \tab \tab a vector of Tracy–Widom statistics.\cr
##' \code{alternative} \tab \tab \tab description of the alternative hypothesis.\cr
##' \code{method} \tab \tab \tab test name.\cr
##' \code{data.name} \tab \tab \tab name of the input data.\cr
##' \code{SigntEigenL} \tab \tab \tab number of significant eigenvalues.\cr
##' }
##'
##'
##' @author
##' Original authors: Lin Wang, Wei Zhang, and Qizhai Li (AssocTests).
##' Documentation adaptation: package maintainers of this project.
##'
##' @references
##' Lin Wang, Wei Zhang, and Qizhai Li. AssocTests: An R Package for Genetic
##' Association Studies. \emph{Journal of Statistical Software}. 2020; 94(5): 1–26.
##'
##' N. Patterson, A. L. Price, and D. Reich. Population Structure and
##' Eigenanalysis. \emph{PLoS Genetics}. 2006; 2(12): 2074–2093.
##'
##' C. A. Tracy and H. Widom. Level-Spacing Distributions and the Airy Kernel.
##' \emph{Communications in Mathematical Physics}. 1994; 159(1): 151–174.
##'
##' A. Bejan. Tracy–Widom and Painlevé II: Computational Aspects and Realisation
##' in S-Plus. In \emph{First Workshop of the ERCIM Working Group on Computing
##' and Statistics}. 2008, Neuchâtel, Switzerland.
##'
##' A. Bejan. Largest eigenvalues and sample covariance matrices. \emph{MSc
##' Dissertation, University of Warwick}. 2005. (This function was written by
##' A. Bejan and made publicly available.)
##'
##' @examples
##' tw(eigenvalues = c(5, 3, 1, 0), eigenL = 4, criticalpoint = 2.0234)
##' @export

tw <- function(eigenvalues, eigenL, criticalpoint = 2.0234)
{
    a <- deparse(substitute(eigenvalues))
    dex <- which(eigenvalues <= 1e-8)
    eigenvalues[dex] <- 1e-8

    L1 <- rev(cumsum(rev(eigenvalues)))
    L2 <- rev(cumsum(rev(eigenvalues^2)))
    N <- eigenL:1
    S2 <- N^2*L2/(L1^2)
    v <- N*(N+2)/(S2-N) # Effective number of markers

    L <- N*eigenvalues/L1

    v.st <- sqrt(v-1)
    N.st <- sqrt(N)

    mu  <- (v.st+N.st)^2/v
    sig <- (v.st+N.st)/v * (1/v.st+1/N.st)^(1/3)

    twstat <- (L - mu) / sig

    d <- which(twstat < criticalpoint)[1]

    if (length(d)==0)
    {
        d <- -100
    }else
    {
        d <- d-1
    }

    structure( 
    list(statistic = c(TW = twstat), 
        alternative = "the eigenvalue is significant", 
        method = "Tracy-Widom test",
        data.name = a,
        SigntEigenL = d
        ), 
    class="htest"
    )
}