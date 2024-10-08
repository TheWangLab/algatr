% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LFMM.R
\name{lfmm_do_everything}
\alias{lfmm_do_everything}
\title{LFMM function to do everything}
\usage{
lfmm_do_everything(
  gen,
  env,
  coords = NULL,
  impute = "structure",
  K_impute = 3,
  entropy = TRUE,
  repetitions = 10,
  project = "new",
  quiet_impute = TRUE,
  save_output = FALSE,
  output_filename = NULL,
  K = NULL,
  lfmm_method = "ridge",
  K_selection = "tracy_widom",
  Kvals = 1:10,
  sig = 0.05,
  p_adj = "fdr",
  calibrate = "gif",
  criticalpoint = 2.0234,
  low = 0.08,
  max.pc = 0.9,
  perc.pca = 90,
  max.n.clust = 10,
  quiet = FALSE
)
}
\arguments{
\item{gen}{genotype dosage matrix (rows = individuals & columns = SNPs) or \code{vcfR} object}

\item{env}{dataframe with environmental data or a Raster* type object from which environmental values for the coordinates can be extracted}

\item{coords}{dataframe with coordinates (only needed if K selection is performed with TESS or if environmental values are not provided)}

\item{impute}{if NAs in \code{gen}, imputation will be performed on missing values; options are "structure" which uses the \code{str_impute()} function to impute based on population structure inferred with \code{LEA::snmf} (default); other option is "simple" based on \code{simple_impute()} which imputes to the median}

\item{K_impute}{if \code{impute = "structure"}, an integer vector (range or single value) corresponding to the number of ancestral populations for which the sNMF algorithm estimates have to be calculated (defaults to 3)}

\item{entropy}{
    A boolean value. If true, the cross-entropy criterion is calculated 
    (see \code{\link[LEA]{create.dataset}} and 
    \code{\link[LEA]{cross.entropy.estimation}}).
}

\item{repetitions}{
    An integer corresponding with the number of repetitions for each value of 
    \code{K}.
}

\item{project}{
    A character string among "continue", "new", and "force". If "continue",
    the results are stored in the current project. If "new", the current 
    project is removed and a new one is created to store the result. If 
    "force", the results are stored in the current project even if the input 
    file has been modified since the creation of the project.
}

\item{quiet_impute}{if \code{impute = "structure"}, whether to suppress the results of cross-entropy scores (defaults to TRUE; only does so if K is range of values); only displays run with minimum cross-entropy}

\item{save_output}{if \code{impute = "structure"}, if TRUE, saves SNP GDS and ped (plink) files with retained SNPs in new directory; if FALSE returns object (defaults to FALSE)}

\item{output_filename}{if \code{impute = "structure"} and \code{save_output = TRUE}, name prefix for saved .geno file, SNMF project file, and SNMF output file results (defaults to FALSE, in which no files are saved)}

\item{K}{number of latent factors (if left as NULL (default), K value selection will be conducted)}

\item{lfmm_method}{lfmm method (either \code{"ridge"} (default) or \code{"lasso"})}

\item{K_selection}{method for performing k selection (can either by "tracy_widom" (default), "quick_elbow", "tess", or "find_clusters")}

\item{Kvals}{values of K to test for "tess"}

\item{sig}{alpha level for determining candidate SNPs (defaults to 0.05)}

\item{p_adj}{method to use for p-value correction (defaults to "fdr"); other options can be found in \code{\link{p.adjust}}}

\item{calibrate}{a character string, "gif" or "median+MAD". If the "gif" option is set (default),
significance values are calibrated by using the genomic control method. Genomic control
uses a robust estimate of the variance of z-scores called "genomic inflation factor".
If the "median+MAD" option is set, the pvalues are calibrated by computing the median and MAD of the zscores. If \code{NULL}, the
pvalues are not calibrated.}

\item{criticalpoint}{if \code{K_selection = "tracy_widom"}, a numeric value corresponding to the significance level. If the significance level is 0.05, 0.01, 0.005, or 0.001, the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, respectively (defaults to 2.0234)}

\item{low}{if \code{K_selection = "quick_elbow"}, numeric, between zero and one, the threshold that defines whether a principal component explains 'much' of the variance (defaults to 0.08).}

\item{max.pc}{if \code{K_selection = "quick_elbow"}, maximum percentage of the variance to capture before the elbow (cumulative sum to PC 'n'; defaults to 0.90).}

\item{perc.pca}{if \code{K_selection = "find_clusters"}, a numeric value between 0 and 100 indicating the minimal percentage of the total variance of the data to be expressed by the retained axes of PCA (defaults to 90).}

\item{max.n.clust}{if \code{K_selection = "find_clusters"}, an integer indicating the maximum number of clusters to try. Values of 'k' will be picked up between 1 and max.n.clust (defaults to 10)}

\item{quiet}{whether to operate quietly and suppress the output of tables and figures (defaults to FALSE)}
}
\value{
list with candidate SNPs, model results, and K-value
}
\description{
LFMM function to do everything
}
\details{
LFMM is run using the lfmm package: Jumentier, B. (2021). lfmm: Latent Factor Mixed Models. R package version 1.1.
See also: Caye, K., Jumentier, B., Lepeule, J., & François, O. (2019). LFMM 2: Fast and accurate inference of gene-environment associations in genome-wide studies. Mol. Biol. Evol. 36(4):852-860. doi: https://doi.org/10.1093/molbev/msz008
}
\seealso{
Other LFMM functions: 
\code{\link{lfmm_df}()},
\code{\link{lfmm_manhattanplot}()},
\code{\link{lfmm_qqplot}()},
\code{\link{lfmm_run}()},
\code{\link{lfmm_table}()},
\code{\link{lfmm_test_tidy}()},
\code{\link{quick_elbow}()},
\code{\link{select_K}()}
}
\concept{LFMM functions}
