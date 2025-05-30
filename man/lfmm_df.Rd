% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LFMM.R
\name{lfmm_df}
\alias{lfmm_df}
\title{Convert LFMM results into a tidy dataframe for downstream processing}
\usage{
lfmm_df(x)
}
\arguments{
\item{x}{lfmm_test_result element from \code{\link{lfmm_run}} results}
}
\value{
tidy dataframe with LFMM results with each SNP, its p-value, association with env var and other relevant statistics
}
\description{
Convert LFMM results into a tidy dataframe for downstream processing
}
\seealso{
Other LFMM functions: 
\code{\link{lfmm_do_everything}()},
\code{\link{lfmm_manhattanplot}()},
\code{\link{lfmm_qqplot}()},
\code{\link{lfmm_run}()},
\code{\link{lfmm_table}()},
\code{\link{lfmm_test_tidy}()},
\code{\link{quick_elbow}()},
\code{\link{select_K}()}
}
\concept{LFMM functions}
