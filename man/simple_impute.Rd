% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/imputation.R
\name{simple_impute}
\alias{simple_impute}
\title{Impute NA values
NOTE: use extreme caution when using this form of simplistic imputation. We mainly provide this code for creating test datasets and highly discourage its use in analyses.}
\usage{
simple_impute(x, FUN = median)
}
\arguments{
\item{x}{matrix}

\item{f}{function to use for imputation (defaults to median)}
}
\value{
matrix of values with missing values imputed
}
\description{
Impute NA values
NOTE: use extreme caution when using this form of simplistic imputation. We mainly provide this code for creating test datasets and highly discourage its use in analyses.
}
\seealso{
Other Imputation functions: 
\code{\link{gen_to_geno}()},
\code{\link{geno_to_dosage}()},
\code{\link{snmf_bestK}()},
\code{\link{snmf_crossent_helper}()},
\code{\link{str_impute}()}
}
\concept{Imputation functions}
