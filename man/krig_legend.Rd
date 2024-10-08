% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TESS.R
\name{krig_legend}
\alias{krig_legend}
\title{Helper function to make a custom legend for TESS maps}
\usage{
krig_legend(gg_df, plot_method, ggplot_fill, minQ)
}
\arguments{
\item{gg_df}{dataframe in tidy format of Q values from \link[algatr]{tess_ggplot}}

\item{plot_method}{method for making rainbow map of kriged layers (options: "maxQ" to only plot the max Q value for each cell (default), "allQ" to plot all Q values greater than \code{minQ}, "maxQ_poly" or "allQ_poly" to create the plots as previously described, but as polygons for each K instead of continuous Q values)}

\item{ggplot_fill}{any ggplot2 scale fill discrete function (default: \link[algatr]{scale_fill_viridis_d}, \code{option = "turbo"})}

\item{minQ}{threshold for minimum Q-value for rainbow plotting if \code{plot_method = "allQ"} or \code{plot_method = "allQ_poly"} is used (defaults to 0.10)}
}
\value{
legend for kriged map from TESS
}
\description{
Helper function to make a custom legend for TESS maps
}
\seealso{
Other TESS functions: 
\code{\link{allK_plot_helper}()},
\code{\link{bestK}()},
\code{\link{ggbarplot_helper}()},
\code{\link{pops_helper}()},
\code{\link{tess_barplot}()},
\code{\link{tess_col_default}()},
\code{\link{tess_do_everything}()},
\code{\link{tess_ggbarplot}()},
\code{\link{tess_ggplot}()},
\code{\link{tess_krig}()},
\code{\link{tess_ktest}()},
\code{\link{tess_plot_allK}()}

Other TESS functions: 
\code{\link{allK_plot_helper}()},
\code{\link{bestK}()},
\code{\link{ggbarplot_helper}()},
\code{\link{pops_helper}()},
\code{\link{tess_barplot}()},
\code{\link{tess_col_default}()},
\code{\link{tess_do_everything}()},
\code{\link{tess_ggbarplot}()},
\code{\link{tess_ggplot}()},
\code{\link{tess_krig}()},
\code{\link{tess_ktest}()},
\code{\link{tess_plot_allK}()}
}
\concept{TESS functions}
\keyword{internal}
