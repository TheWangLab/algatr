% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gen_dist.R
\name{gen_dist_corr}
\alias{gen_dist_corr}
\title{Plot the relationship between two distance metrics}
\usage{
gen_dist_corr(dist_x, dist_y, metric_name_x, metric_name_y)
}
\arguments{
\item{dist_x}{df containing square distance matrix for x axis}

\item{dist_y}{df containing square distance matrix for y axis}

\item{metric_name_x}{name of distance metric for x axis; if DPS used, must be \code{"dps"}}

\item{metric_name_y}{name of distance metric for y axis; if DPS used, must be \code{"dps"}}
}
\value{
scatterplot comparing two user-defined genetic distance metrics
}
\description{
Plot the relationship between two distance metrics
}
