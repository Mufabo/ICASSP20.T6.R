% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{my_kmeanspp}
\alias{my_kmeanspp}
\title{Kmeans with kmeans++ initialization and cityblock distance}
\usage{
my_kmeanspp(data, ll, iter = 25, rep = 50)
}
\arguments{
\item{data}{matrix. Data}

\item{ll}{int. Number of clusters}

\item{iter}{int. Number of iterations}

\item{rep}{int. Number of replications}
}
\value{
list
\enumerate{
\item centroids matrix. Cluster centers
\item clusters vector. Cluster membership of data samples
}
}
\description{
Kmeans with kmeans++ initialization and cityblock distance
}
