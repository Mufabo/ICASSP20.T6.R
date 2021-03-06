% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_31.R
\name{data_31}
\alias{data_31}
\title{Creates three gaussian clusters with replacement outliers}
\usage{
data_31(N_k, epsilon)
}
\arguments{
\item{N_k}{int. Number of samples in each cluster.}

\item{epsilon}{scalar. Percentage of outliers}
}
\value{
list
\enumerate{
\item data matrxi[3*N_k, 2] Samples
\item labels vector[N] Cluster memberships of the samples
\item r int. Number of dimensions in generated data. Set to 2
\item N int. Total number of samples. Set to 3*N_k
\item K_true int. Number of clusters. Set to 3
\item mu_true matrix[r, K_true] Cluster means
\item scatter_true array[r, r, K_true] Cluster scatter matrices
}
}
\description{
Creates three gaussian clusters with replacement outliers
}
\note{
"Robust M-Estimation based Bayesian Cluster Enumeration for Real Elliptically Symmetric Distributions"
Christian A. Schroth and Michael Muma, Signal Processing Group, Technische UniversitÃ¤t Darmstadt
submitted to IEEE Transactions on Signal Processing
}
\examples{

}
