% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/supermap.R
\name{getmetacell_RNA}
\alias{getmetacell_RNA}
\title{getmetacell_RNA}
\usage{
getmetacell_RNA(object, cluster.name, size_factor = 10^6)
}
\arguments{
\item{object}{A SeuratObject}

\item{cluster.name}{The cluster name in the SeuratObject}

\item{size}{size factor for normalization}
}
\value{
A gene expression matrix, where rows represent cells and columns represent genes
}
\description{
This function allows us to obtain the gene expression matrix of metacells
}
