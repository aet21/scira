% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sciraEstRegAct.R
\name{sciraEstRegAct}
\alias{sciraEstRegAct}
\title{Estimates Regulatory Activity}
\usage{
sciraEstRegAct(data, regnet.m, norm = c("z", "c"), ncores = 4)
}
\arguments{
\item{data}{Either a vector with names labeling genes (Entrez gene IDs) or a matrix with rows labeling genes (Entrez gene IDs) and columns labeling samples (cells).}

\item{regnet.m}{The regulon matrix with rows labeling genes (Entrez IDs) and columns labeling the transcription factors (also Entrez gene IDs) and entries either 1 (positive), -1 (inhibitory) or 0 (no association). Typically, the output of either 'sciraSelReg' or 'sciraInfReg'}

\item{norm}{Type of normalization to use if data is provided as a matrix. "z" stands for z-score transformation, "c" just row centres the data matrix.}

\item{ncores}{Number of parallel cores to use}
}
\value{
A matrix of regulatory activity scores, with rows labeling the transcription factors and columns labeling samples (cells).
}
\description{
This function takes as input a set of regulons, and a matrix of RNA-Seq profiles (single cell or bulk) and
outputs a regulatory activity matrix over the transcription factors and samples
}
\examples{
pheno.v <- c(rep(1,50),rep(2,50));
data.m <- matrix(rnorm(10000,0,1),nrow=1000,ncol=100);
data.m[1:100,1:50] <- matrix(rnorm(5000,2,1),nrow=100,ncol=50);
out.o <- LimmaFn(pheno.v,data.m);

}
