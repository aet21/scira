\docType{data}
\name{scdataYML}
\alias{scdataYML}
\title{Single-cell RNA-Seq data of a liver differentiation timecourse}
\format{
A total of 3 objects containing the normalized single-cell RNA-Seq data
(Fluidigm C1) of a timecourse differentiation of hepatoblasts into
hepatocytes and cholangiocytes. A total of 447 cells over 7 timepoints
(mouse embryonic day to mouse embryonic day 17) were measured.

}
\usage{
data(scdataYML)
}
\description{
A total of 3 objects containing the normalized single-cell RNA-Seq data
(Fluidigm C1) of a timecourse differentiation of hepatoblasts into
hepatocytes and cholangiocytes. A total of 447 cells over 7 timepoints
(mouse embryonic day to mouse embryonic day 17) were measured.
}
\details{
When loaded in, there are three separate objects:
\itemize{
  \item \code{avlscHEIDyml.m}: the matrix representing the normalized
  expression values for 447 cells (columns), with rows labeling genes
  annotated to Entrez gene IDs. Note that although the experiment was
  performed in mice, that we provide the data matrix mapped to human
  homolog genes, since the GTEX regulons we use were derived from human tissues.
  \item \code{timept.v}: a vector specifying the distinct embryonic day timepoints at
  which cells were measured.
  \item \code{phenoYML.v}: a vector specifying the timepoint (indexed as
  a number between 1 and 7) of the 447 cells. The index increases in
  ascending order in line with the timepoint as specified by \code{timept.v}.
}
}
\references{
Teschendorff AE, Wang, N.
\emph{Improved detection of tumor suppressor events in single-cell
  RNA-Seq data}
Submitted.
}
\keyword{datasets}
