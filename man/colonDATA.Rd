\docType{data}
\name{colonDATA}
\alias{colonDATA}
\title{Single-cell RNA-Seq data of normal colon and colon cancer cells}
\format{
A total of 3 objects containing the normalized single-cell RNA-Seq data
(Fluidigm C1) of 432 normal colon and colon cancer cells, the
normal/cancer status information, and the colon-specific regulatory
network inferred from GTEX.
}
\usage{
data(colonDATA)
}
\description{
A total of 3 objects containing the normalized single-cell RNA-Seq data
(Fluidigm C1) of 432 normal colon and colon cancer cells, the
normal/cancer status information, and the colon-specific regulatory
network inferred from GTEX.
}
\details{
When loaded in, there are three separate objects:
\itemize{
  \item \code{avlfpkmEID.m}: the matrix representing the normalized
  expression values for 432 cells (columns), with rows labeling genes
  annotated to Entrez gene IDs. 
  \item \code{statusNC.v}: a vector specifying whether the cell is
  normal (0) or cancer (1).
  \item \code{netCOL.m}: the colon-specific regulon matrix with columns
  labeling 56 TFs, and rows labeling gene targets, with entries either
  1, -1 or 0, to indicate positive, negative and no-interactions, respectively.
}

}
\references{
Teschendorff AE, Wang, N.
\emph{Improved detection of tumor suppressor events in single-cell
  RNA-Seq data}
Submitted.
}
\keyword{datasets}
