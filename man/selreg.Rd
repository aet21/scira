\docType{data}
\name{selreg}
\alias{selreg}
\title{Liver-specific transcription factors and their regulons}
\format{
A list containing the liver specific regulons, summary information and
top-ranked differentially expressed genes between liver and all other
tissues and between liver and blood.
}
\usage{
data(selreg)
}
\description{
The liver-specific regulons are specified as a matrix with columns
labeling 22 liver-specific TFs and rows labeling gene targets.
}
\details{
The list object \code{selreg.o} contains the following elements:
\itemize{
  \item \code{netTOI}: the matrix specifying the regulons with columns
  labeling the 22 liver-specific TFs and rows labeling the gene-targets,
  with entries either 1, -1 or 0 for positive, negative and no
  interactions.
  \item \code{sumnet}: a matrix with rows labeling the 22 liver-specific
  TFs and columns labeling the number of gene-targets in regulon, and
  their distribution in terms of positive and negative interactions.
  \item \code{top}: a list of two matrices ranking the top
  differentially expressed genes between liver tissue and all other
  tissues (first element) and between liver and blood tissue (2nd
  element). The latter comparison is done to remove the confounding
  effect of immune-cell infiltrates in liver tissue.
}
}
\references{
Teschendorff AE, Wang, N.
\emph{Improved detection of tumor suppressor events in single-cell
  RNA-Seq data}
Submitted.
}
\keyword{datasets}
