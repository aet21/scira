\docType{data}
\name{PrtAtlasExp}
\alias{PrtAtlasExp}
\title{Protein Atlas Normal Tissue RNA-Seq dataset}
\format{
A matrix representing the normalized RNA-Seq data of normal tissue-types
from the Protein Atlas project.
}
\usage{
data(PrtAtlasExp)
}
\description{
The rownames of the matrix label the Entrez gene IDs, wheras the columns
label the samples representing normal tissue-types. Data has already been normalized and log-transformed.
}
\details{
\itemize{
  \item \code{PrtAtlasExp.m}: the matrix representing the normalized
  expression values for normal tissue-types of the Protein Atlas
  project. The column names of the matrix label the tissue-type.
}
}
\references{
Teschendorff AE, Wang, N.
\emph{Improved detection of tumor suppressor events in single-cell
  RNA-Seq data}
Submitted.
}
\keyword{datasets}
