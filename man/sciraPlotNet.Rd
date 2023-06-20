\name{sciraPlotNet}
\alias{sciraPlotNet}
\title{Plots Regulatory Network}
\usage{
sciraPlotNet(net)
}
\arguments{
\item{net}{The regulatory network in matrix format with rows labeling
  gene targets (Entrez gene IDs) and columns labeling transcription
  factors (Entrez gene IDs). Entries are either 1, -1 or 0, indicating
  positive, negative and no interactions, respectively}
}
\value{
A ggplot2 object
}
\description{
This function provides visualisation of the regulons
}
\examples{
data(colonDATA);
netplot <- sciraPlotNet(netCOL.m);
netplot
}
