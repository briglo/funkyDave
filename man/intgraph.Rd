% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.r
\name{intgraph}
\alias{intgraph}
\title{intgraph}
\usage{
intgraph(interactionData, annotation, scoreCut = 0.3, numberCut = 20,
  numberSplit = 60)
}
\arguments{
\item{interactionData}{a dataframe with source target and countAboveMean[numeric]}

\item{annotation}{a dataframe of fields name, type, ID, count, optional x+y coordinates for manual graph layout}

\item{scoreCut}{a numeric value pointing to a column of interactionData (default 0.3)}

\item{numberCut}{numeric value for minimum number of interactions to plot (default 20)}

\item{numberSplit}{numeric value for spliting numbers of interactions (i.e high vs low, default 60)}
}
\value{
a plot, invisibly the graph
}
\description{
takes a summarized cellphoneDB object and plots it better than my other attempts
}
\examples{
data("190801_newNetworkGraph")
IG<-intgraph(interactionData=cellphoneDB_data,annotation=cluster_anno,scoreCut=0.3, numberCut=20, numberSplit=60)


}
