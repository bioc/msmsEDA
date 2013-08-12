\name{counts.heatmap}
\alias{counts.heatmap}

\title{
Heatmap of an spectral counts matrix.
}
\description{
Heatmap showing the clustering of proteins and samples in a matrix of spectral counts
}
\usage{
counts.heatmap(msnset,etit=NULL,fac=NULL,to.pdf=FALSE)
}

\arguments{
  \item{msnset}{
A MSnSet with spectral counts in the expression matrix.
}
  \item{etit}{
The root name of the pdf file names where the heatmaps are sent.
}
  \item{fac}{
A factor which is used for the column color bar.
}
  \item{to.pdf}{
A logical indicating whether the heatmaps are sent to a pdf file.
}
}
\details{
A heatmap of the \code{msnset} expression matrix is plot.
If \code{to.pdf} is TRUE two heatmaps are plot, the first is fitted on an A4 
page, the second is plotted with 3mm by row, allocating enough height to make
the rownames readable.
If \code{fac} is not NULL then a column color bar will show the levels 
of the factor.
If \code{to.pdf} is TRUE the heatmaps are sent to pdf files whose names
are the concatenation of \code{etit} and "-HeatMap.pdf" and "-FullHeatMap.pdf",
otherwise \code{etit} has no effect.
}
\value{
No value is returned
}

\author{
Josep Gregori
}

\seealso{
\code{\link{MSnSet}}, \code{\link{heatmap}} and \code{\link{heatmap.2}}
}
\examples{
data(msms.dataset)
msnset <- pp.msms.data(msms.dataset)
counts.heatmap(msnset,fac = pData(msnset)$treat)
}
\keyword{ hplot }
\keyword{ multivariate }

