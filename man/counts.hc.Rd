\name{counts.hc}
\alias{counts.hc}
\title{
Hierarchical clustering on an spectral counts matrix.
}
\description{
Hierarchical clustering of samples in an spectral counts matrix, coloring tree branches according to factor levels.
}
\usage{
counts.hc(msnset, do.plot = TRUE, facs = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{msnset}{
A MSnSet with spectral counts in the expression matrix.
}
  \item{do.plot}{
A logical indicating whether to plot the dendrograms.
}
  \item{facs}{
NULL, or a data frame with factors. See details below.
}
}

\details{
The hierarchical clustering is done by means of \code{hclust} with default parameters.
If \code{do.plot} is TRUE, a dendrogram is plotted for each factor, with branches colored as per factor level. If \code{facs} is NULL then the factors are taken 
from \code{pData(msnset)}.
}
\value{
Invisibly returns the the value obtained from \code{hclust}.
}

\author{
Josep Gregori
}

\seealso{
\code{\link{MSnSet}}, \code{\link{hclust}}
}

\examples{
data(msms.dataset)
msnset <- pp.msms.data(msms.dataset)
hc <- counts.hc(msnset)
str(hc)
}

\keyword{ hplot }
\keyword{ multivariate }
