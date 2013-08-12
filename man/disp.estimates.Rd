 \name{disp.estimates}
\alias{disp.estimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Residual dispersion estimates
}
\description{
Estimates the residual dispersion of each row of a spectral counts matrix as the ratio residual variance to mean of mean values by level, for each factor in 
\code{facs}. Different plots are drawn to help in the interpretation of the results.
}
\usage{
disp.estimates(msnset, facs=NULL, do.plot = TRUE, etit = NULL, to.pdf=FALSE)
}
\arguments{
  \item{msnset}{
A MSnSet with spectral counts in the expression matrix.
}
  \item{facs}{
A factor or a data frame with factors.
}
  \item{do.plot}{
A logical indicating whether to produce dispersion distribution plots.
}
  \item{etit}{
Root name of the pdf file where to send the plots.
}
  \item{to.pdf}{
A logical indicating whether a pdf file should be produced.
}
}
\details{
Estimates the residual dispersion of each protein in the spectral counts matrix, for each factor in \code{facs}, and returns the quantiles at \code{c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)} of the distribution of dispersion values for each factor. If \code{facs} is NULL the factors are taken from \code{pData(msnset)}. If \code{do.plot} is TRUE this function produces a density plot of dispersion values, and the scatterplot of residual variance vs mean values, in log10 scale. If \code{do.pdf} is TRUE \code{etit} provides the root name for the pdf file name, ending with "-DispPlots.pdf". If \code{etit} is NULL a default value of "MSMS" is provided. A different set of plots is produced for each factor in \code{facs}.
}
\value{
Silently returns a matrix with the quantiles at \code{c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)} of the residual dispersion estimates. Each row has the residual dispersion values attribuable to each factor in \code{facs}.
}
\author{
Josep Gregori
}

\examples{
data(msms.dataset)
msnset <- pp.msms.data(msms.dataset)
disp.q <- disp.estimates(msnset)
disp.q
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
\keyword{ distribution }
