\name{counts.pca}
\alias{counts.pca}
\title{
Principal components analysis of an spectral counts matrix.
}
\description{
A summary and different plots are given as a result of principal components analysis of an spectral counts matrix.
}
\usage{
counts.pca(msnset,facs=NULL,do.plot=TRUE,snms=NULL)
}
\arguments{
  \item{msnset}{
A MSnSet with spectral counts in the expression matrix.
}
  \item{do.plot}{
A logical indicating whether to plot the PCA PC1/PC2 map.
}
  \item{facs}{
NULL or a data frame with factors. See details below.
}
  \item{snms}{
Character vector with sample short names to be plotted. If NULL then 'Xnn' is plotted where 'nn' is the column number in the datset.
}
}
\details{
The spectral counts matrix is decomposed by means of \code{prcomp}.
If \code{do.plot} is TRUE, a plot is generated for each factor showing the PC1/PC2 samples map, with samples colored as per factor level. If \code{facs} is NULL 
then the factors are taken from \code{pData(msnset)}. 
}
\value{
Invisibly returns a list with values:
  \item{pca }{The return value obtained from \code{prcomp}. }
  \item{pc.vars }{The percentage of variability corresponding to each principal component.}
}
\author{
Josep Gregori
}

\seealso{
\code{\link{MSnSet}}, \code{\link{prcomp}}
}

\examples{
data(msms.dataset)
msnset <- pp.msms.data(msms.dataset)
lst <- counts.pca(msnset)
str(lst)
print(lst$pc.vars[,1:4])
}

\keyword{ hplot }
\keyword{ multivariate }
