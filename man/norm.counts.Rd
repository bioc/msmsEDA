\name{norm.counts}
\alias{norm.counts}
\title{
Spectral counts matrix normalization
}
\description{
An spectral counts matrix is normalized by means of a set of samples divisors.
}
\usage{
norm.counts(msnset, div)
}
\arguments{
  \item{msnset}{
A MSnSet with spectral counts in the expression matrix.
}
  \item{div}{
A vector of divisors by sample
}
}
\details{
Each column in the data matrix is divided by the corresponding divisor to obtain the normalizad matrix.
}
\value{
A MSnSet object with the normalized spectral counts.
}

\author{
Josep Gregori
}

\seealso{
The \code{\link{MSnSet}} class documentation and \code{\link{normalize}}
}

\examples{
data(msms.dataset)
msnset <- pp.msms.data(msms.dataset)
(tspc <- apply(exprs(msnset),2,sum))
div <- tspc/median(tspc)
e.norm <- norm.counts(msnset, div)
apply(exprs(e.norm),2,sum)
e.norm
}

\keyword{ manip }
