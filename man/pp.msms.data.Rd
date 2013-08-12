\name{pp.msms.data}
\alias{pp.msms.data}
\title{ Spectral counts matrix pre-processing }
\description{
Given a MSnSet, possibly subsetted from a bigger dataset, removes the all zero rows, and those whith row names (accessions) ending with '-R' in the corresponding
expression matrix. NAs are replaced by zeroes, as usually a NA in a spectral counts matrix 
corresponds to a proteint not identified in a sample.
}
\usage{
pp.msms.data(msnset)
}
\arguments{
  \item{msnset}{
A MSnSet with spectral counts in the expression matrix.
}
}
\details{
An '-R' protein corresponds to an artefactual identification. \cr
Rows with all zeros are uninformative and may give rise to errors in the 
analysis.\cr
A NA is understood as a unidintified protein in a sample.
}
\value{
Returns an updated MSnSet object.\cr
Its processingData slot shows that the object has been processed by
\code{pp.msms.data}
}
\author{
Josep Gregori
}

\seealso{
\code{\link{MSnSet}}
}

\examples{
data(msms.dataset)
dim(msms.dataset)
msnset <- pp.msms.data(msms.dataset)
dim(msnset)
}
\keyword{ manip }
\keyword{ array }% __ONLY ONE__ keyword per line
