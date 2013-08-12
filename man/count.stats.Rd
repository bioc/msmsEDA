\name{count.stats}
\alias{count.stats}
\title{
  Summary of statistics of spectral counts by sample in the dataset
}
\description{
Computes the number of proteins identified, the total spectral counts, and a summary of each sample
}
\usage{
count.stats(msnset)
}
\arguments{
  \item{msnset}{
A MSnSet with spectral counts in the expression matrix.
}
}
\value{
A data frame with one row by sample and with variables:
  \item{proteins}{Number of identified proteins in sample}
  \item{counts}{Total spectral counts in sample}
  \item{min}{Min spectral counts}
  \item{lwh}{Tukey's lower hinge spectral counts}
  \item{med}{Median spectral counts}
  \item{hgh}{Tukey's upper hinge spectral counts}
  \item{max}{Max spectral counts}
}
\author{
Josep Gregori
}

\seealso{
\code{\link{MSnSet}}, \code{\link{fivenum}}
}
\examples{
data(msms.dataset)
msnset <- pp.msms.data(msms.dataset)
res <- count.stats(msnset)
res
}
\keyword{univar}
