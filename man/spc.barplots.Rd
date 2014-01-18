\name{spc.barplots}
\alias{spc.barplots}
\title{
Set of SpC barplots by sample}
\description{
Draws bars of height proportional to the sample size of each column  
in a SpC matrix. The sizes are scaled to the median of the total SpC by sample.
}
\usage{
spc.barplots(msms.counts,fact=NULL,...)
}
\arguments{
  \item{msms.counts}{
A SpC matrix with proteins in the rows and samples in the columns.
}
  \item{fact}{
NULL or a factor of length equal to the number of columns in the expression 
matrix. If provided the bars are colored by factor level.
}
  \item{...}{
Extra parameters passed to the plot function.
}
}
\details{
.
}

\author{
Josep Gregori
}

\examples{
data(msms.dataset)
spc.barplots(exprs(msms.dataset),fact=pData(msms.dataset)[,1],
      main="UPS1 200fm vs 600fm")
}

\keyword{ plots }
