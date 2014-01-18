\name{spc.densityplots}
\alias{spc.densityplots}
\title{
SpC density plots of a SpC matrix}
\description{
Draws superposed density plots, one for each column (sample) in a SpC matrix. 
The SpC are previosly transformed by log2, with an offset of 0.1.
If a factor is provided the density curves are colored by factor level to better
visualize the differences. 
}
\usage{
spc.densityplots(msms.counts,fact=NULL,minSpC=2,...)
}
\arguments{
  \item{msms.counts}{
A SpC matrix with proteins in the rows and samples in the columns.
}
  \item{minSpC}{
All matrix cells with values below this threshold are excluded.
}
  \item{fact}{
NULL or a factor of length equal to the number of columns in the expression 
matrix. If provided the density curves are colored by factor level.
}
  \item{...}{
Extra parameters passed to the plot function.
}
}
\details{
More informative plots are obtained when excluding the cells with values below
2, the default for \code{minSpC}.
}

\author{
Josep Gregori
}

\examples{
data(msms.dataset)
spc.densityplots(exprs(msms.dataset),fact=pData(msms.dataset)[,1],
      main="UPS1 200fm vs 600fm")
}

\keyword{ plots }
