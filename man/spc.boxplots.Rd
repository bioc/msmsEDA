\name{spc.boxplots}
\alias{spc.boxplots}
\title{
Set of SpC boxplots by sample}
\description{
Draws a boxplot for each column (sample) in a SpC matrix. 
The SpC are previosly transformed by log2, with an offset of 0.1.
If a factor is provided the boxplots are colored by factor level to better
visualize the differences. 
}
\usage{
spc.boxplots(msms.counts,fact=NULL,minSpC=2,...)
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
matrix. If provided the boxplots are colored by factor level.
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
spc.boxplots(exprs(msms.dataset),fact=pData(msms.dataset)[,1],
      main="UPS1 200fm vs 600fm")
}

\keyword{ plots }
