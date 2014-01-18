\name{spc.scatterplot}
\alias{spc.scatterplot}
\title{
Scatterplot of SpC means comparing two conditions}
\description{
Given a SpC matrix and a two levels factor, draws a scatterplot with SpC means
of one condition in the x axis and SpC means of the second condition in the y
axis. 
}
\usage{
spc.scatterplot(msms.counts, treat, trans="log2", minSpC=2, minLFC=1, ...)
}
\arguments{
  \item{msms.counts}{
A SpC matrix with proteins in the rows and samples in the columns.
}
  \item{treat}{
A two level factor of length equal to the number of columns in the expression 
matrix. The two levels represent the conditions to be compared.
}
  \item{trans}{
The transformation made on the means before plotting. One among "log2", "sqrt",
or "none". The default is "log2".
}
  \item{minSpC}{
Used as signal threshold.
}
  \item{minLFC}{
Used as size effect threshold.
}
  \item{...}{
Extra parameters passed to the plot function.
}
}
\details{
The transformed means are plotted, one condition versus the other. The borders
representing absolute log fold change 1 are drawn as dashed lines.
All features with log fold change equal to or greather than \code{minLFC} and with 
mean SpC in the most abundant condition equal to or greather than \code{minSpC} 
are colored in red.
}

\author{
Josep Gregori
}

\examples{
data(msms.dataset)
spc.scatterplot(exprs(msms.dataset),treat=pData(msms.dataset)[,1],trans="log2",
       minSpC=2,minLFC=1,main="UPS1 200fm vs 600fm")
}

\keyword{ plots }
