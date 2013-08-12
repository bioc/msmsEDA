\name{pnms}
\alias{pnms}
\docType{data}
\title{
Accessions and gene symbols
}
\description{
A data frame with accessions in one column, and protein description including gene symbols in the second column.
}
\usage{data(pnms)}
\format{
  A data frame with 1160 observations on the following 2 variables.
  \describe{
    \item{\code{Accession}}{a character vector with the protein accessions}
    \item{\code{Proteins}}{a character vector with a description of each protein, including the gene symbol}
  }
}
\examples{
data(pnms)
str(pnms)
head(pnms)
}
\keyword{datasets}
