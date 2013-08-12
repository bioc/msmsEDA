\name{gene.table}
\alias{gene.table}
\title{
Gene symbols associated to protein accessions
}
\description{
Given a character vector with protein accessions, and a character vector with protein descriptions including gene symbols, returns a character vector with gene symbols whose names are the protein accessions. A character pattern should also be given to match the gene symbols.
}
\usage{
gene.table(Accession, Protein, patt = "GN=[A-Z0-9_]*", off = 3)
}
\arguments{
  \item{Accession}{
A character vector with protein accessions
}
  \item{Protein}{
A character vector of protein descriptions including gene name symbols.
}
  \item{patt}{
A character pattern to match the gene symbol within the protein description.
}
  \item{off}{
Offset from the first character in the pattern corresponding to the gene symbol.
}
}
\details{
NA is inserted where no match is found 
}
\value{
A character vector with gene symbols, whose names are the corresponding protein accessions.
}
\author{
Josep Gregori
}

\examples{
data(pnms)
head(pnms)
gene.smb <- gene.table(pnms$Accession,pnms$Proteins)
head(gene.smb)
}
\keyword{ manip }

