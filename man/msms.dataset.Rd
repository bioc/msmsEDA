\name{msms.dataset}
\alias{msms.dataset}
\docType{data}
\title{
LC-MS/MS dataset
}
\description{
A MSnSet with a spectral counts matrix as expression and two factors in the
phenoData. \cr
The spectral counts matrix has samples in the columns, and proteins in the rows. \cr
The factors give the treatment and batch conditions of each sample in the dataset.
}
\usage{
data(msms.dataset)
}
\format{
  A MSnSet
}

\references{
Josep Gregori, Laura Villarreal, Olga Mendez, Alex Sanchez, Jose Baselga, 
Josep Villanueva,
"Batch effects correction improves the sensitivity of significance tests in spectral counting-based comparative discovery proteomics."
J Proteomics. 2012 Jul 16;75(13):3938-51. doi: 10.1016/j.jprot.2012.05.005. Epub 2012 May 12.

Laurent Gatto and Kathryn S. Lilley, MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation, Bioinformatics 28(2), 288-289 (2012).
}

\seealso{
See \code{\link{MSnSet}} for detail on the class, and the \code{exprs} and \code{pData} accessors.
}

\examples{
data(msms.dataset)
msms.dataset
dim(msms.dataset)
head(exprs(msms.dataset))
head(pData(msms.dataset))
table(pData(msms.dataset)$treat)
table(pData(msms.dataset)$batch)
table(pData(msms.dataset)$treat, pData(msms.dataset)$batch)
}
\keyword{datasets}
