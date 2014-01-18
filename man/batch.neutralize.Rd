\name{batch.neutralize}
\alias{batch.neutralize}
\title{
Batch effects correction
}
\description{
Computes the SpC matrix where the fixed effects of a blocking factor are 
substracted.
}
\usage{
batch.neutralize(dat, fbatch, half=TRUE, sqrt.trans=TRUE)
}
\arguments{
  \item{dat}{
A SpC matrix with proteins in the rows and samples in the columns.
}
  \item{fbatch}{
A blocking factor of length equal to the number of columns in the expression 
matrix.
}
  \item{half}{
When FALSE, the contrast coefficients are of the contr.treatment style.
When TRUE, the contrast coefficients are of the contr.sum style, its 
aim is to distribute equally the effect to each batch level,
instead of having untouched reference levels.
}
  \item{sqrt.trans}{
When TRUE the fit is done on the square root transformed SpC matrix.  
}
}
\details{
A model with intercept and the blocking factor is fitted. The batch effects
corrected SpC matrix is computed by substracting the estimated effect of the
given blocking factor. When there is no clear reference batch level, the 
default option \code{half=TRUE} should be preferred.
The square root transformation is known to stabilize the variance of Poisson
distributed counts (with variance equal to the mean). The linear model fitting 
gives more accurate errors and p-values on the square root transformed SpC 
matrix. Nevertheless with exploratory data analysis purposes, both the raw and 
square root transformed SpC matrix may give good results.
}
\value{
The batch effects corrected SpC matrix.
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
###  Plot the PCA on the two first PC, and colour by treatment level
ftreat <- pData(msnset)$treat
counts.pca(msnset, facs=ftreat, do.plot=TRUE, snms=as.character(ftreat))
###  Correct the batch effects
spcm <- exprs(msnset)
fbatch <- pData(msnset)$batch
spcm2 <- batch.neutralize(spcm, fbatch, half=TRUE, sqrt.trans=TRUE)
###  Plot the PCA on the two first PC, and colour by treatment level
###  to visualize the improvement.
exprs(msnset) <- spcm2
counts.pca(msnset, facs=ftreat, do.plot=TRUE, snms=as.character(ftreat))
###  Incidence of the correction
summary(as.vector(spcm-spcm2))
plot(density(as.vector(spcm-spcm2)))
}

\keyword{ manip }
