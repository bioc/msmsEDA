\name{filter.flags}
\alias{filter.flags}
\title{
Flag proteins with a minimum signal and/or sufficient dispersion.
}
\description{
In general the spectral counts (SpC) matrix of a LC-MS/MS experiment is a 
sparse matrix, where
most of the features have very low signal. Besides, the features with low 
variance to mean ratio (dispersion) will be scarcely informative in a 
biomarker discovery experiment. Given a minimum number of spectral counts 
and/or a fraction of the features to be excluded by low dispersion, this 
function returns a vector of logicals flagging all features with values
above the given thresholds.
}
\usage{
filter.flags(data,minSpC=2,frac.out=0.4)
}
\arguments{
  \item{data}{
A SpC matrix with proteins in the rows and samples in the columns.
}
  \item{minSpC}{
All features with SpC below this threshold will be flagged as FALSE.
}
  \item{frac.out}{
The fraction of features to be excluded, with the lowest observed dispersion.
These will be flagged as FALSE.
}
}
\details{
The less informative features in a SpC matrix are flagged as FALSE. Those with
high enough signal and dispersion are flagged as TRUE. This vector of logicals
may be used to filter the SpC matrix which is used in plots where only the
relevant informattion matters, and where the high number of 0 may distort the
plot and difficult its interpretation.
}
\value{
A vector of logical values.
}

\author{
Josep Gregori
}

\examples{
data(msms.dataset)
fraction <- 0.3
msnset <- pp.msms.data(msms.dataset)
flags <- filter.flags(exprs(msnset),minSpC=2,frac.out=fraction)
cat("\nNumber of informative features:",sum(flags),"\n")
}

\keyword{ manip }
