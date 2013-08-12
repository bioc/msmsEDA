\name{msmsEDA-package}
\alias{msmsEDA-package}
\alias{msmsEDA}
\docType{package}
\title{
Exploratory Data Analysis of label-free LC-MS/MS spectral counts
}
\description{
Exploratory data analysis to assess the quality of a set of label-free LC-MS/MS experiments, quantified by spectral counts, and visualize de influence of the involved factors.
Visualization tools to discover outliers and eventual confounding.
}
\details{
\tabular{ll}{
Package: \tab msmsEDA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-05-24\cr
License: \tab GPL-2\cr
}

\tabular{ll}{
\code{pp.msms.data}: \tab  {data preprocessing} \cr
\code{gene.table}: \tab  {extract gene symbols from protein description} \cr
\code{count.stats}: \tab  {summaries by sample} \cr
\code{counts.pca}: \tab  {principal components analysis} \cr
\code{counts.hc}: \tab  {hierarchical clustering of samples} \cr
\code{norm.counts}: \tab  {normalization of spectral counts matrix} \cr
\code{counts.heatmap}:  \tab{experiment heatmap} \cr
\code{disp.estimates}: \tab  {dispersion analysis and plots} \cr
}
}
\author{
Josep Gregori, Alex Sanchez and Josep Villanueva \cr
Maintainer: Josep Gregori <gregori@aiqs.es>
}
\references{
Gregori J, Villarreal L, Mendez O, Sanchez A, Baselga J, Villanueva J,
"Batch effects correction improves the sensitivity of significance tests in spectral counting-based comparative discovery proteomics."
J Proteomics. 2012 Jul 16;75(13):3938-51. doi: 10.1016/j.jprot.2012.05.005. Epub 2012 May 12.
}
\keyword{ package }
\keyword{ multivariate }
\keyword{ cluster }
\keyword{ hplot }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
%\examples{
%~~ simple examples of the most important functions ~~
%}
