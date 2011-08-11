\name{flowMerge-package}
\alias{flowMerge-package}
\alias{flowMerge}
\docType{package}
\title{
Merging of mixture components for automated gating of flow cytometry data.
}
\description{
Merges mixture components from the \pkg{flowClust} framework based on the entropy of clustering and provides a simple representation of complicated, non-convex cell populations.
}
\details{
\tabular{ll}{
Package: \tab flowMerge\cr
Type: \tab Package\cr
Version: \tab 0.4.1\cr
Date: \tab 2009-09-07\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}
High density, non-convex cell populations in flow cytometry data often require multiple mixture components for a good model fit. The components are often overlapping, resulting in a complicated representation of individual  cell populations. flowMerge merges overlapping mixture components (based on the max BIC \code{flowClust} model fit) in an iterative manner based on an entropy criterion, allowing these cell populations to be represented by individual mixture components while retaining the good model fitting properties of the BIC solution. Estimates of the number of clusters from a \code{flowMerge} model more accurately represent the "true" number of cell populations in the data.
Running \code{flowMerge} is relatively straightforward. A \code{flowClust} object is converted to a \code{flowObj} object, which groups the model and the data (a \code{flowFrame}) into a single object. This is done by a call to \code{flowObj(model, data)} with a call to \code{merge}, which takes a \code{flowObj} object.
The algorithm may be run in parallel on a multi-core machine or a networked cluster of machines. It uses the functionality in the \code{snow} package to achieve this. Parallelized calls to \code{flowClust} are available via the \code{pFlowClust} and \code{pFlowMerge} functions.

\code{flowMerge} has functionality to automatically select the "correct" number of clusters by fitting a piecewise linear model to the entropy of clustering vs number of clusters, and locating the position of the changepoint. The piecewise linear model fitting is invoked by a call to \code{fitPiecewiseLinreg}, which returns the location of the changepoint.

}
\author{
Greg Finak \email{<greg.finak@ircm.qc.ca>}, Raphael Gottardo \email{<raphael.gottardo@ircm.qc.ca>}

Maintainer: Greg Finak \email{<greg.finak@ircm.qc.ca>}
}
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data; Advances in Bioinformatics (To Appear)}
\keyword{ package }
\seealso{
\code{\link[flowClust]{flowClust},\link{flowObj},\link{pFlowMerge},\link{pFlowClust},\link{fitPiecewiseLinreg},\link{merge},\link{getData},link{plot}}
}
\examples{
#data(rituximab)
#data(RituximabFlowClustFit)
#o<-flowObj(flowClust.res[[which.max(flowMerge:::BIC(flowClust.res))]],rituximab);
#m<-merge(o);
#i<-fitPiecewiseLinreg(m);
#m<-m[[i]];
#plot(m,pch=20,level=0.9);
}
