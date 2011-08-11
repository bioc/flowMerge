\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,flowMerge}
\alias{plot,flowObj}
\title{ Methods for plotting flowMerge and flowObj classes }
\description{
 Plots all possible two-dimensional projections of the parameters in a \code{flowMerge} or \code{flowObj} object and does not require specification of the \code{flowFrame} since a pointer to the data is stored in the object. Informative axis names are used, rather than the usual FL1/FL2/FS/SS channel names. This funciton can take most of the usual additional arguments provided to \code{plot} for the \code{flowClust} package, although some, like the axis names and the \code{data} are fixed. In order for \code{flowMerge} objects to display outliers correctly with \code{plot} (following merging), the \code{updateU} method must be called on them first.
}
\section{Methods}{
\describe{

\item{x = "flowMerge", y = "missing"}{ \code{x} is a \code{flowMerge} object. }

\item{x = "flowObj", y = "missing"}{ \code{x} is a \code{flowObj} object. }
}}
\examples{
#data(rituximab)
#data(RituximabFlowClustFit)
#o<-flowObj(flowClust.res[[which.max(flowMerge:::BIC(flowClust.res))]],rituximab);
#m<-merge(o);
#i<-fitPiecewiseLinreg(m);
#m<-m[[i]];
#plot(m,pch=20,level=0.9);
}
\seealso{\link[flowClust]{flowClust}}
\keyword{methods}
