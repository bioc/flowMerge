\name{updateU}
\alias{updateU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Update uncertainties}
\description{
  Updates the uncertainties in a \code{flowMerge} ojbect after merging clusters. This function is now internal and no longer exported. It is called automatically within the cluster merging method.
}
\usage{
updateU(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of type \code{flowMerge} }
}
\details{
  Updates the \code{u} slot of the \code{flowMerge} object following merging. The update is computation intensive, and so, is not automatically performed on each \code{flowMerge} object. Should only be done on objects used in further analysis.
}
\value{
A \code{flowMerge} object with the \code{u} slot updated to reflect the new parameter values.
}
\author{ Greg Finak \email{<greg.finak@ircm.qc.ca>}}
\seealso{ \link{flowMerge-class},\link{merge} }
\examples{
#data(rituximab)
#data(RituximabFlowClustFit)
#o<-flowObj(flowClust.res[[which.max(flowMerge:::BIC(flowClust.res))]],rituximab);
#m<-merge(o);
#i<-fitPiecewiseLinreg(m);
#m<-m[[i]];
#plot(m,pch=20,level=0.9);
}
\keyword{ misc }

