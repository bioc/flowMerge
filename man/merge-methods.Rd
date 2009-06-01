\name{merge-methods}
\docType{methods}
\alias{merge-methods}
\alias{merge,flowObj-methods}
\title{ Merge mixture components }
\description{
 Merge mixture components in a \code{flowObj} derived from a \code{flowClust} result and a \code{flowFrame} using the cluster merging algorithm.
}
\section{Methods}{
\describe{

\item{x = "ANY", y = "ANY"}{ The generic method. Should not be called. }

\item{x = "flowObj", y = "missing"}{ The merge method for a \code{flowObj}. }
}}
\value{
    An unnamed list of \code{flowMerge} objects with the \code{k}th element corresponding to the k-cluster merged solution.
}
\examples{
data(rituximab)
data(RituximabFlowClustFit)
o<-flowObj(flowClust.res[[which.max(flowMerge:::BIC(flowClust.res))]],rituximab);
m<-merge(o);
}
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data (Submitted)}
\keyword{methods}
