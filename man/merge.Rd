\name{merge}
\alias{merge,flowObj}
\alias{merge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Merge clusters in flow cytometry data }
\description{
 Merge the clusters in a flowClust solution using the cluster merging algorithm and entropy criterion.
}
\usage{
merge(x,y,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A \code{flowObj} object created from a \code{flowClust} object and a \code{flowFrame} using the \code{flowObj} constructor. }
  \item{y}{missing}
  \item{\dots}{Additional arguments, not currently used. }
}
\details{
  Run the cluster merging algorithm on the max BIC solution from a call to \code{flowClust}. 
}
\value{
  A list of unnamed \code{flowMerge} objects. The first element of the list corresponds to the 1--cluster merged solution. The second element corresponds to the 2--cluster merged solution, and so on.
}
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data (Submitted)}
\author{ Greg Finak \email{<greg.finak@ircm.qc.cq>} }

\seealso{ \link[flowClust]{flowClust},\link{flowObj} }
\examples{
data(rituximab)
data(RituximabFlowClustFit)
o<-flowObj(flowClust.res[[which.max(BIC(flowClust.res))]],rituximab)
m<-merge(o);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
