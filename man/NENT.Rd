\name{NENT}
\alias{NENT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract the Normalized Entropy }
\description{
 Extracts the normalized entropy from a list of flowMerge objects.
}
\usage{
NENT(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A list of flowMerge objects }
}
\details{
 The normalized entropy is extracted from a flowMerge object by computing \eqn{\frac{E}{K*n}} where \eqn{E} is the entropy, and \eqn{K} and n are the number of clusters and data points, respectively. 
}
\value{
    Returns a vector of normalized entropy values for the flowMerge objects.
  }
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data (Submitted)}
\author{ Greg Finak \email{<greg.finak@ircm.qc.ca>} }
\section{Warning}{This function doesn't do enough error checking and will try to extract the entropy from a list of anything.}
\seealso{ \link{ENT},\link{merge},\link{flowMerge-class} }
\examples{
data(RituximabFlowClustFit)
data(rituximab)
o<-flowObj(flowClust.res[[which.max(flowMerge:::BIC(flowClust.res))]],rituximab);
m<-merge(o);
flowMerge:::ENT(m);
flowMerge:::NENT(m);
}
\keyword{misc}