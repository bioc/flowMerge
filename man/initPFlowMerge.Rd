\name{initPFlowMerge}
\alias{initPFlowMerge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Initialize a SNOW cluster for use with flowMerge }
\description{
  Initializes a snow cluster for use with flowMerge, ensures that the flowMerge library is loaded in all environments. Not meant to be called by the user
}
\usage{
initPFlowMerge(cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl}{ A \code{snow} cluster }
}
\details{
  A valid \code{snow} cluster.
}
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data (Submitted)}
\author{ Greg Finak \email{<greg.finak@ircm.qc.ca>} }

\seealso{ \link{pFlowMerge} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
