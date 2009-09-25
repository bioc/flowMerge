\name{pFlowClust}
\alias{pFlowClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parallelized FlowClust}
\description{
  A parallelized call to flowClust via the snow package and framework. Not called by the user.
}
\usage{
pFlowClust(flowData,cl, K = 1:15, B.init = 100, tol.init = 0.01, tol = 1e-05, B = 1000, randomStart = 50, nu = 4, nu.est = 1, trans = 1, varNames = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flowData}{ The data object, must be a flowFrame, flowSet or list of flowFrames }
  \item{cl}{The snow cluster object}
  \item{K}{ The number of clusters to try for each flowFrame. Can be a vector. This is what is parallelized across processors. }
  \item{B.init}{ See \code{flowClust} documentation }
  \item{tol.init}{ See \code{flowClust} documentation }
  \item{tol}{ See \code{flowClust} documentation }
  \item{B}{ See \code{flowClust} documentation }
  \item{randomStart}{ See \code{flowClust} documentation }
  \item{nu}{ See \code{flowClust} documentation }
  \item{nu.est}{ See \code{flowClust} documentation }
  \item{trans}{ See \code{flowClust} documentation }
  \item{varNames}{ See \code{flowClust} documentation }
}
\details{
  Calls flowClust via the \code{clusterMap} method of the \code{snow} package. Parallelizes the computation of multiple components for a single flowFrame in a loop over multiple flowFrames. If the snow cluster is NULL, will make the call via mapply. 
  }
\value{
  Returns a list of lists of flowClust objects
  The outer list corresponds to the flowFrames passed into the method.
  The inner list corresponds to the \code{K} cluster solutions passed into the method, for each \code{flowFrame}
  (ie If the input is a list of two flowFrames, and K=1:10, then the result is a list of length 2. Each element of the list is itself a list of length 10. The kth element of the inner list is the flowClust k cluster solution.)
 }
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data (Submitted)}
\author{ Greg Finak \email{<greg.finak@ircm.qc.ca>} }
\seealso{ \link[flowClust]{flowClust},\pkg{snow} }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
