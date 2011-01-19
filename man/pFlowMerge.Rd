\name{pFlowMerge}
\alias{pFlowMerge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Parellel call to flowMerge }
\description{
  Calls the flowMerge methods to compute the merged solution from a flowClust object or set of objects in a parallelized manner using the \code{snow} framework.
}
\usage{
pFlowMerge(flowData, cl, K = 1:15, B.init = 100, tol.init = 0.01, tol = 1e-05, B = 500, randomStart = 10, nu = 4, nu.est = 0, trans = 1, varNames = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flowData}{ The data to be fit. A list of \code{flowFrame}s, a \code{flowSet} or a \code{flowFrame} }
  \item{cl}{ The \code{snow} cluster object. Can be NULL to call the non-parallel version of flowClust }
  \item{K}{ See \code{flowClust} documentation }
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
  Makes a parallelized call to \code{flowClust}. Parses the results to extract the max BIC solution, merges clusters, finds the optimal k-cluster solution using the entropy and returns it. If \code{cl} is \code{NULL}, a non-parallel call is made to the flowClust function.
}
\value{
  A list of flowMerge objects. One per flowFrame passed into the method.
  }
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data (Submitted)}
\author{ Greg Finak \email{<greg.finak@ircm.qc.ca>}}
\section{Warning}{This function does not do any special memory management. A large data set will likely cause it to run out of memory and start swapping incessantly. If you have lots of data, it's best to feed it piecewise to pFlowClust.
}
\seealso{
  \link{pFlowClust},\link[flowClust]{flowClust},\link{merge},\pkg{snow},
  \link{fitPiecewiseLinreg}}  
\examples{
library(multicore)
data(rituximab)
#Parallelized call below:
\dontrun{cl<-makeSOCKcluster(rep("finakg@localhost",7))}
\dontrun{result<-pFlowMerge(rituximab,cl,varNames=c("FSC.H","SSC.H"))}
\dontrun{plot(result)}
cl<-NULL;
result<-pFlowMerge(rituximab,cl=NULL,varNames=c("FSC.H","SSC.H"),K=1:8);
plot(result);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc}
