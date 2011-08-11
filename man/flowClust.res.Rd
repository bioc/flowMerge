\name{flowClust.res}
\alias{flowClust.res}
\alias{RituximabFlowClustFit}
\docType{data}
\title{ A flowClust model fitted to the rituximab data for 1:10 clusters.}
\description{
  The Rituximab data set accessible via \code{data(rituximab)} in the \code{flowClust} package fitted to a flowClust model containing from one to ten components. The results are in the object \code{flowClust.res}.
}
\usage{data(RituximabFlowClustFit)}
\format{
  The format is:
  \code{flowClust.res} is a \code{flowClustList}, where each element of the list is a \code{flowClust} model of the \code{rituximab} data, for \code{K=1} through \code{K=10} components, respectively. The structure of \code{flowClustList} and \code{flowClust} can be found in the corresponding documentation of the \code{flowClust} package. The format of the \code{rituximab} data is found in the documentation for that data set.
}
\details{
  The models have been precomputed for use in \code{flowMerge} examples to save computation time. \code{flowClust} was called on the rituximab data to generate these models with the following command: \code{flowClust.res<-flowClust(rituximab,K=1:10,B=1000,B.init=100,tol=1e-5,tol.init=1e-2,nu=4,randomStart=50,trans=1,nu.est=1)}.
}
\source{
   Gasparetto, M., Gentry, T., Sebti, S., O'Bryan, E., Nimmanapalli,
     R., Blaskovich, M. A., Bhalla, K., Rizzieri, D., Haaland, P.,
     Dunne, J. and Smith, C. (2004) Identification of compounds that
     enhance the anti-lymphoma activity of rituximab using flow
     cytometric high-content screening. \emph{J. Immunol. Methods} \bold{292},
     59-71.}
%\references{
  
%}
\examples{
#data(RituximabFlowClustFit)
#summary(flowClust.res);
}
\keyword{datasets}
