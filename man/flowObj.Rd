\name{flowObj}
\alias{flowObj}

\title{ Create a flowObj object from a flowClust and flowFrame object}
\description{
  Convenience method that creates a \code{flowObj} object from a \code{flowClust} and \code{flowFrame} object, so as to group the model and data together. Useful for high-throughput analysis where one may want to access the data to compute other statistics.
}
\usage{
flowObj(flowC = NULL, flowF = NULL)
}

\arguments{
  \item{flowC}{ A \code{flowClust} object representing the model fit }
  \item{flowF}{ A \code{flowFrame} object on which the flowClust model is based. }
}
\details{
  Calls the \code{new("flowObj",..)} constructor.
}
\value{
    An object of class \code{flowObj-class}
}
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data; Advances in Bioinformatics (To Appear)}
\author{ Greg Finak \email{<greg.finak@ircm.qc.ca>}, Raphael Gottardo \email{<raphael.gottardo@ircm.qc.ca>} }
\seealso{ \link{flowObj-class} }
\examples{
data(rituximab)
data(RituximabFlowClustFit)
o<-flowObj(flowClust.res[[which.max(flowMerge:::BIC(flowClust.res))]],rituximab);
m<-merge(o);

}
\keyword{misc}

