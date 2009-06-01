\name{checkForRemoteErrors}
\alias{checkForRemoteErrors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Check output of snow clusters for errors }
\description{
 Overrides the snow \code{checkForRemoteErrors} function. Try errors are returned when cluster nodes produce errors, rather than completely aborting the computation. Not meant to be called by the user.
}
\usage{
checkForRemoteErrors(val)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{val}{ The result returned from an individual cluster node. }
}
\details{
  This function is meant to be called internally, but must be exported so that it can hide the native \code{checkForRemoteErrors} function in the snow package.
}
\value{
  The result from the snow cluster node, or an object of type \code{try-error} if there was an error.
  }
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data (Submitted)}
\author{ Greg Finak \email{<greg.finak@ircm.qc.ca>}}
\seealso{ \link[snow]{checkForRemoteErrors} }

\keyword{ misc }
