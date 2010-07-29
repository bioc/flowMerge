\name{fitPiecewiseLinreg}
\alias{fitPiecewiseLinreg}
\title{ Fit Piecewise Linear Regression for a list of flowMerge Objects}
\description{
 Fits a two--component piecewise linear regression to the entropy vs number of clusters for a list of merged cluster solutions.
 }
 \usage{
fitPiecewiseLinreg(x, plot=FALSE, normalized=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A "list" of \code{flowMerge} objects for \code{1} through \code{K} clusters derived from a single max BIC \code{flowObj} or \code{flowClust} object.  }
  \item{plot}{ A \code{logical} indicating whether to plot the fit or not. Default is FALSE. }
  \item{normalized}{ A \code{logical} indicating whether the merging should be done using the normalized or unnormalized entropy}
  \item{\dots}{ Additional arguments not currently used. }
}
\details{
  An S4 method that takes a list of \code{flowMerge} objects output by the \code{merge} method, extracts the entropy and fits a piecwise linear regression to the entropy vs number of clusters in order to find the postion of the changepoint. The location of the changepoint corresponds to the optimal merged cluster solution. The piecewise linear regression now is fitted to the entropy vs cumulative sum of merged observations at each number of clusters. This normalizes the change in entropy for the number of data points as described in Baudry et al.
}
\value{
  An \code{integer} value corresponding to the position of the changepoint.
  }
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data; Advances in Bioinformatics (To Appear)}\author{ Greg Finak \email{<greg.finak@ircm.qc.ca>} }
\seealso{ \link{merge}}
\examples{
data(rituximab)
data(RituximabFlowClustFit)
o<-flowObj(flowClust.res[[which.max(BIC(flowClust.res))]],rituximab);
m<-merge(o)
i<-fitPiecewiseLinreg(m);
}
\keyword{misc}

