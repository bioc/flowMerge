\name{map}
\alias{map}
\title{
 Map matrix of probabilities to class assignments.
}
\description{
Traverse the rows of a matrix of probabilities of size n x k, where the n rows are samples, and the k columns are the probability of assignment of the sample to each of k classes. The most probable class assignment is selected for each row and a vector of classes is returned.
}
\usage{
map(z, ...)
}
\arguments{
  \item{z}{
 A matrix of probabilities.
}
  \item{\dots}{
Additional arguments, not currently used.
}
}

\value{
A vector of class assignments of lenght \code{n}.
}
\author{
Greg Finak <greg.finak@ircm.qc.ca>, Raphael Gottardo <raphael.gottardo@ircm.qc.ca>
}


\examples{
z<-t(apply(t(replicate(100,rgamma(5,0.1,1))),1,function(x)x/sum(x)));
map(z);
}
\keyword{ misc }
