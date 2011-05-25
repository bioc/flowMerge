\name{flowMerge-class}
\docType{class}
\alias{flowMerge-class}
\alias{getData,flowMerge-method}
\alias{plot,flowMerge,missing-method}


\title{Class "flowMerge"}
\description{A class to represent flowMerge objects}
\section{Objects from the Class}{
The object unites the flowMerge model output and the data being modeled and contains additional slots for various characteristics of a merged cluster solution, including the entropy of clustering.
}
\section{Slots}{
	 \describe{
    \item{\code{merged}:}{The number of observations merged at the current step in the algorithm.}
	\item{\code{mtree}:}{A tree--structured graph representing the order of merged components in the model. Inspired by SPADE. (Bendall et al.)}
    \item{\code{entropy}:}{The entropy of clustering of the current solution.  }
    \item{\code{DATA}:}{An environment whose first element contains the flowFrame with the data modeled by this flowMerge object  }
    \item{\code{expName}:}{See the \code{flowClust} package for details  }
    \item{\code{varNames}:}{See the \code{flowClust} package for details  }
    \item{\code{K}:}{The number of clusters in the merged solution. See the \code{flowClust} package for details }
    \item{\code{w}:}{The proportions for each component in the merged solution. See the \code{flowClust} package for details}
    \item{\code{mu}:}{The means of the components in the merged solution. See the \code{flowClust} package for details }
    \item{\code{sigma}:}{The covraiances of the components in the merged solution. See the \code{flowClust} package for details }
    \item{\code{lambda}:}{See the \code{flowClust} package for details }
    \item{\code{nu}:}{See the \code{flowClust} package for details }
    \item{\code{z}:}{See the \code{flowClust} package for details }
    \item{\code{u}:}{The uncertainties for each data point.}
    \item{\code{label}:}{See the \code{flowClust} package for details }
    \item{\code{uncertainty}:}{See the \code{flowClust} package for details }
    \item{\code{ruleOutliers}:}{See the \code{flowClust} package for details }
    \item{\code{flagOutliers}:}{See the \code{flowClust} package for details }
    \item{\code{rm.min}:}{See the \code{flowClust} package for details }
    \item{\code{rm.max}:}{See the \code{flowClust} package for details }
    \item{\code{logLike}:}{See the \code{flowClust} package for details }
    \item{\code{BIC}:}{See the \code{flowClust} package for details }
    \item{\code{ICL}:}{See the \code{flowClust} package for details }
  }
}
\section{Extends}{
Class \code{"\linkS4class{flowObj}"}, directly.
Class \code{"\linkS4class{flowClust}"}, by class "flowObj", distance 2.
}
\section{Methods}{
  \describe{
    \item{getData}{\code{signature(obj = "flowMerge")}: Retrieves the \code{flowFrame} in the \code{DATA} environment slot. }
    \item{plot}{\code{signature(x = "flowMerge", y = "missing")}: Plots the clusters in this object. }
     \item{summary}{\code{signature(x="flowMerge")}: Prints a summary of the object.}
     \item{show}{\code{signature(x="flowMerge")}: Prints information about the object.}
	 }
}
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data (Submitted)}
\author{ Greg Finak \email{<greg.finak@ircm.qc.ca>} }


\seealso{\link{flowObj-class}
}
\keyword{classes}
