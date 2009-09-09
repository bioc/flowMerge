\name{flowObj-class}
\docType{class}
\alias{flowObj-class}
\alias{getData,flowObj-method}
\alias{merge,flowObj,missing-method}
\alias{plot,flowObj,missing-method}

\title{Class "flowObj" }
\description{A class inheriting from \code{flowClust} that groups the model and data in a single object.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("flowObj", ...). Has a convenience method \code{flowObj(flowClustObj, flowFrameObj)} for creating instances of the class.
}
	
}
\section{Slots}{
	 \describe{
    \item{\code{DATA}:}{An \code{"environment"} that holds a pointer to the \code{flowFrame} data in position \code{[[1]]}. }
    \item{\code{expName}:}{As described in the \code{flowClust} documentation }
    \item{\code{varNames}:}{As described in the \code{flowClust} documentation }
    \item{\code{K}:}{As described in the \code{flowClust} documentation }
    \item{\code{w}:}{As described in the \code{flowClust} documentation }
    \item{\code{mu}:}{As described in the \code{flowClust} documentation }
    \item{\code{sigma}:}{As described in the \code{flowClust} documentation }
    \item{\code{lambda}:}{As described in the \code{flowClust} documentation }
    \item{\code{nu}:}{As described in the \code{flowClust} documentation }
    \item{\code{z}:}{As described in the \code{flowClust} documentation }
    \item{\code{u}:}{As described in the \code{flowClust} documentation }
    \item{\code{label}:}{As described in the \code{flowClust} documentation }
    \item{\code{uncertainty}:}{As described in the \code{flowClust} documentation }
    \item{\code{ruleOutliers}:}{As described in the \code{flowClust} documentation }
    \item{\code{flagOutliers}:}{As described in the \code{flowClust} documentation }
    \item{\code{rm.min}:}{As described in the \code{flowClust} documentation }
    \item{\code{rm.max}:}{As described in the \code{flowClust} documentation }
    \item{\code{logLike}:}{As described in the \code{flowClust} documentation }
    \item{\code{BIC}:}{As described in the \code{flowClust} documentation }
    \item{\code{ICL}:}{As described in the \code{flowClust} documentation }
  }
}
\section{Extends}{
Class \code{"\linkS4class{flowClust}"}, directly.
}
\section{Methods}{
  \describe{
    \item{getData}{\code{signature(obj = "flowObj")}: Retreives the contents of the \code{DATA} environment }
    \item{merge}{\code{signature(x = "flowObj", y = "missing")}: the \code{flowMerge} algorithm is called via this function on objects of type \code{flowObj}. }
    \item{plot}{\code{signature(x = "flowObj", y = "missing")}: A simplified plotting method. Does not require specification of the data since it is contained in the \code{flowObj} object. Takes most of the same parameters as \code{plot.flowClust}, except the \code{data} parameter }
	 }
}
\references{ Finak G, Bashasharti A, Brinkmann R, Gottardo R. Merging Mixture Model Components for Improved Cell Population Identification in High Throughput Flow Cytometry Data; Advances in Bioinformatics (To Appear)}
\author{ Greg Finak \email{<greg.finak@ircm.qc.ca>}, Raphael Gottardo \email{<raphael.gottardo@ircm.qc.ca>} }
 
\seealso{
	\link{flowMerge-class}, \link{flowObj}
}

\keyword{classes}
