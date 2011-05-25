\name{ptree} 
\alias{ptree} 
\title{Generate a Function to Plot The Merging Tree}
\description{ 
This function generates and returns a new function which can be used to plot the merging tree for a flowMerge model, with nodes highlighted based on the expression of different parameters for each cell population.
} 
\usage{ 
ptree(x,y)
} 
\arguments{ 
  \item{x}{A character string of the name of the variable holding the list of merged models returned from flowMerge} 
  \item{y}{The index of the best fitting merged model in that list} 
} 

\details{
	ptree will generate a function that will plot the merging tree from a flowMerge model. Nodes will be colored by the intensity of staining of that population in a given dimension. Calling \code{f<-ptree("model.name",fitPiecewiseLinreg(model.name))} will assign the function to \code{f}. Calling \code{f(3)} will plot the merging tree with nodes hightlighted according to parameter 3, presuming that there are that many parameters in the model.
}
\value{ 
	Returns a function
}
\section{Side Effects}{
  A plot will be drawn on the current device.
}
\seealso{ 
  \code{\link{merge}}. 
}
\author{
	Greg Finak <gfinak@fhcrc.org>
}
