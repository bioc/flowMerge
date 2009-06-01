\name{flowClust-class}
\docType{class}
\alias{flowClust-class}

\title{Class "flowClust" }
\description{The basic class used in the \pkg{flowClust} package}
\section{Slots}{
	 \describe{
	\item{\code{BIC}:}{The Bayesian Information Criterion for the fitted mixture model. }
    \item{\code{ICL}:}{The Integrated Completed Likelihood for the fitted mixture model. }
    \item{\code{expName}:}{A \code{"character"} vector with the name of the experiment. }
    \item{\code{flagOutliers}:}{ A logical vector of size N, showing whether each data point is called an outlier or not based on the rule defined by 'level'/'u.cutoff' and 'z.cutoff'. }
    \item{\code{varNames}:}{A \code{"character"} vector containing the variable names of the data channels. Usually named something like \code{"FL1 Log"}, \code{"FL2 Log"}, \code{"FS Lin"}, code{"SS Lin"}, rather than the names of the immunohistochemical markers. }
    \item{\code{K}:}{The number of components in the model. }
    \item{\code{w}:}{A \code{"vector"} containing the weights or proportions associated with each mixture component. }
    \item{\code{mu}:}{A \code{"matrix"} of K x p dimensions containing the means of the mixture components in the model. }
    \item{\code{sigma}:}{An \code{"array"} of K x p x p containing the sigma parameters for the mixture components in the model. }
    \item{\code{lambda}:}{A \code{"numeric"} with the estimated parameter \code{lambda} used in the model. }
    \item{\code{logLike}:}{The log-likelihood of the fitted mixture model.}
    \item{\code{nu}:}{A \code{"numeric"}. The degrees of freedom used in the model. }
    \item{\code{rm.min}:}{ Number of points filtered from below. }
    \item{\code{rm.max}:}{Number of points filtered from above. }
    \item{\code{z}:}{A \code{"matrix"} of size N x K containing the posterior probabilities of cluster memberships. The probabilities in each row sum up to one. }
    \item{\code{u}:}{A \code{"matrix"} of size N x K, containing the "weights" (the contribution for computing cluster mean and covariance matrix) of each data point in each cluster. Since this quantity decreases monotonically with the Mahalanobis distance, it can also be interpreted as the level of "outlyingness" of a data point.  Note that, when 'nu=Inf', this slot is used to store the Mahalanobis distances instead.
 }
    \item{\code{label}:}{A \code{"vector"} of size N, showing the cluster membership according to the initial partition (i.e., hierarchical clustering if 'randomStart' is 'FALSE').  Filtered observations will be labelled as 'NA'.  Unassigned observations (which may occur since only 1500 observations at maximum are taken for hierarchical clustering) will be labelled as 0. }
    \item{\code{uncertainty}:}{A \code{"vector"} of size N, containing the uncertainty about the cluster assignment.  Uncertainty is defined as 1 minus the posterior probability that a data point belongs to the cluster to which it is assigned.
 }
    \item{\code{ruleOutliers}:}{ A numeric vector of size 3, storing the rule used to call outliers.  The first element is 0 if the criterion is set by the 'level' argument, or 1 if it is set by 'u.cutoff'.  The second element copies the content of either the 'level' or 'u.cutoff' argument.  The third element copies the content of the 'z.cutoff' argument.  For instance, if points are called outliers when they lie outside the 90% quantile region or have assignment probabilities less than 0.5, then 'ruleOutliers' is 'c(0, 0.9, 0.5)'.  If points are called outliers only if their "weights" in the assigned clusters are less than 0.5 regardless of the assignment probabilities, then 'ruleOutliers' becomes 'c(1, 0.5, 0)'. }
    
        
  }
}}
\section{Methods}{
See the \pkg{flowClust} package.
}
\references{ Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytomet
ry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332
. }
\author{ Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}> }
\note{Please see the \code{flowClust} package documentation for a detailed description of the slots in a \code{flowClust} object. This class is the parent class for \code{flowMerge} and \code{flowObj} classes used in the \code{flowMerge} package.
}

\seealso{
	\link[flowClust]{flowClust}
}
\keyword{classes}
