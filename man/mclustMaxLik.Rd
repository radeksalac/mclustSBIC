\name{mclustMaxLik}
\alias{mclustMaxLik}
\alias{EMclust}
\alias{print.mclustMaxLik}

\title{Maximum log-likelihood of Gaussian mixture models}

\description{
Maximum log-likelihood of Gaussian mixture models computed by EM algorithm initialized by model-based hierarchical clustering.
}

\usage{
mclustMaxLik(data, G = NULL, modelNames = NULL, 
          restart=100, 
          initialization = list(hcPairs = NULL, 
                                subset = NULL, 
                                noise = NULL), 
          x = NULL, 
          \dots)
}
\arguments{
  \item{data}{
    A numeric vector, matrix, or data frame of observations. Categorical
    variables are not allowed. If a matrix or data frame, rows
    correspond to observations and columns correspond to variables. 
  }
  \item{G}{
    An integer vector specifying the numbers of mixture components
    (clusters) for which the sBIC is to be calculated. 
    The default is \code{G=1:9}, unless the argument \code{x} is specified, 
    in which case the default is taken from the values associated 
    with \code{x}. 
  }
  \item{modelNames}{
    A vector of character strings indicating the models to be fitted 
    in the EM phase of clustering. The help file for
    \code{\link{mclustModelNames}} describes the available models.
    The default is:
    \describe{
        \item{\code{c("E", "V")}}{for univariate data}
        \item{\code{mclust.options("emModelNames")}}{for multivariate data (n > d)}
        \item{\code{c("EII", "VII", "EEI", "EVI", "VEI", "VVI")}}{the spherical and diagonal models for multivariate data (n <= d)}
     }
    unless the argument \code{x} is specified, in which case
    the default is taken from the values associated with \code{x}. 
   }
   \item{restart}{
      Number of different random restarts of the EM algorithm.  This is
      invoked only when the default optimization used in the function
      \code{mclustBIC} returns \code{NA} for the maximum log-likelihood.
  }
  \item{initialization}{
   A list containing zero or more of the following components:
  \describe{
  \item{\code{hcPairs}}{
    A matrix of merge pairs for hierarchical clustering such as produced
    by function \code{\link{hc}}. \cr
    For multivariate data, the default is to compute a hierarchical 
    agglomerative clustering tree by applying function \code{\link{hc}} with 
    model specified by \code{mclust.options("hcModelName")}, and
    data transformation set by \code{mclust.options("hcUse")}.\cr
    All the input or a subset as indicated by the \code{subset} argument is 
    used for initial clustering.\cr
    The hierarchical clustering results are then used to start the EM
    algorithm from a given partition.\cr
    For univariate data, the default is to use quantiles to start the EM
    algorithm. However, hierarchical clustering could also be used by 
    calling \code{\link{hc}} with model specified as \code{"V"} or \code{"E"}.
  }
  \item{\code{subset}}{
    A logical or numeric vector specifying a subset of the data
    to be used in the initial hierarchical clustering phase.
    By default no subset is used unless the number of observations exceeds 
    the value specified by \code{mclust.options("subset")}. 
    The \code{subset} argument is ignored if \code{hcPairs} are provided.
    Note that to guarantee exact reproducibility of results a seed must be 
    specified (see \code{\link{set.seed}}).
  }
  \item{\code{noise}}{
    A logical or numeric vector indicating an initial guess as to
    which observations are noise in the data. If numeric the entries
    should correspond to row indexes of the data. If supplied, a noise
    term will be added to the model in the estimation.
  }
  }
  }
  \item{x}{
     An object of class \code{'mclustBIC'}. If supplied, \code{mclustBIC}
     will use the settings in \code{x} to produce another object of
     class \code{'mclustBIC'}, but with \code{G} and \code{modelNames}
     as specified in the arguments. Models that have already been computed
     in \code{x} are not recomputed. All arguments to \code{mclustBIC} 
     except \code{data}, \code{G} and \code{modelName} are
     ignored and their values are set as specified in the attributes of
     \code{x}. 
     Defaults for \code{G} and \code{modelNames} are taken from \code{x}.
  }
  \item{\dots}{
    Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}

\value{
Returns an object of class \code{'mclustMaxLik'} containing a matrix of maximum log-likelihood values for the respective number of mixture components (rows) and model (columns). 
Auxiliary information is returned as attributes.

The corresponding \code{print} method shows the matrix as it was described above with the correponding row and column names.
}

\references{
Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, \emph{The R Journal}, 8/1, pp. 289-317.

Drton M. and Plummer M. (2017) A Bayesian information criterion for singular models, \emph{Journal of the Royal Statistical Society}, Series B, 79, Part 2, pp. 323-380.

Fraley C. and Raftery A. E. (2002) Model-based clustering, discriminant analysis and density estimation, \emph{Journal of the American Statistical Association}, 97/458, pp. 611-631.

Fraley C., Raftery A. E., Murphy T. B. and Scrucca L. (2012) mclust Version 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification, and Density Estimation. \emph{Technical Report} No. 597, Department of Statistics, University of Washington.
}
\seealso{
  \code{\link{mclustSBIC}},
  \code{\link{mclustSBIClearnCoeff}},
  \code{\link{nMclustParamsComp}},
  \code{\link{nMclustParamsShared}},
  \code{\link{priorControl}}, 
  \code{\link{emControl}}, 
  \code{\link{mclustModel}},
  \code{\link{summary.mclustMaxLik}}
  \code{\link{summary.mclustSBIC}}, 
  \code{\link{hc}},
  \code{\link{me}},
  \code{\link{mclustModelNames}},
  \code{\link{mclust.options}}
}
\examples{
irisMaxLikValues = mclustMaxLik(iris[,-5], G=seq(from=1,to=9,by=2), 
                    modelNames=c("EII", "EEI", "EEE"))

}
\keyword{cluster}
% docclass is function
