\name{nMclustParamsShared}
\alias{nMclustParamsShared}
\alias{EMclust}

\title{Number of Parameters Shared Across Gaussian Mixture Components}

\description{
  Number of parameters shared across components of a given parameterized Gaussian 
  mixture model, that is not counting the mixture parameters
}

\usage{
  nMclustParamsShared <- function (modelName, d, noise = FALSE,
                                 \dots)
}
\arguments{
  \item{modelName}{
    A character string indicating the model to be fitted 
    in the EM phase of clustering. The help file for
    \code{\link{mclustModelNames}} describes the available models.
  }
  \item{d}{
    Number of variables in the data.
  }
  \item{noise}{
    A logical or numeric vector indicating an initial guess as to
    which observations are noise in the data. If numeric the entries
    should correspond to row indexes of the data. If supplied, a noise
    term will be added to the model in the estimation.
  }
  \item{\dots}{
    Catches unused arguments in indirect or list calls via \code{do.call}.
  }
}

\value{
  Return an object of class \code{'numeric'}, which is the number of parameters shared across components excluding the mixture parameters.
  
}

\references{
  Scrucca L., Fop M., Murphy T. B. and Raftery A. E. (2016) mclust 5: clustering, classification and density estimation using Gaussian finite mixture models, \emph{The R Journal}, 8/1, pp. 289-317.
  
  Drton M. and Plummer M. (2017) A Bayesian information criterion for singular models, \emph{Journal of the Royal Statistical Society}, Series B, 79, Part 2, pp. 323-380.
  
  Fraley C. and Raftery A. E. (2002) Model-based clustering, discriminant analysis and density estimation, \emph{Journal of the American Statistical Association}, 97/458, pp. 611-631.
  
  Fraley C., Raftery A. E., Murphy T. B. and Scrucca L. (2012) mclust Version 4 for R: Normal Mixture Modeling for Model-Based Clustering, Classification, and Density Estimation. \emph{Technical Report} No. 597, Department of Statistics, University of Washington.
}
\seealso{
  \code{\link{mclustSBIC}},
  \code{\link{mclustMaxLik}},
  \code{\link{mclustSBIClearnCoeff}},
  \code{\link{nMclustParamsComp}},
  \code{\link{priorControl}}, 
  \code{\link{emControl}}, 
  \code{\link{mclustModel}}, 
  \code{\link{summary.mclustBIC}}, 
  \code{\link{hc}},
  \code{\link{me}},
  \code{\link{mclustModelNames}},
  \code{\link{mclust.options}}
}
\examples{
  irisParamsShared = nMclustParamsShared(modelName = 'VEV', d = ncol(iris[,-5]))
}
    
  \keyword{cluster}
    % docclass is function
    