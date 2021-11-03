\name{mclustSBIClearnCoeff}
\alias{mclustSBIClearnCoeff}
\alias{EMclust}

\title{Learning Coefficient of Gaussian mixture models}

\description{
Computes the learning coefficient for a parameterized Gaussian mixture model with \eqn{G} mixture components, while assuming 
that the 'true' number of components is \eqn{G0}. Learning coefficients are together with maximum likelihood estimate and multiplicity the main determinants of the Singular BIC.
}

\usage{
mclustSBIClearnCoeff(modelName, d, G, G0, dir = 0,
          noise = FALSE,
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
  \item{G}{
    An integer specifying the number of mixture components
    (clusters) for which the sBIC is to be calculated.
  }
  \item{G0}{
    An integer specifying the 'true' number of of mixture components
    (clusters) such that \eqn{G0 \le G}. Note, that learning coefficients are computed for pairs of integers \eqn{(G,G0)},
    which corresponds to the following situation: model with \eqn{G} components is fitted, 
    while the actual true number of components is \eqn{G0})
  }
  \item{dir}{
    Determines Dirichlet prior parameter to \eqn{\phi = dir * (r/2 + 1) + r/2}, 
    where \eqn{r} is the number of parameters for each mixture component. It is a prior for 
    the distribution of the random vector of mixture weights.
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
Return an object of class \code{'numeric'}, which is the learning coefficient.
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
  \code{\link{nMclustParamsComp}},
  \code{\link{nMclustParamsShared}},
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
irisLearnCoeff = mclustSBIClearnCoeff(modelName = 'VEV', d = ncol(iris[,-5]), G = 3, G0 = 2, dir = .1, noise = FALSE, equalPro = FALSE)



}

\keyword{cluster}
% docclass is function
