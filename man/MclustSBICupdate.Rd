\name{MclustSBICupdate}
\alias{MclustSBICupdate}
\alias{EMclust}

\title{SBIC Update Function for the BIC Based Mclust}

\description{
Updates the \code{'Mclust'} output such that the optimal model is selected in accordance to the Singular BIC (sBIC) criterion. The class remains unchanged so that other functions and methods from  package \pkg{mclust} can be applied. Note, that some details such as plot titles or axis labels might not be adjusted.}

\usage{
MclustSBICupdate(object,\dots)
}
\arguments{
  \item{object}{
    An object of class \code{'Mclust'} resulting from a call to \code{\link{Mclust}}.
  }
}

\value{
Return an object of class \code{'Mclust'} providing the optimal (according to sBIC) mixture model estimation.
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
  \code{\link{nMclustParamsShared}},
  \code{\link{priorControl}}, 
  \code{\link{emControl}}, 
  \code{\link{mclustModel}}, 
  \code{\link{summary.mclustSBIC}}, 
  \code{\link{hc}},
  \code{\link{me}},
  \code{\link{mclustModelNames}},
  \code{\link{mclust.options}}
}
\examples{
X <- iris[,-length(iris)]
b <- mclustBIC(X)
modelBIC <- Mclust(X, x = b)
plot(modelBIC, what="BIC")

modelSBIC <- MclustSBICupdate(modelBIC)
# now plot the sBIC values
# note, using the class 'Mclust' implies the y-axis label 'BIC'
# however, since MclustSBICupdate was applied, the plotted values truly are the sBIC
plot(modelSBIC, what="BIC")

modDR = MclustDR(modelSBIC)
plot(modDR,what = "classification") 
}
\keyword{cluster}
% docclass is function
