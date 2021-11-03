#' @export
mclustSBIC <- function (data, G = NULL, modelNames = NULL, restart=10, initialization = list(hcPairs = NULL, subset = NULL, noise = NULL), x = NULL, dir = 0.1, ...)
{
  if (!is.null(x)){
    if(class(x) != 'mclustBIC'){
      stop("If provided, argument x must be an object of class 'mclustBIC'")
    }
  }
  call <- match.call()
  data <- data.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name("mclustMaxLik")
  mc[[2]] <- data
  ML <- eval(mc, parent.frame())*2
  Ml = ML  ##Radek: not used!
  class(ML) <- "mclustMaxLik"
  G <- attr(ML, "G")
  
  modelNames <- attr(ML, "modelNames")
  
  sBIC <- matrix(NA, nrow = length(G), ncol = length(modelNames))
  
  mostattributes(sBIC) <- attributes(ML)
  
  #print(ML)
  
  if(any(is.na(ML))){ "some models have NA maximum likelihood"}
  ## need a check whether G has at least one element?
  for (m in 1:ncol(sBIC)) {
    if(!is.na(ML[1,m])){
      sBIC[1,m] <-
        ML[1,m] - log(n)*2*
        mclustSBIClearnCoeff(modelNames[m], d, G[1], G[1], dir = dir,
                             #noise = exists("Vinv"))
                             noise = FALSE)
    }
    ## HELP:  How to deal with the noise option here?  noise = exists...) etc. did not work, again below
  }
  logsum <- function(lx){
    ## compute log(sum x[i]) from lx=log(x)
    c <- max(lx)
    return(c+log(sum(exp(lx-c))))
  }
  logdiff <- function(lx,ly){
    ## compute log(x  - y) from lx=log(x) and ly=log(y)
    ## pay attention to sign
    c <- max(lx,ly)
    if(lx>ly){
      return(list(abs=c+log(exp(lx-c)-exp(ly-c)), sign=1))
    }
    if(lx<ly){
      return(list(abs=c+log(exp(ly-c)-exp(lx-c)), sign=-1))
    }
    if(lx==ly){
      return(list(abs=-Inf, sign=0))
    }
  }
  
  if(nrow(sBIC)>1){
    for (i in 2:nrow(sBIC)) {
      for (m in 1:ncol(sBIC)) {
        if (is.na(ML[i,m])) 
          (next)()
        logLi <- rep(ML[i,m],i)
        for(j in 1:i){
          ## HELP:  How to deal with the noise option here?  noise = exists("Vinv")) etc. did not work, again below
          logLi[j] <- logLi[j] - 2*mclustSBIClearnCoeff(modelNames[m], d, G[i], G[j], dir = dir, noise = exists("Vinv"))*log(n)
        }
        #browser()
        #print(sBIC)
        logbi <- logdiff( logsum(sBIC[1:(i-1),m]), logLi[i] )
        
        logci <- logsum( logLi[1:(i-1)] + sBIC[1:(i-1),m] )
        ## Numerically stable formula for the root is either standard one or 2c/(b+sqrt(b^2+4c))
        if(logbi$sign==0){
          sBIC[i,m] <- 1/2*logci
        }
        if(logbi$sign<0){
          ## now logbi is already absolute value of bi
          sBIC[i,m] <-
            logsum( c( 1/2*logsum(c(2*logbi$abs,log(4)+logci) ),
                       logbi$abs) ) - log(2)
        }
        if(logbi$sign>0){
          sBIC[i,m] <- log(2)+logci-
            logsum(c(logbi$abs,1/2*logsum(c(2*logbi$abs,log(4)+logci))))
        }
      }
    }
  }
  class(sBIC) <- "mclustSBIC" # "mclustBIC"
  attr(sBIC, "criterion") <- "sBIC"
  return(sBIC)
}

#' @export
print.mclustSBIC <- function (x, pick = 3, ...) 
{
  subset <- !is.null(attr(x, "subset"))
  oldClass(x) <- attr(x, "args") <- NULL
  attr(x, "criterion") <- NULL
  attr(x, "control") <- attr(x, "initialization") <- NULL
  attr(x, "oneD") <- attr(x, "warn") <- attr(x, "Vinv") <- NULL
  attr(x, "prior") <- attr(x, "G") <- attr(x, "modelNames") <- NULL
  ret <- attr(x, "returnCodes") == -3
  n <- attr(x, "n")
  d <- attr(x, "d")
  attr(x, "returnCodes") <- attr(x, "n") <- attr(x, "d") <- NULL
  
  oldClass(x) <- attr(x, "args") <- attr(x, "criterion") <- NULL 
  cat("sBIC:\n")
  print(x, ...)
  cat("\n")
  cat("Top", pick, "models based on the sBIC criterion:\n")
  print(pickBIC(x, pick), ...)
  invisible()
}
#' @export
plot.mclustSBIC <- function(x, ylab = "sBIC", ...) 
{
  plot.mclustBIC(x, ylab = ylab, ...)  
}



#' @export
mclustMaxLik <- function (data, G = NULL, modelNames = NULL, restart=10, initialization = list(hcPairs = NULL, subset = NULL, noise = NULL), x = NULL, ...)
  ## this is max log-likelihood
  ## from mclustBIC.Rd: All arguments to mclustBIC except data, G and modelName 
  # are ignored and their values are set as specified in the attributes of x. 
  # Defaults for G and modelNames are taken from x.
{
  call <- match.call()
  data <- data.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name("mclustBIC")
  mc[[2]] <- data
  BIC <- eval(mc, parent.frame()) ##Radek: appllies mclustBIC, not efficient, sollte eigentlich in else sein
  ## print(BIC)
  class(BIC) <- "mclustBIC"
  G <- attr(BIC, "G")
  
  modelNames <- attr(BIC, "modelNames")
  MaxLik <- matrix(NA, nrow = length(G), ncol = length(modelNames))
  mostattributes(MaxLik) <- attributes(BIC)
  # if(!is.null(x)){
  #   for (i in 1:nrow(MaxLik)){
  #     for(j in 1:ncol(MaxLik)){
  #       MaxLik[i,j]<- x[i, j] + log(n)*nMclustParams(modelNames[j], d, G[i], noise=!is.null(initialization$noise) )
  #       if (is.na(MaxLik[i,j])){
  #         ##Radek: if NA: since they do it just once we want to try it more times
  #         MaxLik[i, j]<- max_loglik(data,n,i,modelname = modelNames[j],restart = restart)
  #       }
  #     }
  #   }
  # }
  # else{
  for (i in 1:nrow(MaxLik)) {
    for (j in 1:ncol(MaxLik)) {
      if (is.na(BIC[i, j])){
        MaxLik[i, j]<- max_loglik(data,n,i,modelname = modelNames[j],restart = restart)
      }
      else{
        # Sumry <- summary(BIC, data, G = G[i], modelNames = modelNames[j])
        # MaxLik[i, j] <- maxlik.Mclust(Sumry)
        MaxLik[i,j]<- BIC[i, j] + log(n)*nMclustParams(modelNames[j], d, G[i], noise=!is.null(initialization$noise) )
      }
    }
  }
  # }
  class(MaxLik) <- "mclustMaxLik"
  attr(MaxLik, "criterion") <- "MaxLik"
  return(MaxLik/2)
}


max_loglik<- function(data,n,i,modelname,restart){
  model.loglik=-Inf
  cluster.save<-matrix(rep(0,n*i),n,i)
  
  while(model.loglik==(-Inf)){
    for (j in 1:restart){
      cluster.density=NULL
      for (k in 1:n){
        gene.cluster=c(0,sort(runif(i-1)))
        gene.cluster=c(gene.cluster[-1]-gene.cluster[-i],1-gene.cluster[i])
        cluster.density=rbind(cluster.density,gene.cluster)
      }
      emEst <- me(modelName=modelname, data=data, z=cluster.density )
      
      if(!is.na(emEst$loglik)){
        if(emEst$loglik>model.loglik){
          model.loglik<- emEst$loglik
        }
      }
    }
  }
  
  return(2*model.loglik)
} 







#' @export
summary.mclustMaxLik <- function(object, G, modelNames, ...)
{
  if(!missing(G)) 
    object <- object[rownames(object) %in% G,,drop=FALSE]
  if(!missing(modelNames)) 
    object <- object[,colnames(object) %in% modelNames,drop=FALSE]
  structure(pickBIC(object, ...),
            class = "summary.mclustMaxLik")
}


maxlik.Mclust <- function(object, ...)
{
  attr(object$bic, "Vinv")
  n <- object$n
  object$bic + log(n)*nMclustParams(object$modelName, object$d, object$G, noise = !is.null(object$Vinv))
}
#' @export
print.mclustMaxLik <- function (x, pick = 3, ...) 
{
  subset <- !is.null(attr(x, "subset"))
  oldClass(x) <- attr(x, "args") <- NULL
  attr(x, "criterion") <- NULL
  attr(x, "control") <- attr(x, "initialization") <- NULL
  attr(x, "oneD") <- attr(x, "warn") <- attr(x, "Vinv") <- NULL
  attr(x, "prior") <- attr(x, "G") <- attr(x, "modelNames") <- NULL
  ret <- attr(x, "returnCodes") == -3
  n <- attr(x, "n")
  d <- attr(x, "d")
  attr(x, "returnCodes") <- attr(x, "n") <- attr(x, "d") <- NULL
  
  oldClass(x) <- attr(x, "args") <- attr(x, "criterion") <- NULL 
  cat("Maximum Log-Likelihood:\n")
  print(x, ...)
  invisible()
}


#' @export
mclustSBIClearnCoeff <- function (modelName, d, G, G0, dir = 0, noise = FALSE, ...)
  ## d : dimension
  ## G : number of components  (i in sBIC paper)
  ## G0 : number of 'true' components  (j in sBIC paper)
  ## dir : determines Dirichlet prior parameter to
  ## phi = dir*(r/2+1)+r/2 where r is the
  ## number of parameters for each component
  ## dir = 1 recovers the BIC solution
  ## output is the learning coefficient
{
  if(G == 0){
    browser()
    if (!noise) 
      stop("undefined model")
    twolambda <- 1
  }
  else{
    if(G0>G){
      G0 = G ## submodels only
    }
    r = nMclustParamsComp(modelName, d = d)
    phi = dir*(r/2+1)+r/2
    twolambda = G0*r + G0 - 1 + (G - G0)*phi +
      nMclustParamsShared(modelName, d = d, noise = noise)
  }
  return(twolambda/2)
}

#' @export
nMclustParamsComp <- function (modelName, d, ...) 
  ## d : dimension
  ## G : number of components
  ## output : number of parameters specific to each component
  ## this is not counting the mixture parameters
{
  modelName <- switch(EXPR = modelName, X = "E", XII = "EII", 
                      XXI = "EEI", XXX = "EEE", modelName)
  checkModelName(modelName)
  nparams <- nVarParamsComp(modelName, d = d) + d
  ## d mean parameters added
  return(nparams)
}
#' @export
nMclustParamsShared <- function (modelName, d, noise = FALSE, ...) 
  ## d : dimension
  ## G : number of components
  ## output : number of parameters shared across components
  ## this is not counting the mixture parameters
{
  modelName <- switch(EXPR = modelName, X = "E", XII = "EII", 
                      XXI = "EEI", XXX = "EEE", modelName)
  checkModelName(modelName)
  nparams <- nVarParamsShared(modelName, d = d) 
  if (noise) 
    nparams <- nparams + 2
  return(nparams)
}

nVarParamsComp <- function (modelName, d, ...)
  ## d : dimension
  ## output : number of parameters for each component distribution
{
  modelName <- switch(EXPR = modelName, X = "E", XII = "EII", 
                      XXI = "EEI", XXX = "EEE", modelName)
  switch(EXPR = modelName, E = 0, V = 1, EII = 0, VII = 1, 
         EEI = 0, VEI = 1, EVI = (d - 1), VVI = d,
         EEE = 0, EVE = (d - 1), VEE = 1, VVE = d,
         EEV = d * (d - 1)/2, VEV = 1 + d * (d - 1)/2,
         EVV = - 1 + d * (d + 1)/2, VVV = d * (d + 1)/2,
         stop("invalid model name"))
}

nVarParamsShared <- function (modelName, d, ...)
  ## d : dimension
  ## output : number of parameters shared across components
{
  modelName <- switch(EXPR = modelName, X = "E", XII = "EII", 
                      XXI = "EEI", XXX = "EEE", modelName)
  switch(EXPR = modelName, E = 1, V = 0, EII = 1, VII = 0, 
         EEI = d, VEI = (d - 1), EVI = 1, VVI = 0,
         EEE = d * (d + 1)/2, EVE = 1 + d * (d - 1)/2,
         VEE = (d - 1) + d * (d - 1)/2, VVE = d * (d - 1)/2,
         EEV = 1 + (d - 1), VEV = (d - 1),
         EVV = 1, VVV = 0, stop("invalid model name"))
}


#################### RADEK:
#' @export
MclustSBICupdate <- function (object, ...)
  ## object of class 'Mclust'
{
  SBIC = mclustSBIC(X,x=modBIC$BIC)
  class(SBIC) = "mclustBIC"
  modSBIC = Mclust(X, x = SBIC)
  return(modSBIC)
}


#################### RADEK: sbic.R (from icl.R)
## plot.mclustSBIC and print.mclustSBIC are already above
#' @export
summary.mclustSBIC <- function(object, G, modelNames, ...)
{
  if(!missing(G)) 
    object <- object[rownames(object) %in% G,,drop=FALSE]
  if(!missing(modelNames)) 
    object <- object[,colnames(object) %in% modelNames,drop=FALSE]
  structure(pickBIC(object, ...),
            class = "summary.mclustSBIC")
}
#' @export
print.summary.mclustSBIC <- function(x, digits = getOption("digits"), ...)
{
  cat("Best sBIC values:\n")
  x <- drop(as.matrix(x))
  x <- rbind(sBIC = x)
  print(x, digits = digits)
  invisible()
}


