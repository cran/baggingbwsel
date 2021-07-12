#' Bagged CV bandwidth selector for Nadaraya-Watson estimator
#' 
#' @param x Covariate vector.
#' @param y Response vector.
#' @param r Positive integer. Size of the subsamples.
#' @param s Positive integer. Number of subsamples.
#' @param h0 Positive real number. Range over which to minimize, left bound.
#' @param h1 Positive real number. Range over which to minimize, right bound.
#' @param nb Positive integer. Number of bins to use in cross-validation.
#' @param ncores Positive integer. Number of cores with which to parallelize the computations.
#' 
#' @details
#' Bagged cross-validation bandwidth selector for the Nadaraya-Watson estimator.
#' 
#' @return Bagged CV bandwidth.
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(10^5)
#' y <- 2*x+rnorm(1e5,0,0.5)
#' bagreg(x, y, 1000, 10, 0.01, 1, 1000, 2)
#' 
#' @importFrom foreach %dopar%
#' @export
bagreg <- function(x,y,r,s,h0,h1,nb=r,ncores=parallel::detectCores())
{
  n <- length(x)
  sx.lst = list()
  sy.lst = list()
  for(i in 1:s)
  {
    idx = sample(1:n,r)
    sx.lst[[i]] = x[idx]
    sy.lst[[i]] = y[idx]
  }
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  paroutput <- foreach::foreach(i=1:s,.combine=c,.noexport=c("x","y")) %dopar%{
    subx = sx.lst[[i]]
    suby = sy.lst[[i]]
    return(h.select(subx,suby,lower=h0,upper=h1,method="cv",poly.index=0,nbins=nb))
  }
  parallel::stopCluster(cl)
  hmean <- mean(paroutput)*(r/n)^0.2
  return(hmean)
}