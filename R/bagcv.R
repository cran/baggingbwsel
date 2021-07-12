#' Bagged CV bandwidth selector for Parzen-Rosenblatt estimator
#' 
#' @param x Vector. Sample.
#' @param r Positive integer. Size of the subsamples.
#' @param s Positive integer. Number of subsamples.
#' @param h0 Positive real number. Range over which to minimize, left bound.
#' @param h1 Positive real number. Range over which to minimize, right bound.
#' @param nb Positive integer. Number of bins.
#' @param ncores Positive integer. Number of cores with which to parallelize the computations.
#' 
#' @details
#' Bagged cross-validation bandwidth selector for the Parzen-Rosenblatt estimator.
#' 
#' @return Bagged CV bandwidth.
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(10^6)
#' bagcv(x, 5000, 100, 0.01, 1, 1000, 2)
#' 
#' @importFrom foreach %dopar%
#' @export
bagcv <-
function(x,r,s,h0,h1,nb=r,ncores=parallel::detectCores())
{
  n <- length(x)
  sx.lst = list()
  for(i in 1:s)
  {
    sx.lst[[i]] = sample(x,r)
  }
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  paroutput <- foreach::foreach(i=1:s,.combine=c,.noexport="x") %dopar%{
    subx = sx.lst[[i]]
    return(stats::bw.ucv(subx,nb,h0,h1))
  }
  parallel::stopCluster(cl)
  hmean <- mean(paroutput)*(r/n)^0.2
  return(hmean)
}
