#' Second order bagging CV bandwidth selector for Parzen-Rosenblatt estimator
#' 
#' @param x Vector. Sample.
#' @param r Vector. The two subsample sizes.
#' @param s Positive integer. Number of subsamples.
#' @param h0 Positive real number. Range over which to minimize, left bound.
#' @param h1 Positive real number. Range over which to minimize, right bound.
#' @param nb Positive integer. Number of bins.
#' @param ncores Positive integer. Number of cores with which to parallelize the computations.
#' 
#' @details
#' Second order bagging cross-validation bandwidth selector for the Parzen-Rosenblatt estimator.
#' 
#' @return Second order bagging CV bandwidth.
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(10^5)
#' tss_dens(x, 5000, 10, 0.01, 1, 1000, 2)
#' 
#' @export
tss_dens <- function(x,r,s,h0,h1,nb=1000,ncores=1)
{
  n <- length(x)
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  r1 <- max(r)
  r2 <- min(r)
  rvec <- rep(c(r1,r2),each=s)
  sx.lst = list()
  for(i in 1:(2*s))
  {
    sx.lst[[i]] <- sample(x,rvec[i],FALSE)
  }
  i<-NULL
  paroutput <- foreach::foreach(i=1:(2*s),.combine=cbind,.noexport="x") %dopar% {
    sx <- sx.lst[[i]]
    return(list(bw=stats::bw.ucv(sx,nb=nb,lower=h0,upper=h1),id=rvec[i]))
  }
  parallel::stopCluster(cl)
  hr1 <- mean(unlist(paroutput[1,which(paroutput[2,]==r1)]))
  hr2 <- mean(unlist(paroutput[1,which(paroutput[2,]==r2)]))
  c0 <- hr1*r1^0.2
  c1 <- (hr2-c0*r2^-0.2)*r2^0.6
  return(c0*n^-0.2+c1*n^-0.6)
}