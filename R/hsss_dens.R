#' Generalized bagging CV bandwidth selector for Parzen-Rosenblatt estimator
#' 
#' @param x Vector. Sample.
#' @param r Positive integer. Size of the subsamples.
#' @param s Positive integer. Number of subsamples.
#' @param nb Positive integer. Number of bins.
#' @param h0 Positive real number. Range over which to minimize, left bound.
#' @param h1 Positive real number. Range over which to minimize, right bound.
#' @param ncores Positive integer. Number of cores with which to parallelize the computations.
#' 
#' @details
#' Generalized bagging cross-validation bandwidth selector for the Parzen-Rosenblatt estimator.
#' 
#' @return Bagged CV bandwidth.
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(10^5)
#' hsss_dens(x, 5000, 100, 1000, 0.001, 1, 2)
#' 
#' @export
hsss_dens <- function(x,r,s,nb=r,h0,h1,ncores=parallel::detectCores())
{
  n <- length(x)
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  lr <- length(r)
  slr <- sort(r)
  rvec <- rep(slr,each=s)
  sx.lst = list()
  for(i in 1:(lr*s))
  {
    sx.lst[[i]] <- sample(x,rvec[i],FALSE)
  }
  
  j<-NULL
  paroutput <- foreach::foreach(j=1:(lr*s),.combine=cbind,.noexport="x") %dopar% {
    sx <- sx.lst[[j]]
    return(list(bw=stats::bw.ucv(sx,nb=nb,lower=h0,upper=h1),id=rvec[j]))
  }

  parallel::stopCluster(cl)
  hmeans <- sapply(slr,function(ri)mean(unlist(paroutput[1,which(paroutput[2,]==ri)])))
  
  linmod2 <- stats::lm(log(hmeans)~log(slr))$coefficients
  c0_2 <- exp(linmod2[1])
  p_2 <- linmod2[2]
  h_lm_2 <- c0_2*n^p_2
  
  
  return(hmean=h_lm_2)
  
}