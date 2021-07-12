#' Bagging bootstrap bandwidth selector for Parzen-Rosenblatt estimator
#' 
#' @param x Vector. Sample.
#' @param m Positive integer. Size of the subsamples.
#' @param N Positive integer. Number of subsamples.
#' @param nb Positive integer. Number of bins.
#' @param g Positive real number. Pilot bandwidth.
#' @param lower Positive real number. Range over which to minimize, left bound.
#' @param upper Positive real number. Range over which to minimize, right bound.
#' @param ncores Positive integer. Number of cores with which to parallelize the computations.
#' 
#' @details
#' Bagging bootstrap bandwidth selector for the Parzen-Rosenblatt estimator.
#' 
#' @return Bagged CV bandwidth.
#' 
#' @examples
#' set.seed(1)
#' x <- rnorm(10^5)
#' hboot_bag(x, 5000, 10, 1000, lower=0.001, upper=1, ncores=2)
#'
#' @useDynLib baggingbwsel
#' @importFrom Rcpp sourceCpp
#' 
#' @export
hboot_bag <- function(x, m=n, N=1, nb=1000L, g, lower, upper, ncores=parallel::detectCores(logical=FALSE))
{
  n = length(x)
  storage.mode(x) = "double"
  sx.lst = list()
  for(i in 1:N)
  {
    sx.lst[[i]] = sample(x,m)
  }
  
  if(missing(lower) | missing(upper))
  {
    upper = 1.144*stats::sd(sx.lst[[1]])*m^(-0.2)
    lower = 0.1*upper
  }
  
  tol = 0.1*lower
  
  if(missing(g))
  {
    Rf3 = Rf3_est(sx.lst[[1]])
    g = (0.2115711/(m*Rf3))^(1/7)
  }
  cl = parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  paroutput = foreach::foreach(i=1:N,.combine=c, .noexport="x") %dopar%{
                sx = sx.lst[[i]]
                Z = bw_pair_cnts(sx, nb, m > nb/2)
                d = Z[[1L]]
                cnt = Z[[2L]]
                fboot = function(h) Cbw_boot(m, d, cnt, h, g)
                h = stats::optimize(fboot, c(lower, upper), tol=tol)$minimum
                if(h < lower + tol | h > upper - tol)
                  warning("minimum occurred at one end of the range")
                h
  }
  parallel::stopCluster(cl)
  h0 = mean(paroutput)*(m/n)^0.2
  return(h0)
}