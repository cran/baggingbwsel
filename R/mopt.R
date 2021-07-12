#' Estimation of the optimal subsample size for bagged CV bandwidth for Parzen-Rosenblatt estimator
#'
#' @param x Vector. Sample.
#' @param N Positive integer. Number of subsamples for the bagged bandwidth.
#' @param r Positive integer. Size of the subsamples.
#' @param s Positive integer. Number of subsamples.
#' @param ncores Positive integer. Number of cores with which to parallelize the computations.
#' 
#' @details
#' Estimates the optimal size of the subsamples for the bagged CV bandwidth selector for the Parzen-Rosenblatt estimator.
#' 
#' @return Estimate of the optimal subsample size.
#' 
#' @examples
#' set.seed(1)
#' x <- rt(10^5, 5)
#' mopt(x, 500, 500, 10, 2)
#' 
#' @export
mopt <-
function(x,N,r=1000,s=100,ncores=parallel::detectCores())
{
  
  n <- length(x)
  
  sx.lst = list()
  for(i in 1:s)
  {
    sx.lst[[i]] = sample(x,r)
  }

  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  paroutput <- foreach::foreach(i=1:s,.combine=cbind,.noexport="x") %dopar%{
    subx = sx.lst[[i]]
    fit = mclust::Mclust(subx)
    mu <- fit$parameters$mean
    sigma <- sqrt(fit$parameters$variance$sigmasq)
    w <- fit$parameters$pro
    fhat <- Vectorize(function(x)sum(w/sigma*stats::dnorm((x-mu)/sigma)))
    fhat2 <- Vectorize(function(x)sum(w/sigma^3*kedd::kernel.fun((x-mu)/sigma,2)$kx))
    fhat3 <- Vectorize(function(x)sum(w/sigma^4*kedd::kernel.fun((x-mu)/sigma,3)$kx))
    Af <- stats::integrate(function(x)fhat(x)^2,-Inf,Inf)$value
    Af2 <- stats::integrate(function(x)fhat2(x)^2,-Inf,Inf)$value
    Af3 <- stats::integrate(function(x)fhat3(x)^2,-Inf,Inf)$value
    c1 <- 0.07019858*Af3/(Af2^(8/5))
    Astar <- 0.179569*Af/(Af2^(3/5))
    muCV <- -0.3469317*Af/Af2^(2/5)
    return(list(Astar, muCV, c1))
  }
  parallel::stopCluster(cl)
  
  
  Astar = mean(unlist(paroutput[1,]))
  muCV = mean(unlist(paroutput[2,]))
  c1 = mean(unlist(paroutput[3,]))
  
  mse <- function(m)Astar*m^(-0.2)*n^(-0.4)*(1/N+(N-1)/N*(m/n)^2)+m^(-0.4)*n^(-0.4)*(muCV+c1*m^(-0.2))^2
  m0 <- floor(stats::optimize(mse,c(1,n))$minimum)
  return(m0)

}