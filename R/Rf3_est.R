#' @import mclust
Rf3_est = function(x)
{
  n = length(x)
  fit = mclust::Mclust(x,verbose=FALSE)$parameters
  mu = fit$mean
  sig = sqrt(fit$variance$sigmasq)
  w = fit$pro
  mix = nor1mix::norMix(mu=mu,sigma=sig,w=w)
  q1 = nor1mix::qnorMix(0.001,mix)
  q2 = nor1mix::qnorMix(0.999,mix)
  f3 = Vectorize(function(x)sum(w/sig^4*kedd::kernel.fun((x-mu)/sig,3)$kx))
  Rf3 = stats::integrate(function(x)f3(x)^2,q1,q2)$value
  return(Rf3)
}