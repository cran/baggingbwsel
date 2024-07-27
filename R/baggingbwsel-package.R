#' baggingbwsel: Bagging bandwidth selection in kernel density and 
#' regression estimation
#'
#' This package implements bagging bandwidth selection methods for the 
#' Parzen-Rosenblatt kernel density estimator, and for the Nadaraya-Watson and 
#' local polynomial kernel regression estimators. 
#' These bandwidth selectors can achieve greater statistical precision than their 
#' non-bagged counterparts while being computationally fast. 
#' See Barreiro-Ures et al. (2021a) and Barreiro-Ures et al. (2021b).
#' @name baggingbwsel-package
# @aliases baggingbwsel
# @docType package
#' @useDynLib baggingbwsel
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach %dopar%
#' @importFrom grDevices topo.colors
#' @importFrom graphics box contour image lines par persp points polygon rug
#' @importFrom stats dnorm quantile var
#' @references
#' Barreiro-Ures, D., Cao, R., Francisco-Fernández, M., & Hart, J. D. (2021a). 
#' Bagging cross-validated bandwidths with application to big data. 
#' \emph{Biometrika}, \bold{108}(4), 981-988, \doi{10.1093/biomet/asaa092}.
#' 
#' Barreiro-Ures, D., Cao, R., & Francisco-Fernández, M. (2021b).  
#' Bagging cross-validated bandwidth selection in nonparametric regression 
#' estimation with applications to large-sized samples.
#' \emph{arXiv preprint}, \doi{10.48550/arXiv.2105.04134}.
"_PACKAGE"


