
<!-- 
README.md is generated from README.Rmd. 
Please edit that file 
-->

# baggingbwsel: Bagging bandwidth selection in kernel density and regression estimation

### Version 1.1

This package implements bagging bandwidth selection methods for the
Parzen-Rosenblatt kernel density estimator, and for the Nadaraya-Watson
and local polynomial kernel regression estimators. These bandwidth
selectors can achieve greater statistical precision than their
non-bagged counterparts while being computationally fast. See
Barreiro-Ures et al. (2021a) and Barreiro-Ures et al. (2021b).

## Installation

`baggingbwsel` is not yet available from CRAN, but you can install the
development version from github with:

``` r
# install.packages("remotes")
remotes::install_github("rubenfcasal/baggingbwsel")
```

Note also that, as this package requires compilation, Windows users need
to have previously installed the appropriate version of
[Rtools](https://cran.r-project.org/bin/windows/Rtools/), and OS X users
need to have installed
[Xcode](https://apps.apple.com/us/app/xcode/id497799835).

## Authors

- Daniel Barreiro-Ures (<daniel.barreiro.ures@udc.es>)

- Ruben Fernandez-Casal (<rubenfcasal@gmail.com>)

- Jeffrey Hart

- Ricardo Cao

- Mario Francisco-Fernandez

**Maintainer**: [Ruben Fernandez-Casal](https://rubenfcasal.github.io)
(Dep. Mathematics, University of A Coruña, Spain). Please send comments,
error reports or suggestions to <rubenfcasal@gmail.com>.

## References

- Barreiro-Ures, D., Cao, R., Francisco-Fernández, M., & Hart, J. D.
  (2021a). [Bagging cross-validated bandwidths with application to big
  data](https://doi.org/10.1093/biomet/asaa092). *Biometrika*,
  **108**(4), 981-988, .

- Barreiro-Ures, D., Cao, R., & Francisco-Fernández, M. (2021b).
  [Bagging cross-validated bandwidth selection in nonparametric
  regression estimation with applications to large-sized
  samples](https://arxiv.org/abs/2105.04134). *arXiv preprint*.
