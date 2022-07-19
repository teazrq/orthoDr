
# orthoDr

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/orthoDr)](https://CRAN.R-project.org/package=orthoDr)
[![](https://cranlogs.r-pkg.org/badges/orthoDr)](https://cran.r-project.org/package=orthoDr)
<!-- badges: end -->
  
The goal of `orthoDr` is to use an orthogonality constrained optimization
algorithm to solve a variety of dimension reduction problems in the
semiparametric framework.

## Installation

You can install the released version of `orthoDr` from [CRAN](https://CRAN.R-project.org/package=orthoDr) with:

``` r
  install.packages("orthoDr")
```

## Implemented Methods

This package implements the orthogonality constrained (Stiefel manifold) optimization approach proposed by [Wen & Yin (2013)](https://link.springer.com/article/10.1007/s10107-012-0584-1). A drop-in solver `ortho_optim()` works just the same as the `optim()` function. Relying on this optimization approach, we also implemented a collection of dimension reduction models for survival analysis, regression, and personalized medicine. 

  * CP-SIR, Forward, IR-CP and IR-Semi methods in [Sun, Zhu, Wang & Zeng (2019)](https://arxiv.org/abs/1704.05046)
  * semi-SIR and semi-PHD in [Ma & Zhu (2012)](https://www.tandfonline.com/doi/full/10.1080/01621459.2011.646925)
  * Direct and pseudo-Direct methods in [Zhou, Zhu & Zeng (2021)](https://arxiv.org/abs/1802.06156)

We also implemented several methods and functions for comparison, testing and utilization purposes  
    
  * `hMave`: This is a direct `R` translation of the hMave `MATLAB` code by [Xia, Zhang & Xu (2010)](https://www.tandfonline.com/doi/abs/10.1198/jasa.2009.tm09372)
  * `pSAVE`: partial-SAVE in [Feng, Wen, Yu & Zhu (2013)](https://www.tandfonline.com/doi/full/10.1080/01621459.2012.746065)
  * `dist_cross()`: kernel distances matrix between two sets of data, as an extension of `dist()`
  * `distance()`: distance correlation between two linear spaces
  * `silverman()`: Silverman's rule of thumb bandwidth
  
