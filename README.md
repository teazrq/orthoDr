# orthoDr

[![CRAN Version](https://www.r-pkg.org/badges/version/orthoDr)](https://CRAN.R-project.org/package=orthoDr)

Semi-parametric dimension reduction using orthogonality constrained optimization. Provides methods for regression, survival, and personalized dose finding, as well as a general-purpose solver for problems with orthogonality constraints.

## Features

- **Regression**: SIR, SAVE, PHD, locally efficient, and efficient estimation
- **Survival**: dimension reduction for censored outcomes
- **Dose finding**: personalized dose estimation via partial SAVE and direct/pseudo-direct learning
- **General optimizer**: `ortho_optim()` for any orthogonality-constrained objective
- **OpenMP** parallel gradient approximation

## Installation

```r
# From CRAN
install.packages("orthoDr")

# Development version
remotes::install_github("teazrq/orthoDr")
```

CRAN distributes version 0.6.8. This checkout is the revised development version
2.0.0; numerical results can differ between versions.

## Getting started

- [Get Started](https://teazrq.github.io/orthoDr/articles/get-started.html): installation
  and tutorials for regression, survival, personalized dose finding, and the general optimizer.
- [Runnable guide](https://teazrq.github.io/orthoDr/examples/getting-started.R) and
  [full survival example](https://teazrq.github.io/orthoDr/examples/survival-optimizer.R), including corrected
  2D and 3D plotting calls.
- [Development function reference](https://teazrq.github.io/orthoDr/reference/index.html),
  [release notes](https://teazrq.github.io/orthoDr/news/index.html), and
  [CRAN 0.6.8 reference manual](https://cran.r-project.org/web/packages/orthoDr/refman/orthoDr.html).

For documentation matching your installed version, use `help(package = "orthoDr")`.
The downloadable guide is executed by `R CMD check` through `tests/examples.R`.

## Website maintenance

The package site uses pkgdown with the same Bootstrap 5 Cosmo theme and
navigation structure as [RLT](https://teazrq.github.io/RLT/). To build it locally,
install `pkgdown` in your usual R user library, install the current package
source, then run `Rscript tools/build-site.R` from the repository root.
GitHub Actions rebuilds `docs/` and deploys it to GitHub Pages on pushes to
`master`. Pull requests build without publishing.

## References

- Wen, Z. & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. *Mathematical Programming*, 142(1-2), 397–434.
- Ma, Y. & Zhu, L. (2012). A semiparametric approach to dimension reduction. *JASA*, 107(497), 168–179.
- Ma, Y. & Zhu, L. (2013). Efficient estimation in sufficient dimension reduction. *Annals of Statistics*, 41(1), 250–268.
- Sun, Q., Zhu, R., Wang, T. & Zeng, D. (2017). Counting process based dimension reduction for censored outcomes. arXiv:1704.05046.
- Zhou, W. & Zhu, R. (2018+). Semiparametric efficient dimension reduction. arXiv:1802.06156.
