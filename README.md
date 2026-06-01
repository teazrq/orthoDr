# orthoDr

[![CRAN Version](https://www.r-pkg.org/badges/version/orthoDr)](https://CRAN.R-project.org/package=orthoDr)
[![R-CMD-check](https://github.com/rqzhu-aide/orthoDr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rqzhu-aide/orthoDr/actions/workflows/R-CMD-check.yaml)

Semi-parametric dimension reduction using orthogonality constrained optimization. Provides methods for regression, survival, and personalized dose finding, as well as a general-purpose solver for problems with orthogonality constraints.

## Features

- **Regression**: SIR, SAVE, pHd, MAVE, and semiparametric efficient estimation (SEff)
- **Survival**: dimension reduction for censored outcomes
- **Dose finding**: personalized dose estimation via partial SAVE and direct/pseudo-direct learning
- **General optimizer**: `ortho_optim()` for any orthogonality-constrained objective
- **OpenMP** parallel gradient approximation

## Installation

```r
# From CRAN
install.packages("orthoDr")

# Development version
remotes::install_github("rqzhu-aide/orthoDr")
```

## References

- Wen, Z. & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. *Mathematical Programming*, 142(1-2), 397–434.
- Ma, Y. & Zhu, L. (2012). A semiparametric approach to dimension reduction. *JASA*, 107(497), 168–179.
- Ma, Y. & Zhu, L. (2013). Efficient estimation in sufficient dimension reduction. *Annals of Statistics*, 41(1), 250–268.
- Sun, Q., Zhu, R., Wang, T. & Zeng, D. (2017). Counting process based dimension reduction for censored outcomes. arXiv:1704.05046.
- Zhou, W. & Zhu, R. (2018+). Semiparametric efficient dimension reduction. arXiv:1802.06156.
