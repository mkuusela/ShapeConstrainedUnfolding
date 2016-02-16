# Shape Constrained Unfolding

Matlab code for the paper M. Kuusela and P. B. Stark. Shape-constrained uncertainty quantification in unfolding steeply falling elementary particle spectra. [arXiv:1512.00905](http://arxiv.org/abs/1512.00905)

Uses the following external functions:

-[parfor_progress](http://ch.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor)
-[subplot_tight](http://ch.mathworks.com/matlabcentral/fileexchange/30884-controllable-tight-subplot/)
-[dashline](http://ch.mathworks.com/matlabcentral/fileexchange/1892-dashline)

## Contents
### Main files

| File | Description |
| --- | --- |
| generateDataIncJets.m | Generates the data for the inclusive jet case study |
| plotPaper.m | Generates the figures for the paper |
| splineDemo.m | Spline fit for Figure 1 |
| unfoldDAgostini.m | D'Agostini iteration |
| unfoldDAgostiniCVVar.m | D'Agostini iteration with weighted crossvalidation |
| unfoldIncJetsDAgostiniCoverage.m | Coverage study for D'Agostini with fixed stopping point |
| unfoldIncJetsDAgostiniCoverageCVVar.m | Coverage study for D'Agostini when the stopping point is chosen using weighted crossvalidation |
| unfoldIncJetsSVDCoverage.m | Coverage study for SVD unfolding with fixed regularization strength |
| unfoldIncJetsSVDCoverageCVVar.m | Coverage study for SVD unfolding when regularization strength is chosen using weighted crossvalidation |
| unfoldIncJetsStrictBounds.m | Shape-constrained strict bounds for the inclusive jet spectrum |
| unfoldIncJetsStrictBoundsCoverage.m | Coverage study for shape-constrained strict bounds |
| unfoldIncJetsStrictBoundsInit.m | Precomputes values for the strict bounds, needs to be run before the main script |
| unfoldSVD.m | SVD unfolding |
| unfoldStrictBounds.m | Unfolded confidence intervals using shape-constrained strict bounds |
| unfoldStrictBoundsNoConOc.m | As above, but omits computing the conservatively discretized convex intervals |

### Helper files

| File | Description |
| --- | --- |
| CVVarObjectiveSVD.m | Weighted crossvalidation objective function for SVD |
| computeBetaInit.m | Compute initial point for the spline fit |
| constraintConLbTaylor2.m | Nonlinear constraint for the convex lower bound |
| constraintConUbTaylor2.m | Nonlinear constraint for the convex upper bound |
| funCon.m | Objective function for fmincon |
| incJetsTheory.m | The true jet spectrum |
| logLikelihood.m | Log-likelihood for the spline fit |
| makeFeasible.m | Adjusts the solution iteratively to make it feasible |
| unfoldDAgostiniOneIter.m | One iteration of D'Agostini |

Copyright (c) 2016 Mikael Kuusela and Philip B. Stark