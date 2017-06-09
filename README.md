# Shape-constrained unfolding

Matlab code for the paper M. Kuusela and P. B. Stark (2017), Shape-constrained uncertainty quantification in unfolding steeply falling elementary particle spectra, [arXiv:1512.00905](http://arxiv.org/abs/1512.00905) [stat.AP], to appear in the Annals of Applied Statistics.

Tested on Matlab R2014a. Requires the Optimization Toolbox and the Curve Fitting Toolbox.

Note: The interior-point algorithm in `linprog` and the SQP algorithm in `fmincon` were changed in Matlab R2015b and R2016b, respectively. The present code may not be fully compatible with the new solvers, but should work with the legacy implementations of these algorithms.

Uses the following external functions: [parfor_progress](http://ch.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor--progress-bar--that-works-with-parfor), [subplot_tight](http://ch.mathworks.com/matlabcentral/fileexchange/30884-controllable-tight-subplot/) and [dashline](http://ch.mathworks.com/matlabcentral/fileexchange/1892-dashline)

## Contents
### Main files

| File | Description |
| --- | --- |
| coverageTable.m | Produces Table 1 in the paper |
| generateDataFlat.m | Generates the data for the constant spectrum |
| generateDataIncJets.m | Generates the data for the inclusive jet case study |
| generateDataLinear.m | Generates the data for the linearly decreasing spectrum |
| plotFlat.m | Plots strict bounds intervals for the constant spectrum |
| plotLinear.m | Plots strict bounds intervals for the linearly decreasing spectrum |
| plotPaper.m | Produces the figures for the paper |
| splineDemo.m | Shape-constrained spline fit, *used only for the demo in Figure 1* |
| unfoldDAgostini.m | D'Agostini iteration |
| unfoldDAgostiniCVVar.m | D'Agostini iteration with weighted cross-validation |
| unfoldFlatStrictBounds.m | Shape-constrained strict bounds for the constant spectrum |
| unfoldFlatStrictBoundsCoverage.m | Coverage study of shape-constrained strict bounds for the constant spectrum |
| unfoldIncJetsDAgostiniCoverage.m | Coverage study for D'Agostini with fixed stopping point |
| unfoldIncJetsDAgostiniCoverageCVVar.m | Coverage study for D'Agostini when the stopping point is chosen using weighted cross-validation |
| unfoldIncJetsSVDCoverage.m | Coverage study for SVD unfolding with fixed regularization strength |
| unfoldIncJetsSVDCoverageCVVar.m | Coverage study for SVD unfolding when the regularization strength is chosen using weighted cross-validation |
| unfoldIncJetsStrictBounds.m | Shape-constrained strict bounds for the inclusive jet pT spectrum |
| unfoldIncJetsStrictBoundsCoverage.m | Coverage study of shape-constrained strict bounds for the inclusive jet pT spectrum |
| unfoldIncJetsStrictBoundsInit.m | Precomputation for the strict bounds, *needs to be run before the other strict bounds scripts* |
| unfoldLinearStrictBounds.m | Shape-constrained strict bounds for the linearly decreasing spectrum |
| unfoldLinearStrictBoundsCoverage.m | Coverage study of shape-constrained strict bounds for the linearly decreasing spectrum |
| unfoldSVD.m | SVD unfolding |
| unfoldStrictBounds.m | Unfolded confidence intervals using shape-constrained strict bounds |
| unfoldStrictBoundsNoConOc.m | As above, but omits computing the conservatively discretized convex intervals |

### Helper files

| File | Description |
| --- | --- |
| CVVarObjectiveSVD.m | Weighted cross-validation objective function for SVD |
| computeBetaInit.m | Computes starting point for the spline fit |
| constraintConLbTaylor2.m | Nonlinear constraints for the convex lower bound |
| constraintConUbTaylor2.m | Nonlinear constraints for the convex upper bound |
| funCon.m | Objective function for fmincon |
| incJetsTheory.m | The true inclusive jet pT spectrum |
| logLikelihood.m | Log-likelihood for the spline fit |
| makeFeasible.m | Adjusts the solution iteratively to make it feasible |
| unfoldDAgostiniOneIter.m | One iteration of D'Agostini |

### Directories

| Directory | Description |
| --- | --- |
| data | Simulated datasets, generated using generateData*.m  |
| figures | Figures for the paper, produced using plotPaper.m |
| results | Saved results from running the scripts |

### Copyright (c) 2016-2017 Mikael Kuusela and Philip B. Stark
