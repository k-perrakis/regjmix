# regjmix
R package for fitting Regularized Joint Mixture (RJM) models for continuous response and predictor variables via an EM algorithm. RJMs provide simultaneous Gaussian modeling of the predictors conditional on the latent groups (X|z) and of the response variable conditional on the predictors and the latent groups (y|X,z). Regularization and sparsity is achieved in both components of the model via graphical lasso for the covariances of the predictors and either (i) lasso-type shrinkage or (ii) normal-Jeffreys shrinkage for the coefficients of the regression component.

To install and load regjmix:
```{r}
library(devtools)
install_github('k-perrakis/regjmix')
library(regjmix)
```
