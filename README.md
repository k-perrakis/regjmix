# regjmix

R package for fitting Regularized Joint Mixture (RJM) models for continuous response and predictor variables via an EM algorithm. RJMs provide simultaneous Gaussian modeling of the predictors conditional on the latent groups (X\|z) and of the response variable conditional on the predictors and the latent groups (y\|X,z). Regularization and sparsity is achieved in both components of the model via graphical lasso for the covariances of the predictors and either (i) lasso-type shrinkage or (ii) normal-Jeffreys shrinkage for the coefficients of the regression component.

## Installation

To install and load *regjmix*:

``` r
library(devtools)
install_github('k-perrakis/regjmix')
library(regjmix)
```

## Synthetic example

A simple simulation with two latent groups, total sample size equal to 100 and 10 potential predictors where only the first predictor has an effect on the response, with different regression coefficients in the two groups. 

``` r
n = 100
p = 10
z.latent = c(rep(1, n/2), rep(2, n/2))
set.seed(1)
X = matrix(rnorm(n*p), n, p)
y = c()
y[z.latent==1] =  X[z.latent==1, 1] + rnorm(n/2, 0, 0.1)
y[z.latent==2] = -X[z.latent==2, 1] + rnorm(n/2, 0, 0.1)
```

The function for fitting RJM models is *rjm*. Below we fit a model using lasso regression, using five EM starts, and extract the regression coefficients and the cluster labels.

``` r
rjm.lasso = rjm(y, X, method = 'lasso', outer = 5)
rjm.lasso$model
##  Clusters logLikelihood No.Parameters       AIC       BIC
##         2     -641.4234            85 -1452.847 -1674.286
coef(rjm.lasso)  # alternatively use rjm.lasso$coefficients
## 11 x 2 sparse Matrix of class "dgCMatrix"
##                                     
## Intercept -0.0108038722 0.0004473664
## b1        -1.0005646347 0.9696923453
## b2         .            .           
## b3         .            .           
## b4         .            .           
## b5         .            .           
## b6         .            .           
## b7         .            .           
## b8         .            .           
## b9        -0.0008425582 .           
## b10        .            .
rjm.lasso$z
##   [1] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 2 1 2 2 2 1 2 2 2 2 2
##  [44] 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1
##  [87] 1 1 1 1 1 1 1 1 1 1 1 1 1 1
adjustedRandIndex(z.latent, rjm.lasso$z)
## [1] 0.808062
```

The EM can also be run in parallel as follows.

``` r
cl = detectCores() - 1
cl = makeCluster(cl)
registerDoParallel(cl)
rjm.lasso = rjm(y, X, method = 'lasso', outer = 5, parallel = TRUE)
stopCluster(cl)
```

For additional details:

``` r
?rjm 
```

## Reference
Perrakis, K., Lartigue, T., Dondelinger, F. and Mukherjee, S. (2022). Regularized joint mixture models, [https://arxiv.org/pdf/1908.07869.pdf].
