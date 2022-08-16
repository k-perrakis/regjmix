#' Function rjm
#'
#' Function for regularized joint mixture (RJM) modeling. For details see
#' Perrakis, Lartigue, Dondenlinger and Mukherjee (2022).
#' @author Konstanstinos Perrakis \email{konstantinos.perrakis@durham.ac.uk}
#' @references Perrakis, K., Lartigue, T., Dondenlinger, F. and  Mukherjee, S. (2022) Regularized joint mixture models. \url{https://arxiv.org/pdf/1908.07869.pdf}
#' @param y response vector.
#' @param X predictor matrix (do not center or standardize).
#' @param K number of groups (minimum and default is two).
#' @param outer number of EM starts (default is 10).
#' @param inner maximum number of EM iterations (default is 20).
#' @param parallel option for running the EM's in parallel, with FALSE being the default. The number of clusters must be specified by the user before calling the \code{rjm} function.
#' @param epsilon threshold for termination of the EM, based on the relative change criterion. Default is 1e-06.
#' @param method sparse method for the regression component. Available options are 'nj' (default) for the normal-Jeffreys prior and 'lasso'.
#' @param lasso.pen option to be specified when \code{method = 'lasso'} for tuning the penalty term. Available options are 'fixed' (default) which uses a two-step cross-validation procedure and 'random' which provides the MAP penalty estimate under a Pareto prior.
#' @param constant multiplicative constant (>0) for the lasso penalty when \code{lasso.pen = 'random'}. Default option is min(2p/(3n),1), where p is the number of predictors and n is the sample size.
#' @return An object with S3 class 'regjmix' allowing for call to generic functions such as \code{\link{coef}}.
#' @return An object of class 'regjmix' is a list containing the following components:
#' \item{mu}{matrix of the group-wise means of the predictor variables.}
#' \item{Sigma}{list of the group-wise covariance matrices of the predictor variables.}
#' \item{coefficients}{matrix of the group-wise intercepts and regression coefficients.}
#' \item{lambda}{the group-wise lasso penalties (when \code{method = 'lasso'}).}
#' \item{sigma2}{the group-wise error variances.}
#' \item{tau}{the group probabilities.}
#' \item{z}{the group allocations of the observations.}
#' \item{model}{summary providing number of groups, the maximum log-likelihood value, number of non-zero parameters and values of AIC and BIC.}
#' @export

rjm <-
  function(y,
           X,
           K = 2,
           outer = 10,
           inner = 20,
           parallel = FALSE,
           epsilon = 1e-06,
           method = 'nj',
           lasso.pen = 'fixed',
           constant = min(sqrt(2*ncol(X)/(3*length(y))),1)) {
    n <- length(y)
    p <- ncol(X)
    pcaX <- prcomp(X)
    eigen.values <- pcaX$sdev ^ 2
    project.dim <-
      which.min(abs(cumsum(eigen.values) / sum(eigen.values) - 0.5))
    pcaX <- pcaX$x[, 1:project.dim]
    W <- list()
    W[[1]] <- Mclust(y, G = K, verbose = FALSE)
    W[[1]] <- Mclust(y,
                     G = K,
                     modelNames = W[[1]]$modelName,
                     verbose = FALSE)
    W[[2]] <- Mclust(cbind(y, X), G = K, verbose = FALSE)
    W[[2]] <- Mclust(
      cbind(y, X),
      G = K,
      modelNames = W[[2]]$modelName,
      verbose = FALSE
    )
    W[[3]] <- Mclust(cbind(y, pcaX), G = K, verbose = FALSE)
    W[[3]] <- Mclust(
      cbind(y, pcaX),
      G = K,
      modelNames = W[[3]]$modelName,
      verbose = FALSE
    )
    tau.cand <- matrix(nrow = K, ncol = 3)
    tau.cand[, 1] <- W[[1]]$parameters$pro
    tau.cand[, 2] <- W[[2]]$parameters$pro
    tau.cand[, 3] <- W[[3]]$parameters$pro
    select <- which.max(apply(tau.cand, 2, min))
    W <- W[[select]]
    z.cand <- W$classification
    tau.cand <- W$parameters$pro
    psi <- sqrt(2 * log(p) * n) / 2
    X.ker <- matrix(NA, n, K)
    alpha <- c()
    beta <- list()
    sigma <- c()
    rho <- c()
    chi <- c()
    phi <- list()
    M <- list()
    tau <- c()
    mu <- matrix(NA, p, K)
    Sigma <- list()
    invSigma <- list()
    Q <- rep(NA, inner)
    parameters <- list()
    y.mean <- mean(y)
    y.sd <- sd(y)
    X.mean <- apply(X, 2, mean)
    X.sd <- apply(X, 2, sd)
    y <- scale(y, scale = TRUE)
    X <- scale(X, scale = TRUE)
    tX <- t(X)
    X1 <- cbind(1, X)
    tX1 <- t(X1)

    '%do.choose%' <- ifelse(parallel, `%dopar%`, `%do%`)

    if (method == 'nj') {
      EM <-
        foreach(1:outer,
                .packages = c("glassoFast", "Matrix"),
                .combine = rbind) %do.choose% {
                  # initilization

                  eff.seed <- sample(1:10 ^ 5, 1)
                  for (k in 1:K) {
                    set.seed(eff.seed)
                    tau[k] <- tau.cand[k]
                    set.seed(eff.seed)
                    mu[, k] <-
                      apply(X[z.cand == k,], 2, mean) + rnorm(p, 0, 0.5)
                    set.seed(eff.seed)
                    Sigma[[k]] <-
                      cov(X[z.cand == k,])  + diag(runif(1, 0.001, 0.0015), p)
                    #invSigma[[k]] <- chol2inv(chol(Sigma[[k]]))
                    invSigma[[k]] <- solve(Sigma[[k]], tol = .Machine$double.eps^10)
                    M[[k]] <- diag(W$z[, k])
                    term1 <- tX1 %*% M[[k]]
                    set.seed(eff.seed)
                    all.coef <-
                      c(solve(term1 %*% X1 + diag(p + 1), tol = .Machine$double.eps^10) %*% term1 %*% y) + runif(p +
                                                                                                                   1,-0.1, 0.1)
                    alpha[k] <- all.coef[1]
                    beta[[k]] <- all.coef[-1]
                    term2 <- y - alpha[k] - X %*% beta[[k]]
                    set.seed(eff.seed)
                    sigma[k] <-
                      c((t(term2) %*% M[[k]] %*% term2) / (tau[k] * n + 2) + runif(1, 0.1, 1))
                  }

                  T <- lapply(M, diag)
                  n.k <- sapply(T, sum)
                  M <- sapply(1:K, function(z)
                    Diagonal(x = T[[z]]))
                  ratio <- lapply(beta, function(z)
                    1 / z ^ 2)
                  ind <- lapply(ratio, function(z)
                    which(z != Inf))
                  p.k <- sapply(ind, length)
                  U <-
                    lapply(1:K, function(z)
                      Diagonal(x = 1 / ratio[[z]][ind[[z]]]))
                  V <- lapply(1:K, function(z)
                    Diagonal(x = ratio[[z]][ind[[z]]]))
                  X.k <- lapply(1:K, function(z)
                    X[, ind[[z]]])
                  tX.k <- lapply(X.k, t)
                  non.zero.beta <- lapply(1:K, function(z)
                    beta[[z]][ind[[z]]])


                  for (j in 1:inner) {
                    ### M-step ###

                    tau <- sapply(T, mean)
                    mu <- sapply(1:K, function(z)
                      apply(T[[z]] * X, 2, sum) / n.k[z])
                    Mat <-
                      lapply(1:K, function(z)
                        X - matrix(mu[, z], n, p, byrow = TRUE))
                    Smat <-
                      lapply(1:K, function(z)
                        t(Mat[[z]]) %*% M[[z]] %*% Mat[[z]] / n.k[z])

                    xi <- psi / (tau * n)

                    gL <- lapply(1:K, function(z)
                      glassoFast(as.matrix(Smat[[z]]), rho = xi[z]))
                    Sigma <- lapply(1:K, function(z)
                      Matrix(gL[[z]]$w))
                    invSigma <- lapply(1:K, function(z)
                      Matrix(gL[[z]]$wi))
                    logdSigma <-
                      sapply(1:K, function(z)
                        Re(sum(log(
                          eigen(Sigma[[z]], only.values = TRUE)$values
                        ))))

                    res <-
                      lapply(1:K, function(z)
                        y - alpha[z] - as.matrix(X.k[[z]]) %*% beta[[z]][ind[[z]]])
                    sigma <-
                      sapply(1:K, function(z)
                        as.numeric(t(res[[z]]) %*% M[[z]] %*% res[[z]]) / (n.k[z] + 2))
                    alpha <-
                      sapply(1:K, function(z)
                        (t(y - as.matrix(X.k[[z]]) %*% beta[[z]][ind[[z]]]) %*% T[[z]])/ n.k[z])

                    XU.k <- lapply(1:K, function(z)
                      as.matrix(X.k[[z]] %*% U[[z]]))
                    C.k <- lapply(1:K, function(z) (sigma[z] * Diagonal(x = 1 / T[[z]]) + XU.k[[z]] %*% tX.k[[z]]))
                    check.k <- unlist(lapply(1:K, function(z) sum(C.k[[z]][C.k[[z]]==Inf])))
                    if(any(check.k == Inf)) {
                      status <- TRUE
                      break }
                    check.k <- sapply(1:K, function(z) class(try(solve(C.k[[z]], tol = .Machine$double.eps^10), silent = TRUE)) == class(C.k[[z]]))
                    if(any(check.k == FALSE)) {
                      status <- TRUE
                      break }
                    G.k <-
                      lapply(1:K, function(z)
                        c(1 / sigma[z]) * U[[z]] %*% (
                          Diagonal(p.k[z]) - tX.k[[z]] %*% solve(C.k[[z]], tol = .Machine$double.eps^10) %*% XU.k[[z]]
                        ))
                    non.zero.beta <-
                      sapply(1:K, function(z)
                        c(G.k[[z]] %*% tX.k[[z]] %*% M[[z]] %*% (y - alpha[z])))
                    for (k in 1:K) {
                      if (p.k[k] != 0) {
                        beta[[k]][ind[[k]]]  <- non.zero.beta[[k]]
                        beta[[k]][-ind[[k]]] <- 0
                      } else{
                        beta[[k]] <- rep(0, p)
                      }
                    }

                    ### E-step ###

                    for (i in 1:n) {
                      res <- sapply(1:K, function(z)
                        as.matrix(X[i,] - mu[, z]))
                      X.ker[i,] <-
                        sapply(1:K, function(z)
                          as.matrix(t(res[, z]) %*% invSigma[[z]]) %*% res[, z])
                    }

                    res <-
                      sapply(1:K, function(z)
                        y - alpha[z] - as.matrix(X.k[[z]] %*% non.zero.beta[[z]]))
                    y.ker <-
                      sapply(1:K, function(z)
                        as.matrix(t(res[, z]) %*% M[[z]]) %*% res[, z] / sigma[z] + as.matrix(t(as.matrix(
                          non.zero.beta[[z]]
                        )) %*% V[[z]] %*% as.matrix(non.zero.beta[[z]])))

                    Q.k <-
                      sapply(1:K, function(z)
                        sum(T[[z]] * (
                          log(tau[z]) - 0.5 * X.ker[, z] - 0.5 * logdSigma[z]
                        )) - psi / 2 *
                          sum(abs(invSigma[[z]])) -
                          (n.k[z] + 2) * 0.5 * log(sigma[z]) - 0.5 * y.ker[z])
                    Q.cur <- sum(Q.k)
                    Q[j] <- Q.cur
                    if (j > 1) {
                      if (abs(Q[j] / Q[j - 1] - 1) < epsilon) {
                        break
                      }
                    }

                    y.ker1 <-
                      sapply(1:K, function(z)
                        as.matrix(t(res[, z]) %*% M[[z]]) %*% res[, z] / sigma[z])

                    L.k <-
                      sapply(1:K, function(z)
                        sum(T[[z]] * (
                          log(tau[z]) - 0.5 * X.ker[, z] - 0.5 * logdSigma[z]
                        ))  -
                          n.k[z] * 0.5 * log(sigma[z]) - 0.5 * y.ker1[z])
                    logLik <- sum(L.k)

                    if(K > 1){
                      for (k in 1:K) {
                        T[[k]] <- rep(NA, n)
                        q1 <- 0.5 * (log(sigma[k]) - log(sigma[-k]))
                        q2 <- (log(tau[-k]) - log(tau[k]))
                        a <- alpha[-k]
                        b <- beta[-k]
                        s <- sigma[-k]
                        q3 <-
                          0.5 * (c(y - alpha[k] - X %*% beta[[k]]) ^ 2 / sigma[k] -  sapply(
                            1:(K - 1),
                            FUN = function(z)
                              (y - a[z] - X %*% b[[z]]) ^ 2 / s[z]
                          ))
                        iS <- invSigma[c(1:K)[-k]]
                        for (i in 1:n) {
                          cX <- X[i,] - as.matrix(mu[,-k])
                          q4 <-
                            0.5 * (rep(t(c(X[i,] - mu[, k])) %*% invSigma[[k]] %*% c(X[i,] - mu[, k]), K -
                                         1) - sapply(1:(K - 1), function(z)
                                           as.matrix(t(cX[, z]) %*% iS[[z]] %*% cX[, z])))
                          q <-
                            apply(rbind(q1, q2, q3[i,], q4), 2, sum)
                          q <- sum(exp(q))
                          T[[k]][i] <- 1 / (1 + q)
                        }
                      }
                    }
                    n.k <- sapply(T, sum)

                    status <- min(n.k) < c(n / (10 * K))
                    if (status == TRUE) {
                      break
                    }

                    M <- sapply(1:K, function(z)
                      Diagonal(x = T[[z]]))
                    ratio <- lapply(beta, function(z)
                      1 / z ^ 2)
                    ind <- lapply(ratio, function(z)
                      which(z != Inf))
                    p.k <- sapply(ind, length)
                    U <-
                      lapply(1:K, function(z)
                        Diagonal(x = 1 / ratio[[z]][ind[[z]]]))
                    V <-
                      lapply(1:K, function(z)
                        Diagonal(x = ratio[[z]][ind[[z]]]))
                    X.k <- lapply(1:K, function(z)
                      X[, ind[[z]]])
                    tX.k <- lapply(X.k, t)
                  }
                  list(Q, tau, beta, sigma, mu, Sigma, T, status, alpha, logLik)
                }
    }
    if (method == 'lasso') {
      EM <-
        foreach(
          1:outer,
          .packages = c("glassoFast", "Matrix", "glmnet"),
          .combine = rbind
        ) %do.choose% {
          # initilization

          eff.seed <- sample(1:10 ^ 5, 1)
          for (k in 1:K) {
            set.seed(eff.seed)
            tau[k] <- tau.cand[k]
            set.seed(eff.seed)
            mu[, k] <-
              apply(X[z.cand == k,], 2, mean) + rnorm(p, 0, 0.5)
            set.seed(eff.seed)
            Sigma[[k]] <-
              cov(X[z.cand == k,])  + diag(runif(1, 0.001, 0.0015), p)
            #invSigma[[k]] <- chol2inv(chol(Sigma[[k]]))
            invSigma[[k]] <- solve(Sigma[[k]], tol = .Machine$double.eps^10)
            M[[k]] <- diag(W$z[, k])
            term1 <- tX1 %*% M[[k]]
            set.seed(eff.seed)
            all.coef <-
              c(solve(term1 %*% X1 + diag(p + 1), tol = .Machine$double.eps^10) %*% term1 %*% y) + runif(p +
                                                                                                           1,-0.1, 0.1)
            alpha[k] <- all.coef[1]
            beta[[k]] <- all.coef[-1]
            term2 <- y - alpha[k] - X %*% beta[[k]]
            set.seed(eff.seed)
            sigma[k] <-
              c((t(term2) %*% M[[k]] %*% term2) / (tau[k] * n + 2) + runif(1, 0.1, 1))
            rho[k] <- 1 / sqrt(sigma[k])
            chi[k] <- alpha[k] * rho[k]
            phi[[k]] <- beta[[k]] * rho[k]
          }

          T <- lapply(M, diag)
          n.k <- sapply(T, sum)
          M <- sapply(1:K, function(z)
            Diagonal(x = T[[z]]))
          sqrt2 <- sqrt(2)
          z.status <- matrix(nrow = n, ncol = inner)
          assign.status <- FALSE
          count <- 0

          for (j in 1:inner) {
            ### M-step ###

            tau <- sapply(T, mean)
            mu <- sapply(1:K, function(z)
              apply(T[[z]] * X, 2, sum) / n.k[z])
            Mat <-
              lapply(1:K, function(z)
                X - matrix(mu[, z], n, p, byrow = TRUE))
            Smat <-
              lapply(1:K, function(z)
                t(Mat[[z]]) %*% M[[z]] %*% Mat[[z]] / n.k[z])

            xi <- psi / (tau * n)

            gL <- lapply(1:K, function(z)
              glassoFast(as.matrix(Smat[[z]]), rho = xi[z]))
            Sigma <- lapply(1:K, function(z)
              Matrix(gL[[z]]$w))
            invSigma <- lapply(1:K, function(z)
              Matrix(gL[[z]]$wi))
            logdSigma <-
              sapply(1:K, function(z)
                Re(sum(log(
                  eigen(Sigma[[z]], only.values = TRUE)$values
                ))))

            eff.sample <-
              lapply(1:K, function(z)
                which(T[[z]] != 0))
            eff.n <- sapply(eff.sample, length)

            scale.M <-
              lapply(1:K, function(z)
                sqrt(M[[z]][eff.sample[[z]], eff.sample[[z]]])/2)
            yy <-
              sapply(1:K, function(z)
                scale.M[[z]] %*% y[eff.sample[[z]]] * rho[z] - diag(scale.M[[z]]) * chi[z])
            XX <-
              sapply(1:K, function(z)
                scale.M[[z]] %*% X[eff.sample[[z]],])

            if (lasso.pen == 'fixed') {
              if (j == 1 | (assign.status == TRUE & count == 1)) {
                lambda <- sapply(1:K, function(z)
                  cv.glmnet(
                    x = as.matrix(XX[[z]]) / eff.n[z],
                    y = yy[[z]],
                    family = 'gaussian',
                    standardize = TRUE,
                    intercept = FALSE
                  )$lambda.min)
              }
            }
            if (lasso.pen == 'random') {
              lambda <-
                sapply(1:K, function(z)
                  constant / sum(abs(phi[[z]])) * sqrt2 * sqrt(K*log(p) / n))
              lambda[lambda==Inf] <- n
            }

            y.eff <-
              lapply(1:K, function(z)
                y[eff.sample[[z]]])
            X.eff <-
              lapply(1:K, function(z)
                X[eff.sample[[z]],])

            inner.product.y <- sapply(1:K, function(z)
              t(y.eff[[z]]) %*% M[[z]][eff.sample[[z]],eff.sample[[z]]] %*% y.eff[[z]])
            inner.product.yX <-
              sapply(1:K, function(z)
                t(y.eff[[z]]) %*% M[[z]][eff.sample[[z]],eff.sample[[z]]] %*% (rep(chi[z], eff.n[z]) + X.eff[[z]] %*% phi[[z]]))
            inner.product.y <-
              unlist(lapply(inner.product.y, as.vector))
            inner.product.yX <-
              unlist(lapply(inner.product.yX, as.matrix))
            rho <-
              sapply(1:K, function(z)
                (
                  sqrt(inner.product.yX[z] ^ 2 + 4 * (n.k[z] + p + 2) * inner.product.y[z]) +
                    inner.product.yX[z]
                ) / (2 * inner.product.y[z]))
            sigma <- 1 / rho ^ 2

            chi <-
              sapply(1:K, function(z)
                (t(rho[z] * y - X %*% phi[[z]]) %*% T[[z]]) / n.k[z])
            alpha <- chi / rho

            yy <-
              sapply(1:K, function(z)
                scale.M[[z]] %*% y[eff.sample[[z]]] * rho[z] - diag(scale.M[[z]]) * chi[z])

            phi <-
              lapply(1:K, function(z)
                glmnet(
                  x = as.matrix(XX[[z]]) / eff.n[z],
                  y = yy[[z]],
                  family = 'gaussian',
                  lambda = lambda[z],
                  standardize = TRUE,
                  intercept = FALSE
                )$beta)

            phi <- lapply(phi, as.numeric)
            phi <- lapply(1:K, function(z)
              phi[[z]] / eff.n[z])

            beta <- lapply(1:K, function(z)
              phi[[z]] / rho[z])
            beta <- lapply(beta, as.numeric)

            ind <- lapply(1:K, function(z)
              which(phi[[z]] != 0))
            non.zero.phi <- lapply(1:K, function(z)
              phi[[z]][ind[[z]]])
            X.k <- lapply(1:K, function(z)
              as.matrix(XX[[z]][, ind[[z]]]))

            ### E-step ###

            for (i in 1:n) {
              res <- sapply(1:K, function(z)
                as.matrix(X[i,] - mu[, z]))
              X.ker[i,] <-
                sapply(1:K, function(z)
                  as.matrix(t(res[, z]) %*% invSigma[[z]]) %*% res[, z])
            }

            res <-
              sapply(1:K, function(z)
                yy[[z]] - as.matrix(X.k[[z]] %*% non.zero.phi[[z]]))
            y.ker <-
              sapply(1:K, function(z)
                as.matrix(t(res[[z]]) %*% res[[z]]) + lambda[z] * sum(abs(non.zero.phi[[z]])))

            if (lasso.pen == 'fixed') {
              Q.k <-
                sapply(1:K, function(z)
                  sum(T[[z]] * (
                    log(tau[z]) - 0.5 * X.ker[, z] - 0.5 * logdSigma[z]
                  )) - psi / 2 *
                    sum(abs(invSigma[[z]])) +
                    (n.k[z] + p + 2) * log(rho[z]) -
                    y.ker[z] / 2)
            }

            if (lasso.pen == 'random') {
              Q.k <-
                sapply(1:K, function(z)
                  sum(T[[z]] * (
                    log(tau[z]) - 0.5 * X.ker[, z] - 0.5 * logdSigma[z]
                  )) - psi / 2 *
                    sum(abs(invSigma[[z]])) +
                    (n.k[z] + p + 2) * log(rho[z]) -
                    y.ker[z] /2 + constant * sqrt2 * sqrt(K*log(p) / n) * log(lambda[z]))
            }

            Q.cur <- sum(Q.k)
            Q[j] <- Q.cur
            if (j > 1) {
              if (abs(Q[j] / Q[j - 1] - 1) < epsilon) {
                break
              }
            }

            y.ker1 <-
              sapply(1:K, function(z)
                as.matrix(t(res[[z]]) %*% res[[z]]))

            L.k <-
              sapply(1:K, function(z)
                sum(T[[z]] * (
                  log(tau[z]) - 0.5 * X.ker[, z] - 0.5 * logdSigma[z]
                ))  -
                  n.k[z] * 0.5 * log(rho[z]) - 0.5 * y.ker1[z])
            logLik <- sum(L.k)

            if(K > 1){
              for (k in 1:K) {
                T[[k]] <- rep(NA, n)
                q1 <- 0.5 * (log(sigma[k]) - log(sigma[-k]))
                q2 <- (log(tau[-k]) - log(tau[k]))
                a <- alpha[-k]
                b <- beta[-k]
                s <- sigma[-k]
                q3 <-
                  0.5 * (c(y  - alpha[k] - X %*% beta[[k]]) ^ 2 / sigma[k] -  sapply(
                    1:(K - 1),
                    FUN = function(z)
                      (y - a[z] - X %*% b[[z]]) ^ 2 / s[z]
                  ))
                iS <- invSigma[c(1:K)[-k]]
                for (i in 1:n) {
                  cX <- X[i,] - as.matrix(mu[,-k])
                  q4 <-
                    0.5 * (rep(t(c(X[i,] - mu[, k])) %*% invSigma[[k]] %*% c(X[i,] - mu[, k]), K -
                                 1) - sapply(1:(K - 1), function(z)
                                   as.matrix(t(cX[, z]) %*% iS[[z]] %*% cX[, z])))
                  q <- apply(rbind(q1, q2, q3[i,], q4), 2, sum)
                  q <- sum(exp(q))
                  T[[k]][i] <- 1 / (1 + q)
                }
              }
            }
            z.status[, j] <-
              apply(matrix(unlist(T), n, K), 1, which.max)
            if (j > 3) {
              assign.status <- all.equal(z.status[, j], z.status[, j - 1])
              if (assign.status == TRUE)
                count <- count + 1
            }
            n.k <- sapply(T, sum)

            status <- min(n.k) < c(n / (10 * K))
            if (status == TRUE) {
              break
            }

            M <- sapply(1:K, function(z)
              Diagonal(x = T[[z]]))
          }
          list(Q, tau, beta, sigma, mu, Sigma, T, status, lambda, alpha, logLik)
        }
    }

    Q <- matrix(unlist(EM[, 1]), outer, inner, byrow = TRUE)
    status <- unlist(EM[, 8], use.names = FALSE)
    if (sum(status) != outer) {
      status <- which(status == FALSE)
      Q <- Q[status, ]
      outer <- length(status)
      if (is.matrix(Q) == TRUE) {
        last.iter <- apply(Q, 1, function(z)
          sum(is.na(z) == FALSE))
        Q.last.iter <- c()
        for (i in 1:outer) {
          Q.last.iter[i] <- Q[i, last.iter[i]]
        }
        max.Q <- which.max(Q.last.iter)
      } else {
        max.Q <- status
      }
      if (method == 'nj') {
        alphas <- unlist(EM[max.Q, 9], use.names = FALSE)
        logLik <- unlist(EM[max.Q, 10], use.names = FALSE)
      }
      if (method == 'lasso') {
        alphas <- unlist(EM[max.Q, 10], use.names = FALSE)
        logLik <- unlist(EM[max.Q, 11], use.names = FALSE)
      }
      betas <-
        Matrix(unlist(lapply(EM[max.Q, 3], c)), p, K, sparse = TRUE)
      if (method == 'nj') {
        betas[abs(betas) <= 1e-20] <- 0
        lambdas <- 'No global penalty with normal-Jeffreys'
      } else {
        lambdas <- unlist(EM[max.Q, 9], use.names = FALSE)
      }
      p.beta <- sum(betas!=0)
      b0 <- c()
      for (k in 1:K) {
        b0[k] <-
          y.mean + y.sd * (alphas[k] - sum(betas[, k] * X.mean / X.sd))
      }
      betas <- betas * y.sd / X.sd
      betas <- rbind(b0, betas)
      rownames(betas) <- c('Intercept', paste0("b", c(1:p)))
      Sigmas <- unlist(EM[max.Q, 6], use.names = FALSE)
      p.upperSigma <- lapply(1:k, function(z)
        sum(Sigmas[[z]][upper.tri(Sigmas[[z]])]!=0))
      p.upperSigma <- sum(unlist(p.upperSigma))
      if(method == 'lasso' & lasso.pen == 'random'){
        n.par <- ((4 + 2*p)*K + p.beta + p.upperSigma)
        bic <- 2*logLik - log(n) * n.par
        aic <- 2*logLik - 2 * n.par
      } else {
        n.par <- ((3 + 2*p)*K + p.beta + p.upperSigma)
        bic <- 2*logLik - log(n) * n.par
        aic <- 2*logLik - 2 * n.par
      }
      z <- matrix(unlist(EM[max.Q, 7]), n, K)
      z <- apply(z, 1, which.max)
      model.details = data.frame(round(K), logLik, n.par, aic, bic)
      colnames(model.details) = c('Clusters', 'logLikelihood', 'No.Parameters', 'AIC', 'BIC')
      rownames(model.details) = c('')
      results <-
        list(
          mu = matrix(unlist(lapply(EM[max.Q, 5], c)), p, K),
          Sigma = Sigmas,
          coefficients = betas,
          lambda = lambdas,
          sigma2 = unlist(EM[max.Q, 4], use.names = FALSE),
          tau = unlist(EM[max.Q, 2], use.names = FALSE),
          z = z,
          model = model.details
        )
      class(results) <- "regjmix"
    } else {
      message('Produced NULL output: all EMs were terminated because of small group sample size')
      results <- NULL
    }
    results
  }
