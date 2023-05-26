labpmm2 <- function(y, t, r, numUnit, obsPerUnit, 
                    gammaA = 0.01, gammaB = 0.01, 
                    fixEffVar = 10^4){
  library(Matrix)
  # y: the response variable (count data)
  # t: the predictive variables
  # numUnit: number of groups/units, aka number of random intercepts 
  # in case of one random effect intercept for each group/unit
  # obsPerUnit: number of observations per group/unit
  # numFixEff: number of fixed effect variables
  # gammaA, gammaB: parameters of the prior precision eta ~ Gamma(a,b)
  # fixEffVar: variance of the prior of fixed effect B ~ N(0, fixEffVar*I)
  
  #require r and t being matrices
  n = length(y)
  r = matrix(r, nrow = n)
  t = matrix(t,nrow = n)
  numFixEff = ncol(t)+1
  numRanEff = ncol(r)+1
  
  # Design matrix for the random intercept
  # z1: column vector of ones
  z1 = matrix(rep(1,obsPerUnit), ncol=1)
  # z: design matrix of random effect intercept b0
  z = kronecker(diag(numUnit), z1)
  # u: design matrix of random slope b1
  u = matrix(0, nrow(r), numUnit)
  j1 = 1
  for (i in 1:numUnit){
    j2 = j1 + obsPerUnit - 1
    u[j1:j2, i] = r[j1:j2]
    j1 = j2 + 1
  }
  ones = rep(1, length(y))
  t = cbind(ones, t,z,u)
  t = as.matrix(t)
  
  ## Calculate likelihood of data given eta
  
  get_minus_log_llh <- function(phi, inv.cov.mat){
    # return minus log likelihood of model given eta and data
    # phi: vector of fixed effects and random effects
    # inv.cov.mat: inverse of covariance matrix of phi which is known given eta
    - sum(t %*% phi * y - exp(t%*%phi)) + 1/2 * phi %*% inv.cov.mat %*% phi
  }
  
  ## Gradient and hessian of f
  
  get_gradient <- function(phi, inv.cov.mat){
    # return the gradient vector of log likelihood of model given eta and data
    # phi: vector of fixed effects and random effects
    # inv.cov.mat: inverse of covariance matrix of phi which is known given eta
    - colSums(t*y - t*c(exp(t%*%phi))) + inv.cov.mat%*%phi
  }
  
  get_hessian <- function(phi, inv.cov.mat){
    # return the hessian matrix of log likelihood of model given eta and data
    # phi: vector of fixed effects and random effects
    # inv.cov.mat: inverse of covariance matrix of phi which is known given eta
    t(t)%*%(t*c(exp(t%*%phi))) + inv.cov.mat
  }
  
  ##
  
  get_phi_model <- function(norm.type = "2", grad.tol = 1e-5, 
                            phi.start = rep(1,numFixEff+numUnit*numRanEff), 
                            max.iter = 1e5, inv.cov.mat){
    
    # return the MLE of phi given data and eta
    # norm.type is the type of distance in function "norm"
    # phi.start is our starting guess of phi. Default: vector of ones
    # grad_tol is convergence tolerance of gradient
    # max.iter is maximum number of iterations
    
    invcovmat <- inv.cov.mat
    
    # set phi be our current guess
    x <- phi.start
    
    for (i in 1:max.iter){
      # check local gradient
      j <- get_gradient(phi = x, inv.cov.mat = invcovmat)
      grad.norm <- norm(j, type = norm.type)
      if (grad.norm < grad.tol){
        return(x)
      }
      H <- get_hessian(phi = x, inv.cov.mat = invcovmat)
      inv.H <- chol2inv(chol(H))
      x <- x - inv.H%*%j
    }
    
    return(x)
  }
  
  get_func_for_v <- function(v){
    v1 = v[1]
    v2 = v[2]
    ## Return function of posterior of eta
    .cov.mat = diag(rep(c(fixEffVar, 1/exp(v1), 1/exp(v2)), times = c(numFixEff, numUnit, numUnit)))
    .inv.cov.mat <- chol2inv(chol(.cov.mat))
    .phi.mode <- get_phi_model(inv.cov.mat = .inv.cov.mat)
    .H <- get_hessian(phi = .phi.mode, inv.cov.mat = .inv.cov.mat)
    .la.var <- chol2inv(chol(.H))
    funcval <- get_minus_log_llh(c(.phi.mode), inv.cov.mat = .inv.cov.mat) +
      gammaB * (exp(v1)+exp(v2)) - (numUnit/2+gammaA)*(v1+v2) -
      1/2* log(det(.la.var))
    return(funcval)
  }

  vMode <- optim(c(1,1), fn = get_func_for_v)
  ### Laplace Approximation at etaMode
  vMode <- vMode$par
  # its corresponding covariance matrix
  covMatMode = diag(rep(c(fixEffVar, 1/exp(vMode[1]), 1/exp(vMode[2])), 
                        times = c(numFixEff,numUnit,numUnit)))
  invCovMatMode <- chol2inv(chol(covMatMode))
  # estimation of fixed effects and random effects
  phiModel <- get_phi_model(phi.start = rep(1,numFixEff+numUnit*numRanEff), 
                            inv.cov.mat = invCovMatMode)
  fixEffectModel <- phiModel[1:numFixEff]
  randomEffectModel <- phiModel[numFixEff+1:length(phiModel)]
  
  hessianMode <- get_hessian(phiModel, inv.cov.mat = invCovMatMode)
  la.var.model <- chol2inv(chol(hessianMode))
  
  diag = la.var.model[row(la.var.model) == col(la.var.model)]
  sd_fixEffect = sqrt(diag[1:numFixEff])
  
  results <- list(fix_effect = fixEffectModel, 
                  random_effect = randomEffectModel,
                  variance = 1/exp(vMode), 
                  sd_fixEffect = sd_fixEffect)
  return(results)
}
