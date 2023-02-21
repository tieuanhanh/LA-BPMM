library(microbenchmark)
gammaA = 0.001
gammaB = 0.001
fixEffVar = 10^4
eta = 7

microbenchmark(lapmm(y=seizure.rate, t = t, numUnit=numUnit, obsPerUnit = obsPerUnit, 
                     numFixEff = numFixEff,
                     gammaA = 0.001, gammaB = 0.001, fixEffVar = 10^4))

rjags <- function(){
  jags <- jags.model(file="epilepsy.model.txt",
                     data = my.data,
                     inits = my.inits,
                     n.chains = 3)
  update(jags,1000)
  epi.sim <- coda.samples(model = jags,
                          variable.names = parameters,
                          n.iter=1500, 
                          thin=10)
}

microbenchmark(rjags)

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

get_phi_model <- function(norm.type = "2", grad.tol = 1e-5, 
                          phi.start = rep(1,numFixEff+numUnit), 
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

pre_dat <- function(obsPerUnit, numUnit, y, t){
  z1 = matrix(rep(1,obsPerUnit), ncol=1)
  z = kronecker(diag(numUnit), z1)
  ones = rep(1, length(y))
  t = cbind(ones, t,z)
  t = as.matrix(t)
  return t
}

get_func_for_eta <- function(eta, gamma.a, gamma.b, fix.eff.var){
  ## Return function of posterior of eta 
  .cov.mat = diag(rep(c(fix.eff.var, 1/eta), times = c(numFixEff,numUnit)))
  .inv.cov.mat <- chol2inv(chol(.cov.mat))
  .phi.mode <- get_phi_model(inv.cov.mat = .inv.cov.mat)
  .H <- get_hessian(phi = .phi.mode, inv.cov.mat = .inv.cov.mat)
  .la.var <- chol2inv(chol(.H))
  funcval <- get_minus_log_llh(c(.phi.mode), inv.cov.mat = .inv.cov.mat) + 
    gamma.b * eta - (numUnit/2+gamma.a-1)*log(eta) -1/2* log(det(.la.var))
  
  return(funcval)
}
covMatMode = diag(rep(c(fixEffVar, 1/eta), times = c(numFixEff,numUnit)))
phiModel <- get_phi_model(phi.start = rep(1,numFixEff+numUnit), 
                          inv.cov.mat = invCovMatMode)
hessianMode <- get_hessian(phiModel, inv.cov.mat = invCovMatMode)

microbenchmark(
  pre_dat(obsPerUnit, numUnit, y, t),
  etaMode <- optimize(get_func_for_eta, c(0, 20), tol = 0.0001, 
                      gamma.a = gammaA, gamma.b = gammaB, fix.eff.var = fixEffVar),
  phiModel <- get_phi_model(phi.start = rep(1,numFixEff+numUnit), 
                            inv.cov.mat = invCovMatMode),
  invCovMatMode <- chol2inv(chol(covMatMode)),
  hessianMode <- get_hessian(phiModel, inv.cov.mat = invCovMatMode),
  la.var.model <- chol2inv(chol(hessianMode))
)

######


## Return function of posterior of eta 
cov.mat = diag(rep(c(fixEffVar, 1/eta1, 1/eta2), times = c(numFixEff,numUnit,numUnit)))
inv.cov.mat <- chol2inv(chol(cov.mat))
phi.mode <- get_phi_model(inv.cov.mat = inv.cov.mat)
H <- get_hessian(phi = phi.mode, inv.cov.mat = inv.cov.mat)
la.var <- chol2inv(chol(H))
funcval <- get_minus_log_llh(c(.phi.mode), inv.cov.mat = inv.cov.mat) + 
  gamma.b * eta - (numUnit/2+gamma.a-1)*log(eta) -1/2* log(det(la.var))
etaMode <- optimize(get_func_for_eta, c(0, 20), tol = 0.0001, gamma.a = gammaA, gamma.b = gammaB, fix.eff.var = fixEffVar)

