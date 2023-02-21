lapmm2 <- function(y, t, r, numUnit, obsPerUnit, 
                  gammaA = 0.001, gammaB = 0.001, fixEffVar = 10^4){
  
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
  
  get_func_for_eta <- function(eta){
    eta1 = eta[1]
    eta2 = eta[2]
    ## Return function of posterior of eta 
    .cov.mat = diag(rep(c(fixEffVar, 1/eta1, 1/eta2), times = c(numFixEff, numUnit, numUnit)))
    .inv.cov.mat <- chol2inv(chol(.cov.mat))
    .phi.mode <- get_phi_model(inv.cov.mat = .inv.cov.mat)
    .H <- get_hessian(phi = .phi.mode, inv.cov.mat = .inv.cov.mat)
    .la.var <- chol2inv(chol(.H))
    funcval <- - get_minus_log_llh(c(.phi.mode), inv.cov.mat = .inv.cov.mat) + 
      gammaB * eta1 - (numUnit/2+gammaA-1)*log(eta1) - 1/2* log(det(.la.var)) -
      (numUnit/2+gammaA-1)*log(eta2) + gammaB * eta2
    return(funcval)
  }
  
  # etaMode: the eta at which the log posterior of eta is maximized
  #etaMode <- optim( c(1,1), get_func_for_eta, NULL, method = "BFGS")
  etaMode <- optim( c(1,1), get_func_for_eta, NULL, method = "L-BFGS-B", 
                    lower = rep(0,2), upper = c(20,20))
  ### Laplace Approximation at etaMode
  etaMode <- etaMode$par
  etaMode1 <- etaMode[1]
  etaMode2 <- etaMode[2]
  # its corresponding covariance matrix
  covMatMode = diag(rep(c(fixEffVar, 1/etaMode1, 1/etaMode2), times = c(numFixEff,numUnit,numUnit)))
  invCovMatMode <- chol2inv(chol(covMatMode))
  # estimation of fixed effects and random effects
  phiModel <- get_phi_model(phi.start = rep(1,numFixEff+numUnit*numRanEff), 
                            inv.cov.mat = invCovMatMode)
  fixEffectModel <- phiModel[1:numFixEff]
  randomEffectModel <- phiModel[numFixEff+1:numUnit]
  
  hessianMode <- get_hessian(phiModel, inv.cov.mat = invCovMatMode)
  la.var.model <- chol2inv(chol(hessianMode))
  
  results <- list(fix_effect = fixEffectModel, random_effect = randomEffectModel,
                  etaMode = etaMode, la.var.model = la.var.model)
  
  return(results)
}

#example 1: simulation data

generate_data <- function(fix.effect, random.intercept, random.slope, num.unit, obs.per.unit){
  group = rep(1:num.unit,each=obs.per.unit)
  n = num.unit*obs.per.unit
  num.fix.effect = length(fix.effect)
  set.seed(1)
  
  # Design matrix for the random intercept
  z1 = matrix(rep(1,obs.per.unit), ncol=1)
  z = kronecker(diag(num.unit), z1)
  
  # vector of parameter
  phi = c(fix.effect, random.intercept, random.slope)
  
  # Design matrix has first column of 1
  x = matrix(nrow = n, ncol = num.fix.effect)
  x[,1] = rep(1, n) 
  
  for (i in 2:num.fix.effect){
    x[,i] = runif(n, min = 0, max = 1) - 0.5
  }
  
  r = runif(n, min = 0, max = 1) - 0.5
  r = matrix(r)
  u = matrix(0, nrow(r), numUnit)
  j1 = 1
  for (i in 1:numUnit){
    j2 = j1 + obsPerUnit - 1
    u[j1:j2, i] = r[j1:j2]
    j1 = j2 + 1
  }
  
  t = cbind(x, z, u)
  
  lambda = exp(t%*%phi)
  
  #outcome
  
  y = rpois(n, lambda = lambda) 
  
  data = cbind(group, y, x, r)
  return(data)
}

fixEffect = c(1,2,3,4,5)
pRecision = c(1,2)
numUnit = 60
obsPerUnit = 4
numFixEff = length(fixEffect)
set.seed(2)
randomIntercept = rnorm(numUnit, 0, sd=sqrt(1/pRecision[1]))
randomSlope = rnorm(numUnit, 0, sd=sqrt(1/pRecision[2]))

data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, 
                     random.slope = randomSlope,
                     num.unit = numUnit, obs.per.unit = obsPerUnit)
t = data[,c(4:7)]
y = data[, 2]
r = data[, 8]

model1 <- lapmm2(y=y, t=t, r = r, numUnit=numUnit, obsPerUnit=obsPerUnit, 
                gammaA = 0.5, gammaB = 0.0164, fixEffVar = 10^4)
model1$fix_effect
model1$etaMode

mat = model1$la.var.model
diag = mat[row(mat) == col(mat)]
sqrt(diag[1:5])

# Model 2: Episilep
library(HSAUR3)

#load data
data = epilepsy
attach(epilepsy)

n = nrow(data)
numUnit = 59
obsPerUnit = 4
numFixEff = 6

x2 = log(base/4)
x2 = x2-mean(x2)
x3 = ifelse(treatment=="placebo", 0, 1)
x4 = log(base/4)*x3
x4 = x4-mean(x4)
x5 = log(age)
x5 = (x5-mean(x5))
x6 = (as.numeric(period)-2.5)*2/10

y=seizure.rate
t = cbind(x2,x3,x4,x5,x6)

model2 <- lapmm2(y=seizure.rate, t = t, r=x6, numUnit=numUnit, obsPerUnit = obsPerUnit, 
                gammaA = 0.001, gammaB = 0.001, fixEffVar = 10^4)
model2$fix_effect
model2$etaMode

mat = model2$la.var.model
diag = mat[row(mat) == col(mat)]
sqrt(diag[1:6])
