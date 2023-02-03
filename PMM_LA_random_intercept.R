
# generate 100 units, each units has 10 observations/measurements
rm(list = ls())

generate_data <- function(fix.effect, random.intercept, num.unit, obs.per.unit){
  group = rep(1:num.unit,each=obs.per.unit)
  n = num.unit*obs.per.unit
  num.fix.effect = length(fix.effect)
  set.seed(1)
  
  # Design matrix for the random intercept
  z1 = matrix(rep(1,obs.per.unit), ncol=1)
  z = kronecker(diag(num.unit), z1)
  
  # vector of parameter
  phi = c(fix.effect, random.intercept)
  
  # Design matrix has first column of 1
  x = matrix(nrow = n, ncol = num.fix.effect)
  x[,1] = rep(1, n) 
  
  for (i in 2:num.fix.effect){
    x[,i] = runif(n, min = 0, max = 1) - 0.5
  }
  t = cbind(x, z)
  
  lambda = exp(t%*%phi)
  
  #outcome
  
  y = rpois(n, lambda = lambda) 
  
  data = cbind(group, y, t)
  return(data)
}

fixEffect = c(1,2,3,4,5)
pRecision = 10
numUnit = 100
obsPerUnit = 10
numFixEff = length(fixEffect)
set.seed(2)
randomIntercept = rnorm(numUnit, 0, sd=sqrt(1/pRecision))
pHi = c(fixEffect, randomIntercept)
data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, num.unit = numUnit, obs.per.unit = obsPerUnit)
t = data[,-c(1:2)]
y = data[, 2]
data2 = data.frame(data)

### Spaghetti plot over time
library(ggplot2)

time = rep(1:obsPerUnit, numUnit)
data2 = cbind(data2,time)

### 
ggplot(data = data2, aes(x = factor(time), y = y, color = group)) +       
  geom_line(aes(group = group)) + geom_point()

### Frequentist Poisson Mixed Model 
library(Matrix)
library(lme4)
library(tictoc)

tic("Running time for likelihood method")
fpmm <- glmer(y ~ V4 + V5+ V6+V7+(1|group), data=data2, family="poisson")
toc()
summary(fpmm)

## Get the random intercept
getME(fpmm, "b")

### Laplace Approximation

## Calculate likelihood of data given eta

get_minus_log_llh <- function(phi, inv.cov.mat){
  # x is a vector in R2
  - sum(t %*% phi * y - exp(t%*%phi)) + 1/2 * phi %*% inv.cov.mat %*% phi
}

## Gradient and hessian of f

get_gradient <- function(phi, inv.cov.mat){
  - colSums(t*y - t*c(exp(t%*%phi))) + inv.cov.mat%*%phi
}

get_hessian <- function(phi, inv.cov.mat){
  t(t)%*%(t*c(exp(t%*%phi))) + inv.cov.mat
}

##

get_phi_model <- function(norm.type = "2", grad.tol = 1e-5, phi.start, max.iter = 1e5, cov.mat){
  # norm.type is the type of distance in function "norm"
  # phi.start is our starting guess
  # grad_tol is convergence tolerance of gradient
  # max.iter is maximum number of iterations
  ## Gradient and hessian of f
  
  inv.cov.mat <- chol2inv(chol(cov.mat))
  
  get_gradient <- function(phi, inv.cov.mat){
    - colSums(t*y - t*c(exp(t%*%phi))) + inv.cov.mat%*%phi
  }
  
  get_hessian <- function(phi, inv.cov.mat){
    t(t)%*%(t*c(exp(t%*%phi))) + inv.cov.mat
  }
  
  # set your current guess
  x <- phi.start
  
  for (i in 1:max.iter){
    # check local gradient
    j <- get_gradient(phi = x, inv.cov.mat = inv.cov.mat)
    grad.norm <- norm(j, type = norm.type)
    if (grad.norm < grad.tol){
      return(x)
    }
    H <- get_hessian(phi = x, inv.cov.mat = inv.cov.mat)
    inv.H <- chol2inv(chol(H))
    x <- x - inv.H%*%j
  }
  
  return(x)
}

### Try code with eta = 0.0001
etaFix = 0.001
covMat = diag(rep(c(10^4, 1/etaFix), times = c(numFixEff,numUnit)))
invCovMat <- chol2inv(chol(covMat))
phiModel <- get_phi_model(phi.start = rep(1,numFixEff+numUnit), cov.mat = covMat)
print(phiModel[1:numFixEff])
-get_minus_log_llh(c(phiModel), inv.cov.mat = invCovMat)

### Find the log posterior of eta

gammaA = 0.001
gammaB = 0.001
bandWidthEtaGrid = 0.1
etaGrid = seq(0.1, 20, by = bandWidthEtaGrid)
fixEffVar = 10^4
numFixEff = length(fixEffect)

get_log_posteria_eta <- function(eta.grid, gamma.a, gamma.b, fix.eff.var){
  ## Calculate likelihood of data given eta
  
  .get_minus_log_llh <- function(phi, inv.cov.mat){
    # x is a vector in R2
    - sum(t %*% phi * y - exp(t%*%phi)) + 1/2 * phi %*% inv.cov.mat %*% phi
  }
  
  ## Gradient and hessian of f
  
  .get_gradient <- function(phi, inv.cov.mat){
    - colSums(t*y - t*c(exp(t%*%phi))) + inv.cov.mat%*%phi
  }
  
  .get_hessian <- function(phi, inv.cov.mat){
    t(t)%*%(t*c(exp(t%*%phi))) + inv.cov.mat
  }
  
  .get_phi_model <- function(norm.type = "2", grad.tol = 1e-5, phi.start, max.iter = 1e5, inv.cov.mat){
    # norm.type is the type of distance in function "norm"
    # phi.start is our starting guess
    # grad_tol is convergence tolerance of gradient
    # max.iter is maximum number of iterations
    
    # set your current guess
    x <- phi.start
    
    for (i in 1:max.iter){
      # check local gradient
      j <- .get_gradient(phi = x, inv.cov.mat = inv.cov.mat)
      grad_norm <- norm(j, type = norm.type)
      if (grad_norm < grad.tol){
        return(x)
      }
      H <- .get_hessian(phi = x, inv.cov.mat = inv.cov.mat)
      inv.H <- chol2inv(chol(H))
      x <- x - inv.H%*%j
    }
    
    return(x)
  }
  
  i = 0
  .iVec = rep(1,numFixEff+numUnit)
  
  posterior.eta = 0
  ## Calculate posterior of eta 
  for (eta in eta.grid){
    i <- i+1
    .cov.mat = diag(rep(c(fix.eff.var, 1/eta), times = c(numFixEff,numUnit)))
    .inv.cov.mat <- chol2inv(chol(.cov.mat))
    .phi.mode <- .get_phi_model(phi.start = .iVec, inv.cov.mat = .inv.cov.mat)
    .H <- .get_hessian(phi = .phi.mode, inv.cov.mat = .inv.cov.mat)
    .la.var <- chol2inv(chol(.H))
    posterior.eta[i] <- - .get_minus_log_llh(c(.phi.mode), inv.cov.mat = .inv.cov.mat) - 
      gamma.b * eta + (numUnit/2+gamma.a-1)*log(eta) +1/2* log(det(.la.var))
  }
  return(posterior.eta)
}

tic()
logPosEta = get_log_posteria_eta(eta.grid = etaGrid, gamma.a = gammaA, gamma.b = gammaB, fix.eff.var = fixEffVar)
toc()

### Find the eta with maximum posterior

maxLogPosEta <- max(logPosEta)
maxLogPosEta
etaMode <- etaGrid[which.max(logPosEta)]
etaMode

### Plot eta posterior

plot(etaGrid, logPosEta, cex=0.1, col = "red", ylab = "P(eta/D)", 
     main = "Log Posterior distribution of eta grid", type="l")
abline(v = bandWidthEtaGrid* which.max(logPosEta), lty=2, col = "blue")

### LA at eta mode

covMatMode = diag(rep(c(fixEffVar, 1/etaMode), times = c(numFixEff,numUnit)))
invCovMatMode <- chol2inv(chol(covMatMode))
phiModel <- get_phi_model(phi.start = rep(1,numFixEff+numUnit), 
                             cov.mat = covMatMode)
phiModel[1:numFixEff]
hessianMode <- get_hessian(phiModel, inv.cov.mat = invCovMatMode)
laVarMode <- chol2inv(chol(hessianMode))

### Inference for fix intercept 
get_inference_fix_effect <- function(phi.model, la.var.model, i){
  beta.mode = phi.model[i]
  beta.var = la.var.model[i,i]
  x.grid = seq(beta.mode-3,beta.mode+3,by=0.001)
  beta.posterior = dnorm(x.grid, beta.mode, sd=sqrt(beta.var))
  lw.bound.hpd = beta.mode - qnorm(0.975)*sqrt(beta.var)
  up.bound.hpd = beta.mode + qnorm(0.975)*sqrt(beta.var)
  
  print("95 CI for LA of marginal posterior beta0")
  print(lw.bound.hpd)
  print(beta.mode)
  print(up.bound.hpd)
  
  par(mfrow=c(1,1))
  plot(x.grid, beta.posterior, col = "black", type = "l", lty = 1, xlab = "beta",
       ylab="P(beta/eta.mode,D)", main = "LA for Marginal Posterior of beta")
  abline(v=c(lw.bound.hpd, beta.mode, up.bound.hpd), col = c("blue", "red", "blue"), lty= c(3,2,3))
  legend("topright", legend=c("lower bound CI", "true beta0", "upper bound CI"), 
         col = c("blue", "red", "blue"), lty= c(3,2,3), cex = 0.6)
}
get_inference_fix_effect(phi.model = phiModel, la.var.model = laVarMode, i=1)

