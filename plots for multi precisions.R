library(Matrix)
library(lme4)
library(boot)
library(INLA)
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(INLA)
library(coda)
library(rjags)
library(readr)
library(dplyr)
library(runjags)

### Plots for multivariate \eta

fixEffect = c(0.5,1,0.5, 1.5)
pRecision = c(1,2)
numUnit = 50
obsPerUnit = 6
numFixEff = length(fixEffect)
set.seed(2)
randomIntercept = rnorm(numUnit, 0, sd=sqrt(1/pRecision[1]))
randomSlope = rnorm(numUnit, 0, sd=sqrt(1/pRecision[2]))

generate_data <- function(fix.effect, random.intercept, random.slope, num.unit, obs.per.unit){
  group = rep(1:num.unit,each=obs.per.unit)
  n = num.unit*obs.per.unit
  num.fix.effect = length(fix.effect)
  
  # Design matrix for the random intercept
  z1 = matrix(rep(1,obs.per.unit), ncol=1)
  z = kronecker(diag(num.unit), z1)
  
  # vector of parameter
  phi = c(fix.effect, random.intercept, random.slope)
  
  # Design matrix has first column of 1
  x = matrix(nrow = n, ncol = num.fix.effect)
  x[,1] = rep(1, n) 
  x[,2] = runif(n, min = 0, max = 1) - 0.5
  x[,3] = rnorm(n, 0, 1)
  x[,4] = rbinom(n, 1, 0.5) - 0.5
  
  r = rep(1:obsPerUnit/obsPerUnit, numUnit) - 3.5/6
  r = matrix(r)
  u = matrix(0, n, numUnit)
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
  # if y is too large then phi model diverge => should remove outliers of y.
  # In this case, we set constraint for y
  # for (i in 1:n){
  #   y[i] = ifelse(y[i] < 5000, y[i], 5000)
  # }
  # 
  data = cbind(group, y, x[,2:4], r)
  return(data)
}
set.seed(2)
data = generate_data(fix.effect = fixEffect, 
                     random.intercept = randomIntercept,
                     random.slope = randomSlope,
                     num.unit = numUnit, obs.per.unit = obsPerUnit)

t = data[,c(3:5)]
y = data[, 2]
r = data[, 6]

gammaA = 0.01
gammaB = 0.01
fixEffVar = 10^4

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
vMode1 <- vMode[1]
vMode2 <- vMode[2]

# Approach 2
# its corresponding covariance matrix
covMatMode = diag(rep(c(fixEffVar, 1/exp(vMode1), 1/exp(vMode2)), times = c(numFixEff,numUnit,numUnit)))
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

######################
# INLA
data2 = data.frame(data)
colnames(data2)[1:6] = c("group", "y", "RV1", "RV2", "RV3", "RV4")
data2$slopeid = data2$group + numUnit

prior.fixed <- list(mean.intercept = 0, prec.intercept = 0.0001,
                    mean = 0, prec = 0.0001)
prec.prior <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01)))

formula = y ~ RV1 + RV2 + RV3 + 
  f(group, model ="iid", hyper = prec.prior) + 
  f(slopeid, RV4, model ="iid")

model.inla <- inla(formula,
                   data = data2, family = "poisson",
                   control.fixed = prior.fixed)
###MCMC
cat("model
 {
    for (j in 1:J){
      b0[j] ~ dnorm(0, tau1)
      b1[j] ~ dnorm(0, tau2)
    }
    for (i in 1:N){
      rate[i] ~ dpois(phi[i])
      log(phi[i]) <- beta0 + beta1*RV1[i] + beta2*RV2[i] + beta3*RV3[i]
                  + b0[group[i]] + b1[group[i]]*RV4[i]
    }
    tau1 ~ dgamma(0.01, 0.01)
    tau2 ~ dgamma(0.01, 0.01)
    var1 = 1/tau1
    var2 = 1/tau2
    beta0 ~ dnorm(0,1.0E-4)
    beta1 ~ dnorm(0,1.0E-4)
    beta2 ~ dnorm(0,1.0E-4)
    beta3 ~ dnorm(0,1.0E-4)
  }", file="simulation3.model.txt")
t2 = data.frame(data[,-c(1:2)])
colnames(t2)[1:4] = c("RV1", "RV2", "RV3", "RV4")
y = data[, 2]
group = data[,1]
N <- length(y)
J <- n_distinct(group)

my.data <- list(rate=y, RV1 = t2$RV1, RV2 = t2$RV2, 
                RV3 = t2$RV3, RV4 = t2$RV4,
                N = N, J = J, group = group)
# set seed for sample chains
my.inits <- list(
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 111, beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, tau1=1, tau2=1),
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 112, beta0=0.1,beta1=0.1, beta2=0.1,beta3=0.1, tau1=0.5, tau2=1),
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 113, beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, tau1=0.5, tau2=0.5)
)

parameters <- c("beta0","beta1","beta2","beta3", "tau1", "tau2", "var1", "var2")

jags <- jags.model(file="simulation3.model.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 3)
update(jags,1500)
sim.sim <- coda.samples(model = jags,
                        variable.names = parameters,
                        n.iter=5000, 
                        thin=10)
# combine 3 chains
sim.mcmc <- as.mcmc.list(sim.sim)
sim.combined <- combine.mcmc(sim.mcmc)

### Likelihood
fpmm <- glmer(y ~ RV1 + RV2 + RV3 +(1+RV4|group), data=data2, family="poisson")

### compute CI of variance 

grid_n_plot = 50
width_plot = 0.1
i = 0
v1_grid_plot = seq(vMode1 - grid_n_plot/2*width_plot, 
                   vMode1 + grid_n_plot/2*width_plot, width_plot)
log_post_v1 = rep(0, length(v1_grid_plot))
for (k in v1_grid_plot){
  i = i + 1
  log_post_v1[i] = - get_func_for_v(c(k, vMode2))
}
log_post_v1_scale = log_post_v1 - max(log_post_v1)
post_v1 = exp(log_post_v1_scale)
post_v1_norm_plot = post_v1/sum(post_v1)

grid_n = 6
width = 0.25
i = 0
v1_grid = seq(vMode1 - grid_n/2*width, 
              vMode1 + grid_n/2*width, width)
log_post_v1 = rep(0, length(v1_grid))
for (k in v1_grid){
  i = i + 1
  log_post_v1[i] = - get_func_for_v(c(k, vMode2))
}
log_post_v1_scale = log_post_v1 - max(log_post_v1)
post_v1 = exp(log_post_v1_scale)
post_v1_norm = post_v1/sum(post_v1)
mean_v1 = sum(v1_grid*post_v1_norm)
m2_v1 = sum(v1_grid^2*post_v1_norm)
var_v1 = m2_v1 - mean_v1^2
sd_v1 = sqrt(var_v1)

lb_95_v1 = mean_v1 - 1.96*sd_v1
ub_95_v1 = mean_v1 + 1.96*sd_v1
lb_95_var1 = 1/exp(mean_v1 + 1.96*sd_v1)
ub_95_var1 = 1/exp(mean_v1 - 1.96*sd_v1)
ci_95_var1 = c(lb_95_var1, ub_95_var1)

# compute statistics of variance of random intercept (mean, sd)
mean_var1 = sum((1/exp(v1_grid))*post_v1_norm)

# plot posterior of nu 1
# png("/Users/nana/rgit/LA-BPMM/pic/post nu1 (Methods)")
par(mar = c(5,5,3,3))
plot(v1_grid_plot, post_v1_norm_plot, col = "blue", type = "l",
     xlim = c(-1.2,1),
     xlab = expression(paste(nu[1])),
     ylab = expression(paste(hat(p), "(", nu[1], "| D)")))
points(x = v1_grid, y = rep(-0.0003, length(v1_grid)), 
       col = "red", pch = 17, cex = 0.5)
abline(v=c(mean_v1, vMode1, lb_95_v1, ub_95_v1), col = c("blue", "red", "black", "black"))
# dev.off()

# plot posterior of variance 1
# png("/Users/nana/rgit/LA-BPMM/pic/post var2 (Methods)")
par(mar = c(5,5,3,3))
plot(1/exp(v1_grid_plot), post_v1_norm_plot, col = "blue", 
     type = "l", xlim = c(0,3),
     xlab = expression(sigma[0]^2),
     ylab = expression(paste(hat(p),"(", sigma[0]^2, " | D)")))
abline(v=c(mean_var1, 1/exp(vMode1), ci_95_var1), 
       col = c("blue", "red", "black", "black"))
# dev.off()

### the second variance
grid_n_plot = 50
width_plot = 0.1
i = 0
v2_grid_plot = seq(vMode2 - grid_n_plot/2*width_plot, 
              vMode2 + grid_n_plot/2*width_plot, width_plot)
log_post_v2 = rep(0, length(v2_grid_plot))
for (k in v2_grid_plot){
  i = i + 1
  log_post_v2[i] = - get_func_for_v(c(vMode1, k))
}
log_post_v2_scale = log_post_v2 - max(log_post_v2)
post_v2 = exp(log_post_v2_scale)
post_v2_norm_plot = post_v2/sum(post_v2)

grid_n = 6
width = 0.3
i = 0
v2_grid = seq(vMode2 - grid_n/2*width, 
              vMode2 + grid_n/2*width, width)
log_post_v2 = rep(0, length(v2_grid))
for (k in v2_grid){
  i = i + 1
  log_post_v2[i] = - get_func_for_v(c(vMode1, k))
}
log_post_v2_scale = log_post_v2 - max(log_post_v2)
post_v2 = exp(log_post_v2_scale)
post_v2_norm = post_v2/sum(post_v2)
mean_v2 = sum(v2_grid*post_v2_norm)
m2_v2 = sum(v2_grid^2*post_v2_norm)
var_v2 = m2_v2 - mean_v2^2
sd_v2 = sqrt(var_v2)

ub_95_v2 = mean_v2 + 1.96*sd_v2
lb_95_v2 = mean_v2 - 1.96*sd_v2
lb_95_var2 = 1/exp(mean_v2 + 1.96*sd_v2)
ub_95_var2 = 1/exp(mean_v2 - 1.96*sd_v2)
ci_95_var2 = c(lb_95_var2, ub_95_var2)


# statistics of variance 2
mean_var2 = sum((1/exp(v2_grid))*post_v2_norm)

# plot posterior of nu 2
# png("/Users/nana/rgit/LA-BPMM/pic/post nu2 (Methods)")
par(mar = c(5,5,3,3))
plot(v2_grid_plot, post_v2_norm_plot, col = "blue", type = "l",
     xlim = c(-1, 1.7),
     xlab = expression(paste(nu[2])),
     ylab = expression(paste(hat(p),"(", nu[2], "| D)")))
points(x = v2_grid, y = rep(-0.0003, length(v2_grid)), 
       col = "red", pch = 17, cex = 0.5)
abline(v=c(mean_v2, vMode2, lb_95_v2, ub_95_v2), col = c("blue", "red", "black", "black"))
# dev.off()

# plot posterior of variance 2
# png("/Users/nana/rgit/LA-BPMM/pic/post var2 (Methods)")
par(mar = c(5,5,3,3))
plot(1/exp(v2_grid_plot), post_v2_norm_plot, col = "blue", 
     type = "l", xlim = c(0,2.5),
     xlab = expression(sigma[1]^2),
     ylab = expression(paste(hat(p), "(", sigma[1]^2, " | D)")))
abline(v=c(mean_var2, 1/exp(vMode2), ci_95_var2), 
       col = c("blue", "red", "black", "black"))
# dev.off()

# plot posterior of fixed effects

phi.mode = function(v1, v2){
  cov.mat = diag(rep(c(fixEffVar, 1/exp(v1), 1/exp(v2)), 
                     times = c(numFixEff, numUnit, numUnit)))
  inv.cov.mat <- chol2inv(chol(cov.mat))
  phi.mode <- get_phi_model(inv.cov.mat = inv.cov.mat)
  return(phi.mode)
}

var_la = function(v1, v2){
  cov.mat = diag(rep(c(fixEffVar, 1/exp(v1), 1/exp(v2)), 
                     times = c(numFixEff, numUnit, numUnit)))
  inv.cov.mat <- chol2inv(chol(cov.mat))
  phi.mode <- get_phi_model(inv.cov.mat = inv.cov.mat)
  hessianMode <- get_hessian(phi = phi.mode, 
                             inv.cov.mat = inv.cov.mat)
  var_la <- chol2inv(chol(hessianMode))
  return (var_la)
}

# plots for fixed effects
width = 0.3
v1_grid = seq(vMode1 - grid_n/2*width, 
              vMode1 + grid_n/2*width, width)
width = 0.3
v2_grid = seq(vMode2 - grid_n/2*width, 
              vMode2 + grid_n/2*width, width)
phi_grid = matrix(nrow = length(v1_grid)*length(v2_grid),
                  ncol = numFixEff)
var_grid = matrix(nrow = length(v1_grid)*length(v2_grid),
                  ncol = numFixEff)
log_weights = rep(0, length(v1_grid)*length(v2_grid))

i= 0
for (v1 in v1_grid){
  for (v2 in v2_grid){
  i = i+1
  phi_vector = phi.mode(v1,v2)
  var_vector = var_la(v1,v2)
  log_weights[i] = - get_func_for_v(c(v1, v2))
  for (j in 1:numFixEff){
    phi_grid[i, j] = phi_vector[j]
    var_grid[i, j] <- var_vector[j,j]
    }
  }
}
log_weights_scale = log_weights - max(log_weights)
weights_scale = exp(log_weights_scale)
weights = weights_scale/sum(weights_scale)

mean_phi <- colSums(phi_grid*weights)
mean_var <- colSums(var_grid*weights)

xaxis = seq(0, 2, by=0.001)

# png("/Users/nana/rgit/LA-BPMM/pic/LABPMM grid2 beta0 (Methods).png",
#     width = 700, height = 400)
par(mar = c(5,5,3,3))
plot(xaxis, dnorm(xaxis, mean_phi[1], sqrt(mean_var[1])), col = "blue",
     type = "l",
     xlim = c(0,1.5),
     ylim = c(0,4),
     xlab = expression(beta[0]),
     ylab = expression(paste(hat(p),"(", beta[0], " | D)")))
for (i in 1:length(log_weights)){
  lines(xaxis, dnorm(xaxis, phi_grid[i,1], sqrt(var_grid[i,1])), col = "red", type = "l", lty = 2)
  lines(xaxis, dnorm(xaxis, phi_grid[i,1], sqrt(var_grid[i,1]))*weights[i], col = "green", type = "l")
}
legend("topright", legend=c("LABPMM-grid",
                            expression(paste("unweighted ", hat(p),"(", beta[0], " |", nu^{(m)},", D)")),
                            expression(paste("weighted ", hat(p), "(", beta[0], " |", nu^{(m)},", D)"))),
       col=c("blue", "red", "green"), lty=1:2, cex=0.8)
# dev.off()

plot(xaxis, dnorm(xaxis, mean_phi[2], sqrt(mean_var[2])), col = "blue",
     type = "l",
     xlim = c(0.5,1.5),
     ylim = c(0,4),
     xlab = expression(beta[1]),
     ylab = expression(paste(hat(p),"(", beta[1], " | D)")))
for (i in 1:length(log_weights)){
  lines(xaxis, dnorm(xaxis, phi_grid[i,2], sqrt(var_grid[i,2])), col = "red", type = "l", lty = 2)
  lines(xaxis, dnorm(xaxis, phi_grid[i,2], sqrt(var_grid[i,2]))*weights[i], col = "green", type = "l")
}
legend("topright", legend=c("LABPMM-grid",
                            expression(paste("unweighted ", hat(p),"(", beta[1], " |", nu^{(m)},", D)")),
                            expression(paste("weighted ", hat(p), "(", beta[1], " |", nu^{(m)},", D)"))),
       col=c("blue", "red", "green"), lty=1:2, cex=0.8)


plot(xaxis, dnorm(xaxis, mean_phi[3], sqrt(mean_var[3])), col = "blue",
     type = "l",
     xlim = c(0.2,0.8),
     ylim = c(0,13),
     xlab = expression(beta[2]),
     ylab = expression(paste(hat(p),"(", beta[2], " | D)")))
for (i in 1:length(log_weights)){
  lines(xaxis, dnorm(xaxis, phi_grid[i,3], sqrt(var_grid[i,3])), col = "red", type = "l", lty = 2)
  lines(xaxis, dnorm(xaxis, phi_grid[i,3], sqrt(var_grid[i,3]))*weights[i], col = "green", type = "l")
}
legend("topright", legend=c("LABPMM-grid",
                            expression(paste("unweighted ", hat(p),"(", beta[2], " |", nu^{(m)},", D)")),
                            expression(paste("weighted ", hat(p), "(", beta[2], " |", nu^{(m)},", D)"))),
       col=c("blue", "red", "green"), lty=1:2, cex=0.8)


# png("/Users/nana/rgit/LA-BPMM/pic/LABPMM grid2 beta3 (Methods).png",
#     width = 700, height = 400)

plot(xaxis, dnorm(xaxis, mean_phi[4], sqrt(mean_var[4])), col = "blue",
     type = "l",
     xlim = c(1,1.8),
     xlab = expression(beta[3]),
     ylab = expression(paste(hat(p),"(", beta[3], " | D)")))
for (i in 1:length(log_weights)){
  lines(xaxis, dnorm(xaxis, phi_grid[i,4], sqrt(var_grid[i,4])), col = "red", type = "l", lty = 2)
  lines(xaxis, dnorm(xaxis, phi_grid[i,4], sqrt(var_grid[i,4]))*weights[i], col = "green", type = "l")
}
legend("topright", legend=c("LABPMM-grid",
                            expression(paste("unweighted ", hat(p),"(", beta[3], " |", nu^{(m)},", D)")),
                            expression(paste("weighted ", hat(p), "(", beta[3], " |", nu^{(m)},", D)"))),
       col=c("blue", "red", "green"), lty=1:2, cex=0.8)
# dev.off()

phi_est = phi.mode(vMode[1], vMode[2])[1:numFixEff]
var_mat <- var_la(vMode[1], vMode[2])[c(1:numFixEff),c(1:numFixEff)]
var_est <- var_mat[row(var_mat) == col(var_mat)]

### compare only LABPMM
# for (i in 1:numFixEff){
#   xaxis = seq(fixEffect[i]-0.5, fixEffect[i]+0.7, by=0.001)
#   plot(xaxis, dnorm(xaxis, mean_phi[i], sqrt(mean_var[i])), col = "blue",
#        type = "l", 
#        xlab = expression(beta[i]),
#        ylab = expression(paste(hat(p),"(", beta[i], " | D)")))
#   lines(xaxis, dnorm(xaxis, phi_est[i], sqrt(var_est[i])), col = "red", type = "l", lty = 2)
#   legend("topright", legend=c("LABPMM-grid", 
#                               "LABPMM-MAP"),
#          col=c("blue", "red"), lty=c(1,2), cex=0.8)
#   abline(v = c(mean_phi[i], phi_est[i]),col=c("blue", "red"), lty=1:2) 
#   
# }

plot(xaxis, dnorm(xaxis, mean_phi[1], sqrt(mean_var[1])), col = "blue",
     type = "l", 
     xlim = c(0,1.5),
     xlab = expression(beta[0]),
     ylab = expression(paste(hat(p),"(", beta[0], " | D)")))
lines(xaxis, dnorm(xaxis, phi_est[1], sqrt(var_est[1])), col = "red", type = "l", lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP"),
       col=c("blue", "red"), lty=c(1,2), cex=0.8)
abline(v = c(mean_phi[1], phi_est[1]),col=c("blue", "red"), lty=1:2) 

plot(xaxis, dnorm(xaxis, mean_phi[2], sqrt(mean_var[2])), col = "blue",
     type = "l", 
     xlim = c(0.5,1.4),
     xlab = expression(beta[1]),
     ylab = expression(paste(hat(p),"(", beta[1], " | D)")))
lines(xaxis, dnorm(xaxis, phi_est[2], sqrt(var_est[2])), col = "red", type = "l", lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP"),
       col=c("blue", "red"), lty=c(1,2), cex=0.8)
abline(v = c(mean_phi[2], phi_est[2]),col=c("blue", "red"), lty=1:2) 


plot(xaxis, dnorm(xaxis, mean_phi[3], sqrt(mean_var[3])), col = "blue",
     type = "l", 
     xlim = c(0.3,0.7),
     xlab = expression(beta[2]),
     ylab = expression(paste(hat(p),"(", beta[2], " | D)")))
lines(xaxis, dnorm(xaxis, phi_est[3], sqrt(var_est[3])), col = "red", type = "l", lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP"),
       col=c("blue", "red"), lty=c(1,2), cex=0.8)
abline(v = c(mean_phi[3], phi_est[3]),col=c("blue", "red"), lty=1:2) 


plot(xaxis, dnorm(xaxis, mean_phi[4], sqrt(mean_var[4])), col = "blue",
     type = "l", 
     xlim = c(1,1.8),
     xlab = expression(beta[3]),
     ylab = expression(paste(hat(p),"(", beta[3], " | D)")))
lines(xaxis, dnorm(xaxis, phi_est[4], sqrt(var_est[4])), col = "red", type = "l", lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP"),
       col=c("blue", "red"), lty=c(1,2), cex=0.8)
abline(v = c(mean_phi[4], phi_est[4]),col=c("blue", "red"), lty=1:2) 


# compare Bayesian

intercept.inla <- model.inla$marginals.fixed[[1]]
plot(sim.combined[,1], trace = FALSE, density=TRUE,
     type = "l", 
     xlab = expression(beta[0]),
     ylab = expression(paste(hat(p),"(", beta[0], " | D)")),
     main = "")
lines(xaxis, dnorm(xaxis, mean_phi[1], sqrt(mean_var[1])), 
      col = "blue",
      type = "l")
lines(xaxis, dnorm(xaxis, phi_est[1], sqrt(var_est[1])), 
      col = "red", type = "l", lty = 2)
lines(intercept.inla, col = 'black', type = 'l', lty = 2)
abline(v = fixef(fpmm)[1], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "black", "black"), lty=c(1,2,2,1), cex=0.8)

intercept.inla <- model.inla$marginals.fixed[[2]]
plot(sim.combined[,2], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0,4),
     xlab = expression(beta[1]),
     ylab = expression(paste(hat(p),"(", beta[1], " | D)")),
     main = "")
lines(xaxis, dnorm(xaxis, mean_phi[2], sqrt(mean_var[2])), 
      col = "blue",
      type = "l")
lines(xaxis, dnorm(xaxis, phi_est[2], sqrt(var_est[2])), 
      col = "red", type = "l", lty = 2)
lines(intercept.inla, col = 'black', type = 'l', lty = 2)
abline(v = fixef(fpmm)[2], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "black", "black"), lty=c(1,2,2,1), cex=0.8)


intercept.inla <- model.inla$marginals.fixed[[3]]
plot(sim.combined[,3], trace = FALSE, density=TRUE,
     type = "l", 
     xlab = expression(beta[2]),
     ylab = expression(paste(hat(p),"(", beta[2], " | D)")),
     main = "")
lines(xaxis, dnorm(xaxis, mean_phi[3], sqrt(mean_var[3])), 
      col = "blue",
      type = "l")
lines(xaxis, dnorm(xaxis, phi_est[3], sqrt(var_est[3])), 
      col = "red", type = "l", lty = 2)
lines(intercept.inla, col = 'black', type = 'l', lty = 2)
abline(v = fixef(fpmm)[3], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "black", "black"), lty=c(1,2,2,1), cex=0.8)


intercept.inla <- model.inla$marginals.fixed[[4]]
plot(sim.combined[,4], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0,5.5),
     xlab = expression(beta[3]),
     ylab = expression(paste(hat(p),"(", beta[3], " | D)")),
     main = "")
lines(xaxis, dnorm(xaxis, mean_phi[4], sqrt(mean_var[4])), 
      col = "blue",
      type = "l")
lines(xaxis, dnorm(xaxis, phi_est[4], sqrt(var_est[4])), 
      col = "red", type = "l", lty = 2)
lines(intercept.inla, col = 'black', type = 'l', lty = 2)
abline(v = fixef(fpmm)[4], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "black", "black"), lty=c(1,2,2,1), cex=0.8)

# posterior of sigma
inv <- function(x) 1/x
prec1 <- model.inla$marginals.hyperpar$`Precision for group`
marg.var1 <- inla.tmarginal(inv, prec1)

xaxis = seq(-1, 3, by=0.001)
plot(sim.combined[,7], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0, 1.7),
     xlim = c(0.5, 2.5),
     xlab = expression(sigma[0]^2),
     ylab = expression(paste(hat(p),"(", sigma[0]^2, " | D)")),
     main = "")
set.seed(12)
y4 <- 1/exp(rnorm(xaxis, mean_v1, sd_v1))
lines(density(y4),
      col = "red", type = "l", lty = 2)
lines(marg.var1, col = 'green', type = 'l', lty = 2)
abline(v = as.data.frame(VarCorr(fpmm))$vcov[1], lty = 2)
legend("topright", legend=c("LABPMM",
                            "INLA",
                            "MCMC"),
       col=c("red", "green", "black"), lty=c(2,2,1), cex=0.8)

prec2 <- model.inla$marginals.hyperpar$`Precision for slopeid`
marg.var2 <- inla.tmarginal(inv, prec2)
plot(sim.combined[,8], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0, 1.7),
     xlim = c(0, 2.5),
     xlab = expression(sigma[1]^2),
     ylab = expression(paste(hat(p),"(", sigma[1]^2, " | D)")),
     main = "")
set.seed(12)
y4 <- 1/exp(rnorm(xaxis, mean_v2, sd_v2))
lines(density(y4),
      col = "red", type = "l", lty = 2)
lines(marg.var2, col = 'green', type = 'l', lty = 2)
abline(v = as.data.frame(VarCorr(fpmm))$vcov[2], lty = 2)
legend("topright", legend=c("LABPMM",
                            "INLA",
                            "MCMC"),
       col=c("red", "green", "black"), lty=c(2,2,1), cex=0.8)

