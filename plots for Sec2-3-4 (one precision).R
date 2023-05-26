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

### figure to illustrate the laplace approximation (Fig 2.1)

# generate n data point, each ~ Poi (lambda = 2)
n = 50
true.mu = 2
set.seed(1)
d = rpois(n = n, lambda = true.mu)

prior.shape <- 0.01
prior.rate <- 0.01
pmode <- (sum(d) + prior.shape - 1)/(length(d) + prior.rate)
pvar <- (sum(d) + prior.shape - 1)/((length(d) + prior.rate)^2)
exact_mean <- (sum(d) + prior.shape)/(length(d) + prior.rate)
curve(dgamma(x, shape = sum(d) + prior.shape, rate = length(d) + prior.rate),
      xlim = c(0,4), n = 1000,
      lwd = 3, xlab = expression(mu), ylab = expression(paste("p(", mu, " | D)")))
curve(dgamma(x, shape = prior.shape, rate = prior.rate), add = TRUE,
      lty = 2)
curve(dnorm(x, mean = pmode, sd = sqrt(pvar)), 
     n = 1000, add = TRUE, col = 2, lwd = 2)
legend("topright", 
       legend = c("Posterior Density", "Prior Density", "Laplace Approx"), 
       lty = c(1, 2, 1), lwd = c(3, 1, 1), col = c(1, 1, 2), bty = "n")
ci_la <- c(qnorm(0.025, mean = pmode, sd = sqrt(pvar)), 
          qnorm(0.975, mean = pmode, sd = sqrt(pvar)))
ci_conju <- c(qgamma(0.025, shape = sum(d) + prior.shape, rate = length(d) + prior.rate),
              qgamma(0.975, shape = sum(d) + prior.shape, rate = length(d) + prior.rate))
abline(v = c(exact_mean, pmode), col = c(1,2))

### example nested approximation (Figure 2.2.a)

v = 2
a = 1
b = 1

mode_mu <- function(del){
  (sum(d) + v/2 - 1)/(length(d) + v*del/2)
}

var_la <- function(del){
  (sum(d) + v/2 - 1)/(length(d) + v*del/2)^2
}

log_post_del <- function(del){
  result <- (v/2 + a -1)*log(del) - b*del - (sum(d) + v/2)*log(n + v*del/2)
  return(result)
}

mode_del <- optimize(log_post_del, c(0, 100), tol = 0.0001,
                     maximum=TRUE)
del_est <- mode_del$maximum

mu_est <- mode_mu(del_est)
var_est <- var_la(del_est)

del_range <- seq(0.01, 2.5, by = 0.01)
pos_del <- exp(log_post_del(del_range) - max(log_post_del(del_range)))
norm_pos_del <- pos_del/(sum(pos_del))
# Figure 2.2b
plot(x = del_range, y = norm_pos_del, type = "l", 
     xlab = expression(delta),
     ylab = expression(paste("p(", delta, " | D)")))

del_grid <- seq(0.03, 1.5, by = 0.15)
points(x = del_grid, y = rep(-0.0003, length(del_grid)), col = "red", pch = 17)
# abline(v = del_grid)

mu = seq(1.3, 2.8, by=0.001)
mu_grid <- mode_mu(del_grid)
var_grid <- var_la(del_grid)
pos_grid <- exp(log_post_del(del_grid)-min(log_post_del(del_grid)))
weight_grid <- pos_grid/sum(pos_grid)
mean_mu <- sum(mu_grid*weight_grid)
mean_var <- sum(var_grid*weight_grid)
ci1 <- c(qnorm(0.025, mean_mu, sqrt(mean_var)),
         qnorm(0.975, mean_mu, sqrt(mean_var)))

## Figure 2.3: compare 2 approaches of LABPMM
plot(mu, dnorm(mu, mean_mu, sqrt(mean_var)), col = "blue",
     type = "l", 
     xlab = expression(mu),
     ylab = expression(paste("p(", mu, " | D)")))
for (i in 1:10){
  lines(mu, dnorm(mu, mu_grid[i], sqrt(var_grid[i])), col = "red", type = "l", lty = 2)
}

legend("topright", legend=c("Nested Laplace appr.", 
                            expression(paste("p(", mu, " |", delta^{(m)},", D)"))),
       col=c("blue", "red"), lty=1:2, cex=0.8)


plot(mu, dnorm(mu, mean_mu, sqrt(mean_var)), col = "blue",
     type = "l", 
     xlab = expression(mu),
     ylab = expression(paste("p(", mu, " | D)")))
lines(mu, dnorm(mu, mu_est, sqrt(var_est)), col = "red", type = "l", lty = 2)
legend("topright", legend=c("Approach 1", 
                            "Approach 2"),
       col=c("blue", "red"), lty=1:2, cex=0.8)
abline(v = c(mean_mu, mu_est),col=c("blue", "red"), lty=1:2) 
ci2 <- c(qnorm(0.025, mu_est, sqrt(var_est)),
         qnorm(0.975, mu_est, sqrt(var_est)))

func <- function(x, a = 1, b=153, c=-50){
  a*x^2 + x*b + c 
}
del_est = uniroot(func, lower = 0, upper = 3, tol = 1e-9)$root

### INLA method: plot of posterior of nu and fixed effects

rm(list = ls())
fixEffect = c(0.3,1.5,0.5, 1)
pRecision = 1
numUnit = 50
obsPerUnit = 6
numFixEff = length(fixEffect)
set.seed(2)
randomIntercept = rnorm(numUnit, 0, sd=sqrt(1/pRecision))
source("labpmm_map.R")
source("labpmm_ci.R")
source("labpmm_grid.R")

generate_data <- function(fix.effect, random.intercept, num.unit, obs.per.unit){
  group = rep(1:num.unit,each=obs.per.unit)
  n = num.unit*obs.per.unit
  num.fix.effect = length(fix.effect)
  
  # Design matrix for the random intercept
  z1 = matrix(rep(1,obs.per.unit), ncol=1)
  z = kronecker(diag(num.unit), z1)
  
  # vector of parameter
  phi = c(fix.effect, random.intercept)
  
  # Design matrix has first column of 1
  x = matrix(nrow = n, ncol = num.fix.effect)
  x[,1] = rep(1, n) 
  x[,2] = runif(n, min = 0, max = 1) - 0.5
  x[,3] = rnorm(n, 0, 1)
  x[,4] = rbinom(n, 1, 0.5) - 0.5
  
  t = cbind(x, z)
  
  lambda = exp(t%*%phi)
  
  #outcome
  
  y = rpois(n, lambda = lambda) 
  # if y is too large then phi model diverge => should remove outliers of y.
  # In this case, we set constraint for y
  # for (i in 1:n){
  #   y[i] = ifelse(y[i] < 5000, y[i], 5000)
  # }
  # 
  data = cbind(group, y, x[,-c(1)])
  return(data)
}

gammaA=0.01
gammaB=0.01
fixEffVar = 10^4

set.seed(1)
data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, num.unit = numUnit, obs.per.unit = obsPerUnit)
t = data[,-c(1:2)]
y = data[, 2]

z1 = matrix(rep(1,obsPerUnit), ncol=1)
# z: design matrix of random effect intercept b0
z = kronecker(diag(numUnit), z1)
ones = rep(1, length(y))
t = cbind(ones, t,z)
t = as.matrix(t)

### INLA
data2 <- data.frame(data)
colnames(data2)[3:5] = c("RV1", "RV2", "RV3")

prior.fixed <- list(mean.intercept = 0, prec.intercept = 0.0001,
                    mean = 0, prec = 0.0001)
prec.prior <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01)))

model.inla <- inla(y ~ RV1 + RV2 + RV3
                   + f(group, model ="iid", hyper = prec.prior),
                   data = data2, family = "poisson",
                   control.fixed = prior.fixed)

### MCMC

cat("model
  {
    for (j in 1:J){
      b[j] ~ dnorm(0, tau)
    }
    for (i in 1:N){
      rate[i] ~ dpois(phi[i])
      log(phi[i]) <- beta0 + beta1*RV1[i] + beta2*RV2[i] + beta3*RV3[i]
                  + b[group[i]]
    }
    tau ~ dgamma(0.01,0.01)
    variance = 1/tau
    beta0 ~ dnorm(0,1.0E-4)
    beta1 ~ dnorm(0,1.0E-4)
    beta2 ~ dnorm(0,1.0E-4)
    beta3 ~ dnorm(0,1.0E-4)
  }", file="simulation2.model.txt")

t1 = data.frame(data[,-c(1:2)])
colnames(t1)[1:3] = c("RV1", "RV2", "RV3")
y = data[, 2]
group = data[,1]
N <- length(y)
J <- n_distinct(group)

my.data <- list(rate=y, RV1 = t1$RV1, RV2 = t1$RV2, 
                RV3 = t1$RV3, 
                N = N, J = J, group = group)

my.inits <- list(
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 111, beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, tau=1),
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 112, beta0=0.1,beta1=0.1, beta2=0.1,beta3=0.1, tau=1),
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 113, beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, tau=0.1)
)

parameters <- c("beta0","beta1","beta2","beta3", "tau","variance")

jags <- jags.model(file="simulation2.model.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 3)
update(jags,1500)
sim.sim <- coda.samples(model = jags,
                        variable.names = parameters,
                        n.iter=5000, 
                        thin=10)
sim.mcmc <- as.mcmc.list(sim.sim)
sim.combined <- combine.mcmc(sim.mcmc)

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
get_func_for_v <- function(v, gamma.a, gamma.b, fix.eff.var){
  ## Return function of posterior of eta 
  .cov.mat = diag(rep(c(fix.eff.var, 1/exp(v)), times = c(numFixEff,numUnit)))
  .inv.cov.mat <- chol2inv(chol(.cov.mat))
  .phi.mode <- get_phi_model(inv.cov.mat = .inv.cov.mat)
  .H <- get_hessian(phi = .phi.mode, inv.cov.mat = .inv.cov.mat)
  .la.var <- chol2inv(chol(.H))
  funcval <- get_minus_log_llh(c(.phi.mode), inv.cov.mat = .inv.cov.mat) + 
    gamma.b * exp(v) - (numUnit/2+gamma.a)*v -1/2* log(det(.la.var))
  return(funcval)
}
vMode <- optimize(get_func_for_v, c(-20, 20), tol = 0.0001, 
                  gamma.a = gammaA, gamma.b = gammaB, 
                  fix.eff.var = fixEffVar)
vMode <- vMode$minimum
1/exp(vMode)

i = 0
grid_n_plot = 50
v_grid_plot = seq(vMode - grid_n_plot/2*0.1, vMode + grid_n_plot/2*0.1, 0.1)
log_post_v = rep(0, length(v_grid_plot))
for (v in v_grid_plot){
  i = i + 1
  log_post_v[i] = -get_func_for_v(v, gamma.a = gammaA, 
                                  gamma.b = gammaB, 
                                  fix.eff.var = fixEffVar)
}
log_post_v_scale = log_post_v - max(log_post_v)
post_v_scale = exp(log_post_v_scale)
post_v_norm_plot = post_v_scale/sum(post_v_scale)

i = 0
grid_n = 10
v_grid = seq(vMode - grid_n/2*0.15, vMode + grid_n/2*0.15, 0.15)
log_post_v = rep(0, length(v_grid))
for (v in v_grid){
  i = i + 1
  log_post_v[i] = -get_func_for_v(v, gamma.a = gammaA, 
                                  gamma.b = gammaB, 
                                  fix.eff.var = fixEffVar)
}
log_post_v_scale = log_post_v - max(log_post_v)
post_v_scale = exp(log_post_v_scale)
post_v_norm = post_v_scale/sum(post_v_scale)
mean_v = sum(v_grid*post_v_norm)
mode_v = v_grid[which.max(post_v_norm)]
mean_v
mode_v
m2_v = sum(v_grid^2*post_v_norm)
var_v = m2_v - mean_v^2
sd_v = sqrt(var_v)
lb_95_v = mean_v - 1.96*sd_v
ub_95_v = mean_v + 1.96*sd_v
ub_95_var = 1/exp(lb_95_v)
lb_95_var = 1/exp(ub_95_v)
ci_95_var = c(lb_95_var, ub_95_var)

# plot for posterior of \nu (Figure 3.1)
par(mar = c(5,5,3,3))
plot(v_grid_plot, post_v_norm_plot, col = "blue", 
     type = "l", xlim = c(-3, 1.5), xlab = expression(nu),
     ylab = expression(paste("p(", nu, " | D)")))
points(x = v_grid, y = rep(-0.0003, length(v_grid)), 
       col = "red", pch = 17, cex = 0.5)
abline(v=c(mean_v, mode_v, lb_95_v, ub_95_v), col = c("blue", "red", "black", "black"))

# plot for posterior of variance (Figure 3.6)
par(mar = c(5,5,3,3))
plot(1/exp(v_grid_plot), post_v_norm_plot, col = "blue", 
     type = "l", xlim=c(0,5), xlab = expression(sigma^2),
     ylab = expression(paste("p(", sigma^2, " | D)")))
mean_var = sum(1/exp(v_grid)*post_v_norm)
abline(v=c(mean_var, 1/exp(mode_v), ci_95_var), col = c("blue", "red", "black", "black"))


# plot compare multiple methods for variance estimate 

fpmm <- glmer(y ~ RV1 + RV2 + RV3 +(1|group), data=data2, family="poisson")


inv <- function(x) 1/x
prec <- model.inla$marginals.hyperpar$`Precision for group`
marg.var <- inla.tmarginal(inv, prec)
xaxis = seq(-2, 3, by=0.01)
plot(sim.combined[,6], trace = FALSE, density=TRUE,
     type = "l", 
     xlab = expression(sigma^2),
     ylab = expression(paste(hat(p), "(", sigma^2, " | D)")),
     main = "")
lines(marg.var, col = 'red', type = 'l', lty = 2)
set.seed(12)
y2 <- 1/exp(rnorm(xaxis, mean_v, sd_v))
lines(density(y2), 
      col = "blue",
      type = "l")
legend("topright", legend=c("LABPMM", 
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "green"), lty=c(1,2,1), cex=0.8)
abline(v=c(4,as.data.frame(VarCorr(fpmm))$vcov), lty = 1:2)


phi.mode = function(v){
  cov.mat = diag(rep(c(fixEffVar, 1/exp(v)), times = c(numFixEff,numUnit)))
  inv.cov.mat <- chol2inv(chol(cov.mat))
  phi.mode <- get_phi_model(inv.cov.mat = inv.cov.mat)
  return(phi.mode)
}

var_la = function(v){
  cov.mat = diag(rep(c(fixEffVar, 1/exp(v)), times = c(numFixEff,numUnit)))
  inv.cov.mat <- chol2inv(chol(cov.mat))
  phi.mode <- get_phi_model(inv.cov.mat = inv.cov.mat)
  hessianMode <- get_hessian(phi = phi.mode, inv.cov.mat = inv.cov.mat)
  var_la <- chol2inv(chol(hessianMode))
  return (var_la)
}

# plots for fixed effects (Figure 4.1, 3.2, 3.4)
phi_grid = rep(0, length(v_grid))
var_grid = rep(0, length(v_grid))


j= 0
for (v in v_grid){
  j = j+1
  phi_grid[j] = phi.mode(v)[1]
  var_grid[j] <- var_la(v)[1,1]
}
mean_phi <- sum(phi_grid*post_v_norm)
mean_var <- sum(var_grid*post_v_norm)
ci1 <- c(qnorm(0.025, mean_phi, sqrt(mean_var)),
         qnorm(0.975, mean_phi, sqrt(mean_var)))

xaxis = seq(-0.5, 1.5, by=0.001)
plot(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
     type = "l", 
     ylim = c(0,3.5),
     xlab = expression(beta[0]),
     ylab = expression(paste("p(", beta[0], " | D)")))
for (i in 1:length(v_grid)){
  lines(xaxis, dnorm(xaxis, phi_grid[i], sqrt(var_grid[i])), col = "red", type = "l", lty = 2)
  lines(xaxis, dnorm(xaxis, phi_grid[i], sqrt(var_grid[i]))*post_v_norm[i], col = "green", type = "l")
}

legend("topright", legend=c("LABPMM-grid", 
                            expression(paste("unweighted p(", beta[0], " |", nu^{(m)},", D)")),
                            expression(paste("weighted p(", beta[0], " |", nu^{(m)},", D)"))),
       col=c("blue", "red", "green"), lty=1:2, cex=0.8)

phi_est = phi.mode(vMode)[1]
var_est <- var_la(vMode)[1,1]
intercept.inla <- model.inla$marginals.fixed[[1]]

# library(MCMCvis)
# MCMCtrace(sim.sim,
#           params = 'beta0', 
#           type = 'density', 
#           ind = FALSE, 
#           pdf = FALSE,
#           ylim = c(0,3),
#           xlim = c(-1,2),
#           xlab_den = expression(beta[0]),
#           ylab_den = expression(paste("p(", beta[0], " | D)")))
plot(sim.combined[,1], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0,2.5),
     xlab = expression(beta[0]),
     ylab = expression(paste("p(", beta[0], " | D)")),
     main = "")
lines(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
      type = "l", 
      xlab = expression(beta[0]),
      ylab = expression(paste("p(", beta[0], " | D)")))
lines(xaxis, dnorm(xaxis, phi_est, sqrt(var_est)), col = "red", type = "l", lty = 2)
lines(intercept.inla, col = 'black', type = 'l', lty = 2)
abline(v = fixef(fpmm)[1], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "black", "black"), lty=c(1,2,2,1), cex=0.8)
### compare only LABPMM
plot(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
      type = "l", 
      xlab = expression(beta[0]),
      ylab = expression(paste("p(", beta[0], " | D)")))
lines(xaxis, dnorm(xaxis, phi_est, sqrt(var_est)), col = "red", type = "l", lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP"),
       col=c("blue", "red", "black", "black"), lty=c(1,2,2,1), cex=0.8)
abline(v = c(mean_phi, phi_est),col=c("blue", "red"), lty=1:2) 


# plot for beta_1
phi_grid = rep(0, length(v_grid))
var_grid = rep(0, length(v_grid))

j= 0
for (v in v_grid){
  j = j+1
  phi_grid[j] = phi.mode(v)[2]
  var_grid[j] <- var_la(v)[2,2]
}
mean_phi <- sum(phi_grid*post_v_norm)
mean_var <- sum(var_grid*post_v_norm)
# ci1 <- c(qnorm(0.025, mean_phi, sqrt(mean_var)),
#          qnorm(0.975, mean_phi, sqrt(mean_var)))

xaxis = seq(0.8, 2.2, by=0.001)
plot(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
     type = "l", 
     xlab = expression(beta[1]),
     ylab = expression(paste("p(", beta[1], " | D)")))
for (i in 1:length(v_grid)){
  lines(xaxis, dnorm(xaxis, phi_grid[i], sqrt(var_grid[i])), col = "red", type = "l", lty = 2)
  lines(xaxis, dnorm(xaxis, phi_grid[i], sqrt(var_grid[i]))*post_v_norm[i], col = "green", type = "l")
}

legend("topright", legend=c("LABPMM-grid", 
                            expression(paste("unweighted p(", beta[1], " |", nu^{(m)},", D)")),
                            expression(paste("weighted p(", beta[1], " |", nu^{(m)},", D)"))),
       col=c("blue", "red", "green"), lty=1:2, cex=0.8)

phi_est = phi.mode(vMode)[2]
var_est <- var_la(vMode)[2,2]
intercept.inla <- model.inla$marginals.fixed[[2]]

plot(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
     type = "l", 
     xlab = expression(beta[1]),
     ylab = expression(paste("p(", beta[1], " | D)")))
lines(xaxis, dnorm(xaxis, phi_est, sqrt(var_est)), col = "red", type = "l", lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP"),
       col=c("blue", "red","black"), lty=c(1,2,2), cex=0.8)
abline(v = c(mean_phi, phi_est),col=c("blue", "red"), lty=1:2) 


plot(sim.combined[,2], trace = FALSE, density=TRUE,
     type = "l", 
     xlab = expression(beta[1]),
     ylab = expression(paste("p(", beta[1], " | D)")),
     main = "")
lines(model.inla$marginals.fixed[[2]], col = 'black', type = 'l', lty = 2)
lines(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
     type = "l")
lines(xaxis, dnorm(xaxis, phi_est, sqrt(var_est)), col = "red", type = "l", lty = 2)
abline(v = fixef(fpmm)[2], lty=2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red","black", "black"), lty=c(1,2,2,1), cex=0.8)

# plot for beta_2
phi_grid = rep(0, length(v_grid))
var_grid = rep(0, length(v_grid))

j= 0
for (v in v_grid){
  j = j+1
  phi_grid[j] = phi.mode(v)[3]
  var_grid[j] <- var_la(v)[3,3]
}
mean_phi <- sum(phi_grid*post_v_norm)
mean_var <- sum(var_grid*post_v_norm)

xaxis = seq(0.4, 0.7, by=0.001)
plot(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
     type = "l", 
     ylim = c(0,11),
     xlab = expression(beta[2]),
     ylab = expression(paste("p(", beta[2], " | D)")))
for (i in 1:length(v_grid)){
  lines(xaxis, dnorm(xaxis, phi_grid[i], sqrt(var_grid[i])), col = "red", type = "l", lty = 2)
  lines(xaxis, dnorm(xaxis, phi_grid[i], sqrt(var_grid[i]))*post_v_norm[i], col = "green", type = "l")
}

legend("topright", legend=c("LABPMM-grid", 
                            expression(paste("unweighted p(", beta[2], " |", nu^{(m)},", D)")),
                            expression(paste("weighted p(", beta[2], " |", nu^{(m)},", D)"))),
       col=c("blue", "red", "green"), lty=1:2, cex=0.8)

phi_est = phi.mode(vMode)[3]
var_est <- var_la(vMode)[3,3]
intercept.inla <- model.inla$marginals.fixed[[3]]

plot(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
     type = "l", 
     xlab = expression(beta[2]),
     ylab = expression(paste("p(", beta[2], " | D)")))
lines(xaxis, dnorm(xaxis, phi_est, sqrt(var_est)), col = "red", type = "l", lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP"),
       col=c("blue", "red","black"), lty=c(1,2,2), cex=0.8)
abline(v = c(mean_phi, phi_est),col=c("blue", "red"), lty=1:2) 


plot(sim.combined[,3], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0,11),
     xlim = c(0.4, 0.7),
     xlab = expression(beta[2]),
     ylab = expression(paste("p(", beta[2], " | D)")),
     main = "")
lines(intercept.inla, col = 'black', type = 'l', lty = 2)
lines(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
      type = "l")
lines(xaxis, dnorm(xaxis, phi_est, sqrt(var_est)), col = "red", type = "l", lty = 2)
abline(v = fixef(fpmm)[3], lty=2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red","black", "black"), lty=c(1,2,2,1), cex=0.8)

# plot for beta_3
phi_grid = rep(0, length(v_grid))
var_grid = rep(0, length(v_grid))

j= 0
for (v in v_grid){
  j = j+1
  phi_grid[j] = phi.mode(v)[4]
  var_grid[j] <- var_la(v)[4,4]
}
mean_phi <- sum(phi_grid*post_v_norm)
mean_var <- sum(var_grid*post_v_norm)

xaxis = seq(0.6, 1.2, by=0.001)
plot(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
     type = "l", 
     ylim = c(0,6),
     xlab = expression(beta[3]),
     ylab = expression(paste("p(", beta[3], " | D)")))
for (i in 1:length(v_grid)){
  lines(xaxis, dnorm(xaxis, phi_grid[i], sqrt(var_grid[i])), col = "red", type = "l", lty = 2)
  lines(xaxis, dnorm(xaxis, phi_grid[i], sqrt(var_grid[i]))*post_v_norm[i], col = "green", type = "l")
}

legend("topright", legend=c("LABPMM-grid", 
                            expression(paste("unweighted p(", beta[3], " |", nu^{(m)},", D)")),
                            expression(paste("weighted p(", beta[3], " |", nu^{(m)},", D)"))),
       col=c("blue", "red", "green"), lty=1:2, cex=0.8)

phi_est = phi.mode(vMode)[4]
var_est <- var_la(vMode)[4,4]
intercept.inla <- model.inla$marginals.fixed[[4]]

plot(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
     type = "l", 
     xlab = expression(beta[3]),
     ylab = expression(paste("p(", beta[3], " | D)")))
lines(xaxis, dnorm(xaxis, phi_est, sqrt(var_est)), col = "red", type = "l", lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP"),
       col=c("blue", "red"), lty=c(1,2), cex=0.8)
abline(v = c(mean_phi, phi_est),col=c("blue", "red"), lty=1:2) 


plot(sim.combined[,4], trace = FALSE, density=TRUE,
     type = "l", 
     xlab = expression(beta[3]),
     ylab = expression(paste("p(", beta[3], " | D)")),
     main = "")
lines(intercept.inla, col = 'black', type = 'l', lty = 2)
lines(xaxis, dnorm(xaxis, mean_phi, sqrt(mean_var)), col = "blue",
      type = "l")
lines(xaxis, dnorm(xaxis, phi_est, sqrt(var_est)), col = "red", type = "l", lty = 2)
abline(v = fixef(fpmm)[4], lty=2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red","black", "black"), lty=c(1,2,2,1), cex=0.8)



