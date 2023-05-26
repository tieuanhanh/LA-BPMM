# Simulation study 1

library(ggplot2)
library(tictoc)
library(dplyr)
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
library(microbenchmark)

# with N = 50 and m = 6 and variance can be changed 

# set.seed(121)
# write(sample(1:10000000,size=500),file="seeds2.txt",ncolumns=1)

### start from here

rm(list = ls())
seeds<-read.table("seeds2.txt",header=F)$V1

library(tictoc)
run = 500
fixEffect = c(0.3,1.5,0.5, 1)
pRecision = 1
true_var = 1/pRecision
numUnit = 50
obsPerUnit = 6
numFixEff = length(fixEffect)
set.seed(2)
randomIntercept = rnorm(numUnit, 0, sd=sqrt(1/pRecision))
hist(randomIntercept)
plot(density(randomIntercept))
var(randomIntercept)
mean(randomIntercept)

# for (i in 1:1000000){
#   set.seed(i)
#   randomIntercept = rnorm(numUnit, 0, sd=sqrt(1/pRecision))
#   if ((abs(mean(randomIntercept)) < 0.001) &
#       (abs(var(randomIntercept)-4)<0.001)){
#     print(i)
#     break
#   }
# }

source("generate_data.R")
source("labpmm_ci.R")
source("labpmm_map.R")
source("labpmm_grid.R")

### LAPMM

write(c("num", "intercept", "RV1", "RV2", "RV3",
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd", 
        "variance", 
        "ci.time"),
      "results_lapmm_simulation_study2c.txt", 
      ncolumns=100, 
      append=T)

for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, num.unit = numUnit, obs.per.unit = obsPerUnit)
  t = data[,-c(1:2)]
  y = data[, 2]
  tic = proc.time()
  model <- labpmm(y=y, t=t, numUnit=numUnit, obsPerUnit=obsPerUnit, 
                     numFixEff=numFixEff, 
                     gammaA = 0.01, gammaB = 0.01)
  ci.time <- proc.time() - tic
  ci.time <- ci.time[3]
  #write to file
  write(c(i, model$fix_effect, model$sd_fixEffect,
          model$variance, model$sd_var, ci.time),
        "results_lapmm_simulation_study2c.txt", ncolumns=100,
        append=T)
  print(i)
}

dat1<-read.table("results_lapmm_simulation_study2c.txt",header=T) #load results
dim(dat1)


beta = dat1[, c(2:5)]
sd = dat1[, c(6:9)]
median(dat1$ci.time)
median(dat1$run.time)
# [1] 0.027

beta_true = fixEffect
beta_mean = colMeans(beta)
beta_mean

### bias
bias = beta_mean - beta_true
bias

### MSE
res = sweep(beta, 2, beta_true)
mse = sqrt(colMeans(res^2))
mse


### ESE
eres = sweep(beta, 2, beta_mean)
ese = sqrt(colSums(eres^2)/(run-1))
ese

### Coverage of CI 90
lower_bound_90 = beta - sd*1.645
upper_bound_90 = beta + sd*1.645

cp_90 = rep(0,4)
for (j in 1:4){
  cp_90[j] <- mean(ifelse((lower_bound_90[,j] < fixEffect[j]) & 
                            (upper_bound_90[,j] > fixEffect[j]),1,0 ))
}
cp_90

### Coverage of CI 95
lower_bound_95 = beta - sd*1.96
upper_bound_95 = beta + sd*1.96

cp_95 = rep(0,4)
for (j in 1:4){
  cp_95[j] <- mean(ifelse((lower_bound_95[,j] < fixEffect[j]) & 
                            (upper_bound_95[,j] > fixEffect[j]),1,0 ))
}
cp_95

# use mode variance gives the better performance
variance = dat1$variance
# variance = dat1$var.mean
true_var = 1/pRecision
mean_var = mean(variance)
bias_var = mean(variance) - true_var
bias_var

mse_var = sqrt(mean((variance - true_var)^2))
mse_var

ese_var = sqrt(sum((variance - mean(variance))^2)/(run-1))
ese_var

cp_90_var = mean(ifelse((dat1$var.lb.90 < true_var) & 
                          (dat1$var.ub.90 > true_var),1,0 ))
cp_90_var

cp_95_var = mean(ifelse((dat1$var.lb.95 < true_var) & 
                          (dat1$var.ub.95 > true_var),1,0 ))
cp_95_var

datafr <- cbind(bias, ese, mse, cp_90, cp_95)
datafr = rbind(datafr, cbind(bias_var, ese_var, mse_var, cp_90_var, cp_95_var))
row.names(datafr) <- c("intercept", "beta1", "beta2", "beta3", "variance")
colnames(datafr) <- c("Bias", "ESE", "MSE", "CP90", "CP95")
round(datafr,3)

# Variance = 4
# Bias   ESE   MSE  CP90  CP95
# intercept  0.252 0.065 0.260 1.000 1.000
# beta1     -0.003 0.064 0.064 0.892 0.940
# beta2      0.000 0.016 0.016 0.936 0.972
# beta3      0.000 0.038 0.038 0.894 0.948
# variance   1.044 0.409 1.121 0.960 0.998

### LABPMM-ci
write(c("num", "intercept", "RV1", "RV2", "RV3",
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd",
        "variance", "variance.sd", 
        "var.lb.95", "var.ub.95", "var.lb.90", "var.ub.90",
        "model.time"),
      "results_lapmm_ci_simulation_study2.txt",
      ncolumns=100,
      append=T)

for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, num.unit = numUnit, obs.per.unit = obsPerUnit)
  t = data[,-c(1:2)]
  y = data[, 2]
  tic = proc.time()
  model <- labpmm_ci(y=y, t=t, numUnit=numUnit, obsPerUnit=obsPerUnit,
                       numFixEff=numFixEff,
                       gammaA = 0.01, gammaB = 0.01)
  run.time = proc.time() - tic
  run.time = run.time[3]
  #write to file
  write(c(i, model$fix_effect, model$sd_fixEffect,
          model$variance, model$sd_var,
          model$ci_95_var, model$ci_90_var, run.time),
        "results_lapmm_ci_simulation_study2.txt", ncolumns=100,
        append=T)
  print(i)
}

dat6<-read.table("results_lapmm_ci_simulation_study2.txt",header=T) #load results
dim(dat6)

###LABPMM grid
write(c("num", "intercept", "RV1", "RV2", "RV3",
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd",
        "variance", "variance.sd", 
        "var.lb.95", "var.ub.95", "var.lb.90", "var.ub.90",
        "model.time"),
      "results_lapmm_grid_simulation_study2c.txt",
      ncolumns=100,
      append=T)

for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, num.unit = numUnit, obs.per.unit = obsPerUnit)
  t = data[,-c(1:2)]
  y = data[, 2]
  tic = proc.time()
  model <- labpmm_grid(y=y, t=t, numUnit=numUnit, obsPerUnit=obsPerUnit,
                      numFixEff=numFixEff,
                      gammaA = 0.01, gammaB = 0.01)
  run.time = proc.time() - tic
  run.time = run.time[3]
  #write to file
  write(c(i, model$fix_effect, model$sd_fixEffect,
          model$variance, model$sd_var,
          model$ci_95_var, model$ci_90_var, run.time),
        "results_lapmm_grid_simulation_study2c.txt", ncolumns=100,
        append=T)
  print(i)
}

dat5<-read.table("results_lapmm_grid_simulation_study2c.txt",header=T) #load results
dim(dat5)


beta = dat5[, c(2:5)]
sd = dat5[, c(6:9)]
median(dat5$model.time)
# [1] 0.027

beta_true = fixEffect
beta_mean = colMeans(beta)
beta_mean

### bias
bias = beta_mean - beta_true
bias

### MSE
res = sweep(beta, 2, beta_true)
mse = sqrt(colMeans(res^2))
mse


### ESE
eres = sweep(beta, 2, beta_mean)
ese = sqrt(colSums(eres^2)/(run-1))
ese

### Coverage of CI 90
lower_bound_90 = beta - sd*1.645
upper_bound_90 = beta + sd*1.645

cp_90 = rep(0,4)
for (j in 1:4){
  cp_90[j] <- mean(ifelse((lower_bound_90[,j] < fixEffect[j]) &
                            (upper_bound_90[,j] > fixEffect[j]),1,0 ))
}
cp_90

### Coverage of CI 95
lower_bound_95 = beta - sd*1.96
upper_bound_95 = beta + sd*1.96

cp_95 = rep(0,4)
for (j in 1:4){
  cp_95[j] <- mean(ifelse((lower_bound_95[,j] < fixEffect[j]) &
                            (upper_bound_95[,j] > fixEffect[j]),1,0 ))
}
cp_95

# use mode variance gives the better performance
variance = dat5$variance
# variance = dat1$var.mean
true_var = 1/pRecision
mean_var = mean(variance)
bias_var = mean(variance) - true_var
bias_var

mse_var = sqrt(mean((variance - true_var)^2))
mse_var

ese_var = sqrt(sum((variance - mean(variance))^2)/(run-1))
ese_var

cp_90_var = mean(ifelse((dat5$var.lb.90 < true_var) &
                          (dat5$var.ub.90 > true_var),1,0 ))
cp_90_var

cp_95_var = mean(ifelse((dat5$var.lb.95 < true_var) &
                          (dat5$var.ub.95 > true_var),1,0 ))
cp_95_var

datafr <- cbind(bias, ese, mse, cp_90, cp_95)
datafr = rbind(datafr, cbind(bias_var, ese_var, mse_var, cp_90_var, cp_95_var))
row.names(datafr) <- c("intercept", "beta1", "beta2", "beta3", "variance")
colnames(datafr) <- c("Bias", "ESE", "MSE", "CP90", "CP95")
round(datafr,3)


### Likelihood MPP
library(Matrix)
library(lme4)
library(boot)

write(c("num", "intercept", "RV1", "RV2", "RV3","variance",
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd", 
        "var.lb.95", "var.ub.95", 
        "RV0.lb.95", "RV1.lb.95", "RV2.lb.95", "RV3.lb.95", 
        "RV0.ub.95", "RV1.ub.95", "RV2.ub.95", "RV3.ub.95", 
        "var.lb.90", "var.ub.90", 
        "RV0.lb.90", "RV1.lb.90", "RV2.lb.90", "RV3.lb.90", 
        "RV0.ub.90", "RV1.ub.90", "RV2.ub.90", "RV3.ub.90", 
        "model.time", "ci.time"),
      "results_llh_simulation_study2c.txt", 
      ncolumns=100, 
      append=T)

for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect, 
                       random.intercept = randomIntercept, 
                       num.unit = numUnit, 
                       obs.per.unit = obsPerUnit)
  data2 = data.frame(data)
  tic = proc.time()
  fpmm <- glmer(y ~ V3 + V4+ V5 +(1|group), data=data2, family="poisson")
  model.time <- proc.time() - tic
  model.time <- model.time[3]
  vcov.matrix = vcov(fpmm)
  variance = vcov.matrix[row(vcov.matrix) == col(vcov.matrix)]
  sd_fixeff = sqrt(variance)
  variance = as.data.frame(VarCorr(fpmm))$vcov
  llh_ci_95 <- confint.merMod(fpmm, method = "profile")
  llh_ci_90 <- confint.merMod(fpmm, method = "profile", level = 0.9)
  var.95.ci <- llh_ci_95[1,]
  lb.95.ci <- llh_ci_95[-c(1),1]
  ub.95.ci <- llh_ci_95[-c(1),2]
  var.90.ci <- llh_ci_90[1,]
  lb.90.ci <- llh_ci_90[-c(1),1]
  ub.90.ci <- llh_ci_90[-c(1),2]
  ci.time <- proc.time() - tic
  ci.time <- ci.time[3]
  #write to file
  write(c(i, fixef(fpmm), variance, sd_fixeff, 
          var.95.ci, lb.95.ci, ub.95.ci,
          var.90.ci, lb.90.ci, ub.90.ci,
          model.time, ci.time),"results_llh_simulation_study2c.txt", ncolumns=100,
        append=T)
  print(i)
}
dat2 <- read.table("results_llh_simulation_study2c.txt", header = T)
median(dat2$model.time)
median(dat2$ci.time)
#0.041
beta = dat2[, c(2:6)]
sd = dat2[, c(7:10)]

beta_true = c(fixEffect, 1/pRecision)
beta_mean = colMeans(beta)
beta_mean

### bias
bias = beta_mean - beta_true
bias

### MSE
res = sweep(beta, 2, beta_true)
mse = sqrt(colMeans(res^2))
mse

### ESE
eres = sweep(beta, 2, beta_mean)
ese = sqrt(colSums(eres^2)/(run-1))
ese

### Coverage of CI 90
# lower_bound_90 = beta[1:numFixEff] - sd*1.645
# upper_bound_90 = beta[1:numFixEff] + sd*1.645
# 
# cp_90 = rep(0,numFixEff)
# for (j in 1:numFixEff){
#   cp_90[j] <- mean(ifelse((lower_bound_90[,j] < fixEffect[j]) & 
#                             (upper_bound_90[,j] > fixEffect[j]),1,0 ))
# }
# cp_90

### Coverage of CI 95
# lower_bound_95 = beta[1:numFixEff] - sd*1.96
# upper_bound_95 = beta[1:numFixEff] + sd*1.96
# 
# cp_95 = rep(0,numFixEff)
# for (j in 1:numFixEff){
#   cp_95[j] <- mean(ifelse((lower_bound_95[,j] < fixEffect[j]) & 
#                             (upper_bound_95[,j] > fixEffect[j]),1,0 ))
# }
# cp_95


ci.lb.90 = cbind(dat2$RV0.lb.90, dat2$RV1.lb.90, 
                  dat2$RV2.lb.90, dat2$RV3.lb.90, dat2$var.lb.90^2)
ci.ub.90 = cbind(dat2$RV0.ub.90, dat2$RV1.ub.90, 
                 dat2$RV2.ub.90, dat2$RV3.ub.90, dat2$var.ub.90^2)
cp_90 = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  cp_90[j] <- mean(ifelse((ci.lb.90[,j] < beta_true[j]) & 
                                (ci.ub.90[,j] > beta_true[j]),1,0 ))
}
cp_90

ci.lb.95 = cbind(dat2$RV0.lb.95, dat2$RV1.lb.95, 
                 dat2$RV2.lb.95, dat2$RV3.lb.95, (dat2$var.lb.95)^2)
ci.ub.95 = cbind(dat2$RV0.ub.95, dat2$RV1.ub.95, 
                 dat2$RV2.ub.95, dat2$RV3.ub.95, (dat2$var.ub.95)^2)
cp_95 = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  cp_95[j] <- mean(ifelse((ci.lb.95[,j] < beta_true[j]) & 
                            (ci.ub.95[,j] > beta_true[j]),1,0 ))
}
cp_95


datafr <- cbind(bias, ese, mse, cp_90, cp_95)
row.names(datafr) <- c("intercept", "beta1", "beta2", "beta3", "variance")
colnames(datafr) <- c("Bias", "ESE", "MSE", "CP90", "CP95")
round(datafr,3)
# Bias   ESE   MSE  CP90  CP95
# intercept  0.066 0.059 0.088 1.000 1.000
# beta1     -0.003 0.125 0.125 0.896 0.940
# beta2      0.000 0.036 0.036 0.908 0.956
# beta3      0.003 0.078 0.077 0.900 0.948
# variance   0.246 0.125 0.276 0.964 0.992


###INLA
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(INLA)

inv <- function(x) 1/x
write(c("num", "intercept", "RV1", "RV2", "RV3", 
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd",
        "int.lb.95", "int.lb.90", "int.med", "int.ub.90", "int.ub.95",
        "rv1.lb.95", "rv1.lb.90", "rv1.med", "rv1.ub.90", "rv1.ub.95",
        "rv2.lb.95", "rv2.lb.90", "rv2.med", "rv2.ub.90", "rv2.ub.95",
        "rv3.lb.95", "rv3.lb.90", "rv3.med", "rv3.ub.90", "rv3.ub.95",
        "var.mean",
        "var.lb.95", "var.lb.90", "var.med", "var.ub.90", "var.ub.95",
        "var.mode",
        "var.sd",
        "model.time",
        "tot.time"),"results_inla_simulation_study2c.txt", 
      ncolumns=100, 
      append=T)

for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect, 
                       random.intercept = randomIntercept, 
                       num.unit = numUnit, 
                       obs.per.unit = obsPerUnit)
  data2 = data.frame(data)
  colnames(data2)[1:5] = c("group", "y", "RV1", "RV2", "RV3")
  tic <- proc.time()
  prior.fixed <- list(mean.intercept = 0, prec.intercept = 0.0001,
                      mean = 0, prec = 0.0001)
  prec.prior <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01),
                                 initial = 4, fixed = FALSE))
  model.inla <- inla(y ~ RV1 + RV2 + RV3
                     + f(group, model ="iid", hyper = prec.prior),
                     data = data2, family = "poisson",
                     control.fixed = prior.fixed)
  # model.inla <- inla(y ~ RV1 + RV2 + RV3
  #                    + f(group, model ="iid"),
  #                    data = data2, family = "poisson")
  model.time <- proc.time() - tic
  model.time <- model.time[3]
  
  ### get the credible interval of fixed effects:
  marg.int <- model.inla$marginals.fixed$`(Intercept)`
  int.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975),marg.int)
  marg.rv1 <- model.inla$marginals.fixed$RV1
  marg.rv2 <- model.inla$marginals.fixed$RV2
  marg.rv3 <- model.inla$marginals.fixed$RV3
  rv1.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975),marg.rv1)
  rv2.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975),marg.rv2)
  rv3.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975),marg.rv3)
  inv <- function(x) 1/x
  prec <- model.inla$marginals.hyperpar$`Precision for group`
  marg.var <- inla.tmarginal(inv, prec)
  var.mean <- inla.emarginal(function(x) x, marg.var)
  var.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975), marg.var)
  var.mode <- inla.mmarginal(marg.var)
  mm <- inla.emarginal(function(x) x^2, marg.var)
  var.sd <- sqrt(mm - var.mean^2)
  tot.time <- proc.time() - tic
  tot.time <- tot.time[3]
  
  #write to file
  write(c(i, model.inla$summary.fixed[,1], 
          model.inla$summary.fixed[,2], 
          int.quant,
          rv1.quant,
          rv2.quant,
          rv3.quant,
          var.mean,
          var.quant,
          var.mode,
          var.sd,
          model.time,
          tot.time),"results_inla_simulation_study2c.txt", 
        ncolumns=100, 
        append=T)
  print(i)
}

dat4<-read.table("/Users/nana/rgit/LA-BPMM/Sim1_50.6.1/results_inla_simulation_study2_set_prior.txt",header=T) #load results
dim(dat4)

model_time = quantile(dat4$model.time, 0.5)
model_time

tot_time = quantile(dat4$tot.time, 0.5)
tot_time

beta = cbind(dat4[, c(2:5)], dat4$var.mean)
beta_true = c(fixEffect, 1/pRecision)
beta_mean = colMeans(beta)
beta_mean


### bias
bias = beta_mean - beta_true

bias


### MSE
res = sweep(beta, 2, beta_true)
mse = sqrt(colMeans(res^2))
mse


### ESE
eres = sweep(beta, 2, beta_mean)
ese = sqrt(colSums(eres^2)/(run-1))
ese

ci_90 = cbind(mean(ifelse((dat4$int.lb.90 < fixEffect[1]) 
                          & (dat4$int.ub.90 > fixEffect[1]),1,0 )),
              mean(ifelse((dat4$rv1.lb.90 < fixEffect[2]) 
                          & (dat4$rv1.ub.90 > fixEffect[2]),1,0 )),
              mean(ifelse((dat4$rv2.lb.90 < fixEffect[3]) 
                          & (dat4$rv2.ub.90 > fixEffect[3]),1,0 )),
              mean(ifelse((dat4$rv3.lb.90 < fixEffect[4]) 
                          & (dat4$rv3.ub.90 > fixEffect[4]),1,0 )),
              mean(ifelse((dat4$var.lb.90 < true_var) & 
                            (dat4$var.ub.90 > true_var),1,0 ))
)
ci_90


ci_95 = cbind(mean(ifelse((dat4$int.lb.95 < fixEffect[1]) 
                          & (dat4$int.ub.95 > fixEffect[1]),1,0 )),
              mean(ifelse((dat4$rv1.lb.95 < fixEffect[2]) 
                          & (dat4$rv1.ub.95 > fixEffect[2]),1,0 )),
              mean(ifelse((dat4$rv2.lb.95 < fixEffect[3]) 
                          & (dat4$rv2.ub.95 > fixEffect[3]),1,0 )),
              mean(ifelse((dat4$rv3.lb.95 < fixEffect[4]) 
                          & (dat4$rv3.ub.95 > fixEffect[4]),1,0 )),
              mean(ifelse((dat4$var.lb.95 < true_var) & 
                            (dat4$var.ub.95 > true_var),1,0 ))
)
ci_95


datafr <- cbind(bias, ese, mse, as.numeric(ci_90)*100, as.numeric(ci_95)*100)
row.names(datafr) <- c("intercept", "beta1", "beta2", "beta3", "variance")
colnames(datafr) <- c("Bias", "ESE", "MSE", "Cov 90", "Cov 95")

datafr2 = data.frame(datafr)
round(datafr2, 3)

### Rjags
library(coda)
library(rjags)
library(readr)
library(dplyr)
library(runjags)

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

write(c("num", "intercept", "RV1", "RV2", "RV3", "prec.mean", "var.mean",
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd", "prec.sd", "var.sd",

        "int.lb.95.hpd", "rv1.lb.95.hpd", "rv2.lb.95.hpd",
        "rv3.lb.95.hpd", "pre.lb.95.hpd", "var.lb.95.hpd",
        "int.ub.95.hpd", "rv1.ub.95.hpd", "rv2.ub.95.hpd",
        "rv3.ub.95.hpd", "pre.ub.95.hpd", "var.ub.95.hpd",


        "int.lb.90.hpd", "rv1.lb.90.hpd", "rv2.lb.90.hpd",
        "rv3.lb.90.hpd", "pre.lb.90.hpd", "var.lb.90.hpd",
        "int.ub.90.hpd", "rv1.ub.90.hpd", "rv2.ub.90.hpd",
        "rv3.ub.90.hpd", "pre.ub.90.hpd", "var.ub.90.hpd",

        "int.lb.95", "int.lb.90", "int.med", "int.ub.90", "int.ub.95",
        "rv1.lb.95", "rv1.lb.90", "rv1.med", "rv1.ub.90", "rv1.ub.95",
        "rv2.lb.95", "rv2.lb.90", "rv2.med", "rv2.ub.90", "rv2.ub.95",
        "rv3.lb.95", "rv3.lb.90", "rv3.med", "rv3.ub.90", "rv3.ub.95",
        "pre.lb.95", "pre.lb.90", "pre.med", "pre.ub.90", "pre.ub.95",
        "var.lb.95", "var.lb.90", "var.med", "var.ub.90", "var.ub.95",

        "model.time",
        "tot.time"),"results_rjags_simulation_study2b.txt",
      ncolumns=69,
      append=T)
for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, num.unit = numUnit, obs.per.unit = obsPerUnit)
  t = data.frame(data[,-c(1:2)])
  colnames(t)[1:3] = c("RV1", "RV2", "RV3")
  y = data[, 2]
  group = data[,1]
  N <- length(y)
  J <- n_distinct(group)

  my.data <- list(rate=y, RV1 = t$RV1, RV2 = t$RV2,
                  RV3 = t$RV3,
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

  # tic(i)
  tic <- proc.time()
  jags <- jags.model(file="simulation2.model.txt",
                     data = my.data,
                     inits = my.inits,
                     n.chains = 3)
  update(jags,1500)
  sim.sim <- coda.samples(model = jags,
                          variable.names = parameters,
                          n.iter=5000,
                          thin=10)
  model.time <- proc.time() - tic
  model.time <- model.time[3]

  model_obj = summary(sim.sim)

  sim.mcmc <- as.mcmc.list(sim.sim)
  sim.combined <- combine.mcmc(sim.mcmc)
  hpd.90 = HPDinterval(sim.combined, 0.90)
  hpd.95 = HPDinterval(sim.combined, 0.95)
  quant.beta0 = quantile(sim.combined[,1], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.beta1 = quantile(sim.combined[,2], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.beta2 = quantile(sim.combined[,3], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.beta3 = quantile(sim.combined[,4], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.prec = quantile(sim.combined[,5], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.var = quantile(sim.combined[,6], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant = cbind(quant.beta0,
                quant.beta1,
                quant.beta2,
                quant.beta3,
                quant.prec,
                quant.var)

  tot.time <- proc.time() - tic
  tot.time <- tot.time[3]

  #write to file
  write(c(i, model_obj$statistics[,1], model_obj$statistics[,2],
          hpd.95, hpd.90, quant, model.time, tot.time),
        "results_rjags_simulation_study2b.txt", ncolumns=69,
        append=T)
  print(i)
}

########

dat3<-read.table("results_rjags_simulation_study2b.txt",header=T) #load results
dim(dat3)
attach(dat3)

beta = dat3[, c(2:5,7)]

quantile(dat3$model.time, c(0.5))
quantile(dat3$tot.time, c(0.5))


beta_true = c(fixEffect, 1/pRecision)
beta_mean = colMeans(beta)

### bias
bias = beta_mean - beta_true
bias


### MSE
res = sweep(beta, 2, beta_true)
mse = sqrt(colMeans(res^2))
mse


### ESE
eres = sweep(beta, 2, beta_mean)
ese = sqrt(colSums(eres^2)/(run-1))
ese


### Coverage of CI 90 HPD
lower_bound_90.hpd = cbind(dat3$int.lb.90.hpd,
                           dat3$rv1.lb.90.hpd,
                           dat3$rv2.lb.90.hpd,
                           dat3$rv3.lb.90.hpd,
                           dat3$var.lb.90.hpd)

upper_bound_90.hpd = cbind(dat3$int.ub.90.hpd,
                           dat3$rv1.ub.90.hpd,
                           dat3$rv2.ub.90.hpd,
                           dat3$rv3.ub.90.hpd,
                           dat3$var.ub.90.hpd)

ci_90.hpd = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  ci_90.hpd[j] <- mean(ifelse((lower_bound_90.hpd[,j] < beta_true[j]) &
                                (upper_bound_90.hpd[,j] > beta_true[j]),1,0 ))
}
# ci_90 = colMeans(ci_table_90)
ci_90.hpd
# 1.000 0.896 0.912 0.902

### Coverage of CI 90 ETI
lower_bound_90.et = cbind(dat3$int.lb.90,
                          dat3$rv1.lb.90,
                          dat3$rv2.lb.90,
                          dat3$rv3.lb.90,
                          dat3$var.lb.90)

upper_bound_90.et = cbind(dat3$int.ub.90,
                          dat3$rv1.ub.90,
                          dat3$rv2.ub.90,
                          dat3$rv3.ub.90,
                          dat3$var.ub.90)

ci_90.et = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  ci_90.et[j] <- mean(ifelse((lower_bound_90.et[,j] < beta_true[j]) &
                               (upper_bound_90.et[,j] > beta_true[j]),1,0 ))
}
# ci_90 = colMeans(ci_table_90)
ci_90.et

### Coverage of CI 95 HPD
lower_bound_95.hpd = cbind(dat3$int.lb.95.hpd,
                           dat3$rv1.lb.95.hpd,
                           dat3$rv2.lb.95.hpd,
                           dat3$rv3.lb.95.hpd,
                           dat3$var.lb.95.hpd)

upper_bound_95.hpd = cbind(dat3$int.ub.95.hpd,
                           dat3$rv1.ub.95.hpd,
                           dat3$rv2.ub.95.hpd,
                           dat3$rv3.ub.95.hpd,
                           dat3$var.ub.95.hpd)

ci_95.hpd = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  ci_95.hpd[j] <- mean(ifelse((lower_bound_95.hpd[,j] < beta_true[j]) &
                                (upper_bound_95.hpd[,j] > beta_true[j]),1,0 ))
}
ci_95.hpd

### Coverage of CI 95 ETI
lower_bound_95.et = cbind(dat3$int.lb.95,
                          dat3$rv1.lb.95,
                          dat3$rv2.lb.95,
                          dat3$rv3.lb.95,
                          dat3$var.lb.95)

upper_bound_95.et = cbind(dat3$int.ub.95,
                          dat3$rv1.ub.95,
                          dat3$rv2.ub.95,
                          dat3$rv3.ub.95,
                          dat3$var.ub.95)

ci_95.et = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  ci_95.et[j] <- mean(ifelse((lower_bound_95.et[,j] < beta_true[j]) &
                               (upper_bound_95.et[,j] > beta_true[j]),1,0 ))
}
ci_95.et


datafr <- cbind(bias, ese, mse, ci_90.et*100, ci_95.et*100,
                ci_90.hpd*100, ci_95.hpd*100)
row.names(datafr) <- c("intercept", "beta1", "beta2", "beta3", "variance")
colnames(datafr) <- c("Bias", "ESE", "MSE", "CP 90 ETI", "CP 95 ETI",
                      "CP 90 HPD", "CP 95 HPD")
round(datafr,3)
# Bias   ESE   MSE CP 90 ETI CP 95 ETI CP 90 HPD CP 95 HPD
# intercept 0.174 0.134 0.219      95.4      98.2      94.2      98.0
# beta1     0.014 0.336 0.336      91.6      96.4      91.2      95.6
# beta2     0.011 0.108 0.109      86.6      94.2      87.0      94.4
# beta3     0.022 0.197 0.198      92.4      95.6      92.0      95.4
# variance  0.218 0.320 0.387      99.4      99.8      99.6     100.0

