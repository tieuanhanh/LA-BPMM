# Simulation study for one random intercept and one random slope case
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

# set.seed(121)
# write(sample(1:10000000,size=500),file="seeds2.txt",ncolumns=1)

### start from here
rm(list = ls())
seeds<-read.table("seeds2.txt",header=F)$V1

run = 500
fixEffect = c(0.5,1,0.5, 1.5)
pRecision = c(1,2)
numUnit = 100
obsPerUnit = 6
numFixEff = length(fixEffect)
set.seed(2)
randomIntercept = rnorm(numUnit, 0, sd=sqrt(1/pRecision[1]))
randomSlope = rnorm(numUnit, 0, sd=sqrt(1/pRecision[2]))
source("labpmm2_map.R")
source("labpmm2_grid.R")
source("labpmm2_ci.R")
mean(randomIntercept)
var(randomIntercept)
mean(randomSlope)
var(randomSlope)

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

### LAPMM
### the model with gamma(0.01, 0.01) and log transform of eta)
# for (i in 1:1){
#   #generate data
#   set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
#   data = generate_data(fix.effect = fixEffect, 
#                        random.intercept = randomIntercept,
#                        random.slope = randomSlope,
#                        num.unit = numUnit, obs.per.unit = obsPerUnit)
#   t = data[,c(3:5)]
#   y = data[, 2]
#   r = data[, 6]
#   tic(i)
#   model <- lapmm2_log(y=y, t=t, r = r, numUnit=numUnit, obsPerUnit=obsPerUnit,
#                   gammaA = 0.01, gammaB = 0.01, 
#                   fixEffVar = 10^4)
#   toc(log = TRUE, quiet = TRUE)
#   #log.txt <- tic.log(format = TRUE)
#   log.lst <- tic.log(format = FALSE)
#   tic.clearlog()
#   timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
#   #write to file
#   write(c(i, model$fix_effect, model$variance, model$sd_fixEffect, timings),"results_lapmm_simulation_study3.txt", 
#         ncolumns=14, 
#         append=T)
# }
write(c("num", "intercept", "RV1", "RV2", "RV3", 
        "variance1", "variance2", 
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd", 
        "var1.lb.95", "var1.ub.95", "var1.lb.90", "var1.ub.90",
        "var2.lb.95", "var2.ub.95", "var2.lb.90", "var2.ub.90",
        "mean_var1", "mean_var2",
        "model.time", "tot.time"),
      "results_labpmm_simulation_study3.txt", 
      ncolumns=100, 
      append=T)

for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect, 
                       random.intercept = randomIntercept,
                       random.slope = randomSlope,
                       num.unit = numUnit, obs.per.unit = obsPerUnit)
  t = data[,c(3:5)]
  y = data[, 2]
  r = data[, 6]
  tic = proc.time()
  model <- labpmm2(y=y, t=t, r = r, numUnit=numUnit, obsPerUnit=obsPerUnit,
                      gammaA = 0.01, gammaB = 0.01, 
                      fixEffVar = 10^4)
  model.time = proc.time() - tic
  model.time = model.time[3]
  tic = proc.time()
  ci <- labpmm2_ci(y=y, t=t, r=r, numUnit, obsPerUnit=obsPerUnit, 
                   gammaA = 0.01, gammaB = 0.01, 
                   fixEffVar = 10^4)
  tot.time = proc.time() - tic
  tot.time = tot.time[3]
  #write to file
  write(c(i, model$fix_effect, model$variance, model$sd_fixEffect, 
          ci$ci_95_var1, ci$ci_90_var1, ci$ci_95_var2, ci$ci_90_var2,
        ci$mean_var, model.time, tot.time),
        "results_labpmm_simulation_study3.txt", 
        ncolumns=100, 
        append=T)
  print(i)
}

dat1<-read.table("results_labpmm_simulation_study3.txt",header=T) #load results
dim(dat1)

write(c("num", "intercept", "RV1", "RV2", "RV3",
        "variance1", "variance2", 
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd", 
        "var1.lb.95", "var1.ub.95", "var1.lb.90", "var1.ub.90",
        "var2.lb.95", "var2.ub.95", "var2.lb.90", "var2.ub.90",
         "model.time"),
      "results_labpmm_simulation_study3_grid.txt", 
      ncolumns=100, 
      append=T)

for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect, 
                       random.intercept = randomIntercept,
                       random.slope = randomSlope,
                       num.unit = numUnit, obs.per.unit = obsPerUnit)
  t = data[,c(3:5)]
  y = data[, 2]
  r = data[, 6]
  tic = proc.time()
  model <- labpmm2_grid(y=y, t=t, r = r, numUnit=numUnit, 
                        obsPerUnit=obsPerUnit,
                   gammaA = 0.01, gammaB = 0.01, 
                   fixEffVar = 10^4)
  model.time = proc.time() - tic
  model.time = model.time[3]
  #write to file
  write(c(i, model$fix_effect, model$mean_var,
          model$sd_fixEffect, 
          model$ci_95_var1, model$ci_90_var1, 
          model$ci_95_var2, model$ci_90_var2,
          model.time),
        "results_labpmm_simulation_study3_grid.txt", 
        ncolumns=100, 
        append=T)
  print(i)
}

dat1<-read.table("results_labpmm_simulation_study3_grid.txt",header=T) #load results
dim(dat1)

#for labpmm_grid
beta = dat1[, c(2:7)]
sd = dat1[, c(8:11)]
median(dat1$model.time)
# [1] 0.9545
median(dat1$tot.time)
# [1] 2.011
# grid: 2.87s

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
lower_bound_90 = beta[1:numFixEff] - sd*1.645
upper_bound_90 = beta[1:numFixEff] + sd*1.645

cp_90 = rep(0,4)
for (j in 1:4){
  cp_90[j] <- mean(ifelse((lower_bound_90[,j] < fixEffect[j]) & 
                            (upper_bound_90[,j] > fixEffect[j]),1,0 ))
}
cp_90

### Coverage of CI 95
lower_bound_95 = beta[1:numFixEff] - sd*1.96
upper_bound_95 = beta[1:numFixEff] + sd*1.96

cp_95 = rep(0,4)
for (j in 1:4){
  cp_95[j] <- mean(ifelse((lower_bound_95[,j] < fixEffect[j]) & 
                            (upper_bound_95[,j] > fixEffect[j]),1,0 ))
}
cp_95

# variance = dat1$var.mean
true_var = 1/pRecision

cp_90_var1 = mean(ifelse((dat1$var1.lb.90 < true_var[1]) & 
                          (dat1$var1.ub.90 > true_var[1]),1,0 ))
cp_90_var1

cp_95_var1 = mean(ifelse((dat1$var1.lb.95 < true_var[1]) & 
                          (dat1$var1.ub.95 > true_var[1]),1,0 ))
cp_95_var1

cp_90_var2 = mean(ifelse((dat1$var2.lb.90 < true_var[2]) & 
                           (dat1$var2.ub.90 > true_var[2]),1,0 ))
cp_90_var2

cp_95_var2 = mean(ifelse((dat1$var2.lb.95 < true_var[2]) & 
                           (dat1$var2.ub.95 > true_var[2]),1,0 ))
cp_95_var2

datafr <- cbind(bias, ese, mse, c(cp_90, cp_90_var1, cp_90_var2), 
                c(cp_95, cp_95_var1, cp_95_var2))
row.names(datafr) <- c("intercept", "beta1", "beta2", 
                       "beta3", "variance1", "variance2")
colnames(datafr) <- c("Bias", "ESE", "MSE", "CP90", "CP95")
round(datafr,3)

# Grid
# Bias   ESE   MSE  CP90  CP95
# intercept  0.125 0.053 0.136 0.996 1.000
# beta1     -0.008 0.134 0.134 0.868 0.936
# beta2     -0.001 0.038 0.038 0.882 0.942
# beta3     -0.003 0.081 0.081 0.906 0.946
# variance1  0.333 0.122 0.355 0.920 0.984
# variance2  0.180 0.222 0.285 0.914 0.954


### Likelihood MPP
library(Matrix)
library(lme4)
library(boot)

# for (i in 1:run){
#   #generate data
#   set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
#   data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, 
#                        random.slope = randomSlope, num.unit = numUnit, 
#                        obs.per.unit = obsPerUnit)
#   t = data[,-c(1:2)]
#   y = data[, 2]
#   data2 = data.frame(data)
#   tic(i)
#   fpmm <- glmer(y ~ V3 + V4+ V5 +(1+V6|group), data=data2, family="poisson")
#   vcov.matrix = vcov(fpmm)
#   var_fixeff = vcov.matrix[row(vcov.matrix) == col(vcov.matrix)]
#   sd_fixeff = sqrt(var_fixeff)
#   toc(log = TRUE, quiet = TRUE)
#   #log.txt <- tic.log(format = TRUE)
#   log.lst <- tic.log(format = FALSE)
#   tic.clearlog()
#   variance = as.data.frame(VarCorr(fpmm))$vcov
#   timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
#   #write to file
#   write(c(i, fixef(fpmm), variance, sd_fixeff, timings),"results_llh_simulation_study3.txt", ncolumns=13, 
#         append=T)
# }

# write(c("num", "intercept", "RV1", "RV2", "RV3", "var1", "var2",
#         "int.sd", "RV1.sd", "RV2.sd", "RV3.sd", "var1.sd", "var2.sd",
#         "var1.lb.90.basic", "var1.ub.90.basic", 
#         "var1.lb.90.bca", "var1.ub.90.bca",
#         "var1.lb.90.norm", "var1.ub.90.norm", 
#         "var1.lb.90.perc", "var1.ub.90.perc",
#         "var1.lb.95.basic", "var1.ub.95.basic", 
#         "var1.lb.95.bca", "var1.ub.95.bca",
#         "var1.lb.95.norm", "var1.ub.95.norm", 
#         "var1.lb.95.perc", "var1.ub.95.perc",
#         
#         "var2.lb.90.basic", "var2.ub.90.basic", 
#         "var2.lb.90.bca", "var2.ub.90.bca",
#         "var2.lb.90.norm", "var2.ub.90.norm", 
#         "var2.lb.90.perc", "var2.ub.90.perc",
#         "var2.lb.95.basic", "var2.ub.95.basic", 
#         "var2.lb.95.bca", "var2.ub.95.bca",
#         "var2.lb.95.norm", "var2.ub.95.norm", 
#         "var2.lb.95.perc", "var2.ub.95.perc",
#         
#         "model.time",
#         "tot.time"),"results_llh_simulation_study3_bootstrap.txt", 
#       ncolumns=100, 
#       append=T)

write(c("num", "intercept", "RV1", "RV2", "RV3", "var1", "var2",
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd", 
        "var1.lb.90", "var1.ub.90", 
        "var1.lb.95", "var1.ub.95", 
        "var2.lb.90", "var2.ub.90", 
        "var2.lb.95", "var2.ub.95", 
        "beta0.lb.90", "beta1.lb.90","beta2.lb.90","beta3.lb.90",
        "beta0.ub.90", "beta1.ub.90","beta2.ub.90","beta3.ub.90",
        "beta0.lb.95", "beta1.lb.95","beta2.lb.95","beta3.lb.95",
        "beta0.ub.95", "beta1.ub.95","beta2.ub.95","beta3.ub.95",
        "model.time",
        "tot.time"),"results_llh_simulation_study3.txt", 
      ncolumns=100, 
      append=T)

for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect, 
                       random.slope = randomSlope, 
                       random.intercept = randomIntercept, 
                       num.unit = numUnit, obs.per.unit = obsPerUnit)
  data2 = data.frame(data)
  
  tic <- proc.time()
  
  fpmm <- glmer(y ~ V3 + V4+ V5 +(1+V6|group), data=data2, family="poisson")
  vcov.matrix = vcov(fpmm)
  var_fixeff = vcov.matrix[row(vcov.matrix) == col(vcov.matrix)]
  sd_fixeff = sqrt(var_fixeff)
  variance = as.data.frame(VarCorr(fpmm))$vcov[1:2]
  
  model.time <- proc.time() - tic
  model.time <- model.time[3]
  ci.95 <- confint.merMod(fpmm, level = 0.95, method = "Wald")
  var1.ci.95 <- ci.95[1,]^2
  var2.ci.95 <- ci.95[3,]^2
  beta.lb.95 <- ci.95[c(4:7),1]
  beta.ub.95 <- ci.95[c(4:7),2]
  ci.90 <- confint.merMod(fpmm, level = 0.90, method = "Wald")
  var1.ci.90 <- ci.90[1,]^2
  var2.ci.90 <- ci.90[3,]^2
  beta.lb.90 <- ci.90[c(4:7),1]
  beta.ub.90 <- ci.90[c(4:7),2]
  
  ### calculate bootstrap CI
  # lpmm <- function(formula, data, indices) {
  #   d <- data[indices,] # allows boot to select sample
  #   fit <- glmer(formula, data=d, family="poisson")
  #   return(as.data.frame(VarCorr(fit))$vcov[1:2])
  # }
  # 
  # set.seed(1)
  # # bootstrapping with 2000 replications
  # 
  # results <- boot(data=data2, statistic=lpmm,
  #                 R=2000, formula=y ~ V3 + V4+ V5 +(1+V6|group))
  # var1.sd = sd(results$t[,1])
  # var2.sd = sd(results$t[,2])
  # 
  # var1.boot_ci_90 = boot.ci(results, conf = 0.9, 
  #                      type = c("norm", "basic", "perc", "bca"), 
  #                      index = 1)
  # x = var1.boot_ci_90$basic
  # var1.ci.basic.90 = x[(length(x)-1):length(x)]
  # x = var1.boot_ci_90$bca
  # var1.ci.bca.90 = x[(length(x)-1):length(x)]
  # x = var1.boot_ci_90$normal
  # var1.ci.norm.90 = x[(length(x)-1):length(x)]
  # x = var1.boot_ci_90$percent
  # var1.ci.perc.90 = x[(length(x)-1):length(x)]
  # 
  # var1.boot_ci_95 = boot.ci(results, 
  #                           type = c("norm", "basic", "perc", "bca"),
  #                           index = 1)
  # x = var1.boot_ci_95$basic
  # var1.ci.basic.95 = x[(length(x)-1):length(x)]
  # x = var1.boot_ci_95$bca
  # var1.ci.bca.95 = x[(length(x)-1):length(x)]
  # x = var1.boot_ci_95$normal
  # var1.ci.norm.95 = x[(length(x)-1):length(x)]
  # x = var1.boot_ci_95$percent
  # var1.ci.perc.95 = x[(length(x)-1):length(x)]
  # 
  # var1.ci <- cbind(var1.ci.basic.90, var1.ci.bca.90, var1.ci.norm.90,
  #                 var1.ci.perc.90, var1.ci.basic.95, var1.ci.bca.95,
  #                 var1.ci.norm.95, var1.ci.perc.95)
  # 
  # var2.boot_ci_90 = boot.ci(results, conf = 0.9, 
  #                           type = c("norm", "basic", "perc", "bca"), 
  #                           index = 2)
  # x = var2.boot_ci_90$basic
  # var2.ci.basic.90 = x[(length(x)-1):length(x)]
  # x = var2.boot_ci_90$bca
  # var2.ci.bca.90 = x[(length(x)-1):length(x)]
  # x = var2.boot_ci_90$normal
  # var2.ci.norm.90 = x[(length(x)-1):length(x)]
  # x = var2.boot_ci_90$percent
  # var2.ci.perc.90 = x[(length(x)-1):length(x)]
  # 
  # var2.boot_ci_95 = boot.ci(results,  
  #                           type = c("norm", "basic", "perc", "bca"), 
  #                           index = 2)
  # x = var2.boot_ci_95$basic
  # var2.ci.basic.95 = x[(length(x)-1):length(x)]
  # x = var2.boot_ci_95$bca
  # var2.ci.bca.95 = x[(length(x)-1):length(x)]
  # x = var2.boot_ci_95$normal
  # var2.ci.norm.95 = x[(length(x)-1):length(x)]
  # x = var2.boot_ci_95$percent
  # var2.ci.perc.95 = x[(length(x)-1):length(x)]
  # var2.ci <- cbind(var2.ci.basic.90, var2.ci.bca.90, var2.ci.norm.90,
  #                  var2.ci.perc.90, var2.ci.basic.95, var2.ci.bca.95,
  #                  var2.ci.norm.95, var2.ci.perc.95)
  # 
  tot.time <- proc.time() - tic
  tot.time <- tot.time[3]
  # 
  # write to file
  write(c(i, fixef(fpmm), variance, sd_fixeff, 
          var1.ci.90, var1.ci.95, var2.ci.90, var2.ci.95,
          beta.lb.90, beta.ub.90,
          beta.lb.95, beta.ub.95,
          model.time, tot.time),
        "results_llh_simulation_study3.txt",
        ncolumns=100,
        append=T)
  print(i)
}

dat2<-read.table("results_llh_simulation_study3.txt",header=T) #load results
dim(dat2)

beta = dat2[, c(2:7)]
sd = dat2[, c(8:11)]
median(dat2$model.time)
median(dat2$tot.time)

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

ci.lb.90 = cbind(dat2$beta0.lb.90, dat2$beta1.lb.90, 
                 dat2$beta2.lb.90, dat2$beta3.lb.90, 
                 dat2$var1.lb.90, dat2$var2.lb.90)
ci.ub.90 = cbind(dat2$beta0.ub.90, dat2$beta1.ub.90, 
                 dat2$beta2.ub.90, dat2$beta3.ub.90, 
                 dat2$var1.ub.90, dat2$var2.ub.90)
cp_90 = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  cp_90[j] <- mean(ifelse((ci.lb.90[,j] < beta_true[j]) & 
                            (ci.ub.90[,j] > beta_true[j]),1,0 ))
}
cp_90

ci.lb.95 = cbind(dat2$beta0.lb.95, dat2$beta1.lb.95, 
                 dat2$beta2.lb.95, dat2$beta3.lb.95, 
                 dat2$var1.lb.95, dat2$var2.lb.95)
ci.ub.95 = cbind(dat2$beta0.ub.95, dat2$beta1.ub.95, 
                 dat2$beta2.ub.95, dat2$beta3.ub.95, 
                 dat2$var1.ub.95, dat2$var2.ub.95)
cp_95 = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  cp_95[j] <- mean(ifelse((ci.lb.95[,j] < beta_true[j]) & 
                            (ci.ub.95[,j] > beta_true[j]),1,0 ))
}
cp_95

### Coverage of CI 90
# lower_bound_90 = beta[1:numFixEff] - sd*1.645
# upper_bound_90 = beta[1:numFixEff] + sd*1.645
# 
# ci_table_90 = matrix(0,run,numFixEff)
# for (j in 1:numFixEff){
#   for (i in 1:run){
#     ci_table_90[i,j] <- ifelse((lower_bound_90[i,j] < fixEffect[j]) & 
#                                  (upper_bound_90[i,j] > fixEffect[j]),1,0 )
#   }
# }
# cp_90 = colMeans(ci_table_90)
# cp_90
# 
# ### Coverage of CI 95
# lower_bound_95 = beta[1:numFixEff] - sd*1.96
# upper_bound_95 = beta[1:numFixEff] + sd*1.96
# 
# ci_table_95 = matrix(0,run,numFixEff)
# for (j in 1:numFixEff){
#   for (i in 1:run){
#     ci_table_95[i,j] <- ifelse((lower_bound_95[i,j] < fixEffect[j]) & 
#                                  (upper_bound_95[i,j] > fixEffect[j]),1,0 )
#   }
# }
# cp_95 = colMeans(ci_table_95)
# cp_95
# 
# attach(dat2)
# var1.lb.90 = cbind(dat2$var1.lb.90.basic, dat2$var1.lb.90.bca, 
#                dat2$var1.lb.90.norm, dat2$var1.lb.90.perc)
# var1.ub.90 = cbind(dat2$var1.ub.90.basic, dat2$var1.ub.90.bca, 
#                dat2$var1.ub.90.norm, dat2$var1.ub.90.perc)
# cp_90_var1 = rep(0,4)
# for (j in 1:4){
#   cp_90_var1[j] <- mean(ifelse((var1.lb.90[,j] < 1/pRecision[1]) & 
#                                  (var1.ub.90[,j] > 1/pRecision[1]),1,0 ))
# }
# cp_90_var1
# 
# var2.lb.90 = cbind(dat2$var2.lb.90.basic, dat2$var2.lb.90.bca, 
#                    dat2$var2.lb.90.norm, dat2$var2.lb.90.perc)
# var2.ub.90 = cbind(dat2$var2.ub.90.basic, dat2$var2.ub.90.bca, 
#                    dat2$var2.ub.90.norm, dat2$var2.ub.90.perc)
# cp_90_var2 = rep(0,4)
# for (j in 1:4){
#   cp_90_var2[j] <- mean(ifelse((var2.lb.90[,j] < 1/pRecision[2]) & 
#                                (var2.ub.90[,j] > 1/pRecision[2]),1,0 ))
# }
# cp_90_var2
# 
# var1.lb.95 = cbind(dat2$var1.lb.95.basic, dat2$var1.lb.95.bca, 
#                    dat2$var1.lb.95.norm, dat2$var1.lb.95.perc)
# var1.ub.95 = cbind(dat2$var1.ub.95.basic, dat2$var1.ub.95.bca, 
#                    dat2$var1.ub.95.norm, dat2$var1.ub.95.perc)
# cp_95_var1 = rep(0,4)
# for (j in 1:4){
#   cp_95_var1[j] <- mean(ifelse((var1.lb.95[,j] < 1/pRecision[1]) & 
#                                  (var1.ub.95[,j] > 1/pRecision[1]),1,0 ))
# }
# cp_95_var1
# 
# var2.lb.95 = cbind(dat2$var2.lb.95.basic, dat2$var2.lb.95.bca, 
#                    dat2$var2.lb.95.norm, dat2$var2.lb.95.perc)
# var2.ub.95 = cbind(dat2$var2.ub.95.basic, dat2$var2.ub.95.bca, 
#                    dat2$var2.ub.95.norm, dat2$var2.ub.95.perc)
# cp_95_var2 = rep(0,4)
# for (j in 1:4){
#   cp_95_var2[j] <- mean(ifelse((var2.lb.95[,j] < 1/pRecision[2]) & 
#                                  (var2.ub.95[,j] > 1/pRecision[2]),1,0 ))
# }
# cp_95_var2


datafr <- cbind(bias, ese, mse, cp_90, cp_95)
row.names(datafr) <- c("intercept", "beta1", "beta2", "beta3", "var1", "var2")
colnames(datafr) <- c("Bias", "ESE", "MSE", "CP90", "CP95")
round(datafr,3)
# bootstrap
# bias   ese   mse            
# intercept    0.068 0.057 0.089 0.998 1.000
# beta1       -0.005 0.134 0.134 0.866 0.930
# beta2        0.000 0.038 0.038 0.878 0.940
# beta3        0.002 0.081 0.081 0.900 0.948
# variance_ri  0.240 0.112 0.265 0.914 0.962
# variance_rs  0.131 0.218 0.254 0.906 0.970

# profile
# Bias   ESE   MSE  CP90  CP95
# intercept  0.068 0.057 0.089 0.998 1.000
# beta1     -0.005 0.134 0.134 0.870 0.934
# beta2      0.000 0.038 0.038 0.884 0.942
# beta3      0.002 0.081 0.081 0.904 0.948
# var1       0.240 0.112 0.265 0.954 0.990
# var2       0.131 0.218 0.254 0.930 0.958

### Rjags

library(coda)
library(rjags)
library(readr)
library(dplyr)
library(runjags)

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

# for (i in 1:run){
#   #generate data
#   set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
#   data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept,
#                        random.slope = randomSlope,
#                        num.unit = numUnit, obs.per.unit = obsPerUnit)
#   t = data.frame(data[,-c(1:2)])
#   colnames(t)[1:4] = c("RV1", "RV2", "RV3","RV4")
#   y = data[, 2]
#   group = data[,1]
#   N <- length(y)
#   J <- n_distinct(group)
#
#   my.data <- list(rate=y, RV1 = t$RV1, RV2 = t$RV2,
#                   RV3 = t$RV3, RV4 = t$RV4,
#                   N = N, J = J, group = group)
#
#   my.inits <- list(
#     list(beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, tau1=1, tau2=1),
#     list(beta0=0.1,beta1=0.1, beta2=0.1,beta3=0.1, tau1=0.5, tau2=1),
#     list(beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, tau1=0.5, tau2=0.5)
#   )
#
#   parameters <- c("beta0","beta1","beta2","beta3", "tau1", "tau2", "var1", "var2")
#
#   tic(i)
#   jags <- jags.model(file="simulation3.model.txt",
#                      data = my.data,
#                      inits = my.inits,
#                      n.chains = 3)
#   update(jags,1500)
#   sim.sim <- coda.samples(model = jags,
#                           variable.names = parameters,
#                           n.iter=5000,
#                           thin=10)
#   toc(log = TRUE, quiet = TRUE)
#   #log.txt <- tic.log(format = TRUE)
#   log.lst <- tic.log(format = FALSE)
#   tic.clearlog()
#   timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
#   model_obj = summary(sim.sim)
#   #write to file
#   write(c(i, model_obj$statistics[,1], model_obj$statistics[,2], timings),
#         "results_rjags_simulation_study3.txt", ncolumns=18,
#         append=T)
#   print(i)
# }

write(c("num", "intercept", "RV1", "RV2", "RV3",
        "prec1.mean", "prec2.mean", "var1.mean", "var2.mean",
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd",
        "prec1.sd", "prec2.sd", "var1.sd","var2.sd",

        "int.lb.95.hpd", "rv1.lb.95.hpd", "rv2.lb.95.hpd",
        "rv3.lb.95.hpd", "pre1.lb.95.hpd", "pre2.lb.95.hpd",
        "var1.lb.95.hpd", "var2.lb.95.hpd",
        "int.ub.95.hpd", "rv1.ub.95.hpd", "rv2.ub.95.hpd",
        "rv3.ub.95.hpd", "pre1.ub.95.hpd", "pre2.ub.95.hpd",
        "var1.ub.95.hpd", "var2.ub.95.hpd",


        "int.lb.90.hpd", "rv1.lb.90.hpd", "rv2.lb.90.hpd",
        "rv3.lb.90.hpd", "pre1.lb.90.hpd", "pre2.lb.90.hpd",
        "var1.lb.90.hpd", "var2.lb.90.hpd",
        "int.ub.90.hpd", "rv1.ub.90.hpd", "rv2.ub.90.hpd",
        "rv3.ub.90.hpd", "pre1.ub.90.hpd", "pre2.ub.90.hpd",
        "var1.ub.90.hpd", "var2.ub.90.hpd",

        "int.lb.95", "int.lb.90", "int.med", "int.ub.90", "int.ub.95",
        "rv1.lb.95", "rv1.lb.90", "rv1.med", "rv1.ub.90", "rv1.ub.95",
        "rv2.lb.95", "rv2.lb.90", "rv2.med", "rv2.ub.90", "rv2.ub.95",
        "rv3.lb.95", "rv3.lb.90", "rv3.med", "rv3.ub.90", "rv3.ub.95",
        "pre1.lb.95", "pre1.lb.90", "pre1.med", "pre1.ub.90", "pre1.ub.95",
        "pre2.lb.95", "pre2.lb.90", "pre2.med", "pre2.ub.90", "pre2.ub.95",
        "var1.lb.95", "var1.lb.90", "var1.med", "var1.ub.90", "var1.ub.95",
        "var2.lb.95", "var2.lb.90", "var2.med", "var2.ub.90", "var2.ub.95",

        "model.time",
        "tot.time"),"results_rjags_simulation_study3_test.txt",
      ncolumns=100,
      append=T)

for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect,
                       random.intercept = randomIntercept,
                       random.slope = randomSlope,
                       num.unit = numUnit, obs.per.unit = obsPerUnit)
  t = data.frame(data[,-c(1:2)])
  colnames(t)[1:4] = c("RV1", "RV2", "RV3", "RV4")
  y = data[, 2]
  group = data[,1]
  N <- length(y)
  J <- n_distinct(group)

  my.data <- list(rate=y, RV1 = t$RV1, RV2 = t$RV2,
                  RV3 = t$RV3, RV4 = t$RV4,
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

  tic <- proc.time()
  jags <- jags.model(file="simulation3.model.txt",
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
  # compute the statistics
  model_obj = summary(sim.sim)
  # combine 3 chains
  sim.mcmc <- as.mcmc.list(sim.sim)
  sim.combined <- combine.mcmc(sim.mcmc)
  hpd.90 = HPDinterval(sim.combined, 0.90)
  hpd.95 = HPDinterval(sim.combined, 0.95)
  quant.beta0 = quantile(sim.combined[,1], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.beta1 = quantile(sim.combined[,2], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.beta2 = quantile(sim.combined[,3], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.beta3 = quantile(sim.combined[,4], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.prec1 = quantile(sim.combined[,5], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.prec2 = quantile(sim.combined[,6], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.var1 = quantile(sim.combined[,7], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant.var2 = quantile(sim.combined[,8], c(0.025, 0.05, 0.5, 0.95, 0.975))
  quant = cbind(quant.beta0,
                quant.beta1,
                quant.beta2,
                quant.beta3,
                quant.prec1,
                quant.prec2,
                quant.var1,
                quant.var2)

  tot.time <- proc.time() - tic
  tot.time <- tot.time[3]

  #write to file
  write(c(i, model_obj$statistics[,1], model_obj$statistics[,2],
          hpd.95, hpd.90, quant, model.time, tot.time),
        "results_rjags_simulation_study3_test.txt", ncolumns=100,
        append=T)
  print(i)
}

dat3<-read.table("results_rjags_simulation_study3_test.txt",header=T) #load results
dim(dat3)
attach(dat3)

beta = dat3[, c(2:5,8:9)]

quantile(dat3$model.time, 0.5)
# 13.375

quantile(dat3$tot.time, 0.5)
# 13.403

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
lower_bound_90.hpd = cbind(int.lb.90.hpd,
                           rv1.lb.90.hpd,
                           rv2.lb.90.hpd,
                           rv3.lb.90.hpd,
                           var1.lb.90.hpd,
                           var2.lb.90.hpd)

upper_bound_90.hpd = cbind(int.ub.90.hpd,
                           rv1.ub.90.hpd,
                           rv2.ub.90.hpd,
                           rv3.ub.90.hpd,
                           var1.ub.90.hpd,
                           var2.ub.90.hpd)

ci_90.hpd = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  ci_90.hpd[j] <- mean(ifelse((lower_bound_90.hpd[,j] < beta_true[j]) &
                                (upper_bound_90.hpd[,j] > beta_true[j]),1,0 ))
}
ci_90.hpd


### Coverage of CI 90 ETI
lower_bound_90.et = cbind(int.lb.90,
                          rv1.lb.90,
                          rv2.lb.90,
                          rv3.lb.90,
                          var1.lb.90,
                          var2.lb.90)

upper_bound_90.et = cbind(int.ub.90,
                          rv1.ub.90,
                          rv2.ub.90,
                          rv3.ub.90,
                          var1.ub.90,
                          var2.ub.90)

ci_90.et = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  ci_90.et[j] <- mean(ifelse((lower_bound_90.et[,j] < beta_true[j]) &
                               (upper_bound_90.et[,j] > beta_true[j]),1,0 ))
}
# ci_90 = colMeans(ci_table_90)
ci_90.et
# 1.000 0.896 0.912 0.902

### Coverage of CI 95 HPD
lower_bound_95.hpd = cbind(int.lb.95.hpd,
                           rv1.lb.95.hpd,
                           rv2.lb.95.hpd,
                           rv3.lb.95.hpd,
                           var1.lb.95.hpd,
                           var2.lb.95.hpd)

upper_bound_95.hpd = cbind(int.ub.95.hpd,
                           rv1.ub.95.hpd,
                           rv2.ub.95.hpd,
                           rv3.ub.95.hpd,
                           var1.ub.95.hpd,
                           var2.ub.95.hpd)

ci_95.hpd = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  ci_95.hpd[j] <- mean(ifelse((lower_bound_95.hpd[,j] < beta_true[j]) &
                                (upper_bound_95.hpd[,j] > beta_true[j]),1,0 ))
}
ci_95.hpd

### Coverage of CI 95 ETI
lower_bound_95.et = cbind(int.lb.95,
                          rv1.lb.95,
                          rv2.lb.95,
                          rv3.lb.95,
                          var1.lb.95,
                          var2.lb.95)

upper_bound_95.et = cbind(int.ub.95,
                          rv1.ub.95,
                          rv2.ub.95,
                          rv3.ub.95,
                          var1.ub.95,
                          var2.ub.95)

ci_95.et = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  ci_95.et[j] <- mean(ifelse((lower_bound_95.et[,j] < beta_true[j]) &
                               (upper_bound_95.et[,j] > beta_true[j]),1,0 ))
}
ci_95.et

datafr <- cbind(bias, ese, mse, ci_90.et*100, ci_95.et*100,
                ci_90.hpd*100, ci_95.hpd*100)
row.names(datafr) <- c("intercept", "beta1", "beta2", "beta3", "var1", "var2")
colnames(datafr) <- c("Bias", "ESE", "MSE", "CP 90 ETI", "CP 95 ETI",
                      "CP 90 HPD", "CP 95 HPD")
round(datafr,3)

# Bias   ESE   MSE CP 90 ETI CP 95 ETI CP 90 HPD CP 95 HPD
# intercept  0.062 0.059 0.086     100.0     100.0     100.0     100.0
# beta1     -0.003 0.135 0.135      86.6      93.6      86.2      92.6
# beta2      0.001 0.038 0.038      88.2      93.0      87.6      93.6
# beta3      0.005 0.082 0.082      90.2      94.6      90.0      94.4
# var1       0.344 0.124 0.365      88.4      96.4      96.8      99.2
# var2       0.178 0.229 0.290      92.4      96.4      94.0      98.2

# INLA

library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(INLA)

# inv <- function(x) 1/x
# write(c("num", "intercept", "RV1", "RV2", "RV3", 
#         "int.sd", "RV1.sd", "RV2.sd", "RV3.sd",
#         "int.lb.95", "int.lb.90", "int.med", "int.ub.90", "int.ub.95",
#         "rv1.lb.95", "rv1.lb.90", "rv1.med", "rv1.ub.90", "rv1.ub.95",
#         "rv2.lb.95", "rv2.lb.90", "rv2.med", "rv2.ub.90", "rv2.ub.95",
#         "rv3.lb.95", "rv3.lb.90", "rv3.med", "rv3.ub.90", "rv3.ub.95",
#         "var1.mean",
#         "var1.lb.95", "var1.lb.90", "var1.med", "var1.ub.90", "var1.ub.95",
#         "var1.mode",
#         "var1.sd",
#         "var2.mean",
#         "var2.lb.95", "var2.lb.90", "var2.med", "var2.ub.90", "var2.ub.95",
#         "var2.mode",
#         "var2.sd",
#         "model.time",
#         "tot.time"),"results_inla_simulation_study3.txt", 
#       ncolumns=100, 
#       append=T)
# 
# for (i in 1:run){
#   #generate data
#   set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
#   data = generate_data(fix.effect = fixEffect, 
#                        random.intercept = randomIntercept, 
#                        random.slope = randomSlope,
#                        num.unit = numUnit, 
#                        obs.per.unit = obsPerUnit)
#   data2 = data.frame(data)
#   colnames(data2)[1:6] = c("group", "y", "RV1", "RV2", "RV3", "RV4")
#   data2$slopeid = data2$group + numUnit
#   
#   tic <- proc.time()
#   formula = y ~ RV1 + RV2 + RV3 + f(group, model ="iid") + 
#     f(slopeid, RV4, model ="iid")
#   model.inla <- inla(formula,
#                      data = data2, family = "poisson")
#   model.time <- proc.time() - tic
#   model.time <- model.time[3]
#   
#   ### get the credible interval of fixed effects:
#   marg.int <- model.inla$marginals.fixed$`(Intercept)`
#   int.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975),marg.int)
#   marg.rv1 <- model.inla$marginals.fixed$RV1
#   marg.rv2 <- model.inla$marginals.fixed$RV2
#   marg.rv3 <- model.inla$marginals.fixed$RV3
#   rv1.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975),marg.rv1)
#   rv2.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975),marg.rv2)
#   rv3.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975),marg.rv3)
#   inv <- function(x) 1/x
#   prec1 <- model.inla$marginals.hyperpar$`Precision for group`
#   marg.var1 <- inla.tmarginal(inv, prec1)
#   var1.mean <- inla.emarginal(function(x) x, marg.var1)
#   var1.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975), marg.var1)
#   var1.mode <- inla.mmarginal(marg.var1)
#   mm1 <- inla.emarginal(function(x) x^2, marg.var1)
#   var1.sd <- sqrt(mm1 - var1.mean^2)
#   
#   prec2 <- model.inla$marginals.hyperpar$`Precision for slopeid`
#   marg.var2 <- inla.tmarginal(inv, prec2)
#   var2.mean <- inla.emarginal(function(x) x, marg.var2)
#   var2.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975), marg.var2)
#   var2.mode <- inla.mmarginal(marg.var2)
#   mm2 <- inla.emarginal(function(x) x^2, marg.var2)
#   var2.sd <- sqrt(mm2 - var2.mean^2)
#   tot.time <- proc.time() - tic
#   tot.time <- tot.time[3]
#   
#   #write to file
#   write(c(i, model.inla$summary.fixed[,1], 
#           model.inla$summary.fixed[,2], 
#           int.quant,
#           rv1.quant,
#           rv2.quant,
#           rv3.quant,
#           var1.mean,
#           var1.quant,
#           var1.mode,
#           var1.sd,
#           var2.mean,
#           var2.quant,
#           var2.mode,
#           var2.sd,
#           model.time,
#           tot.time),"results_inla_simulation_study3.txt", 
#         ncolumns=100, 
#         append=T)
# }


# dat1<-read.table("results_inla_simulation_study3.txt",header=T) #load results
# dim(dat1)
# attach(dat1)

write(c("num", "intercept", "RV1", "RV2", "RV3", 
        "int.sd", "RV1.sd", "RV2.sd", "RV3.sd",
        "int.lb.95", "int.lb.90", "int.med", "int.ub.90", "int.ub.95",
        "rv1.lb.95", "rv1.lb.90", "rv1.med", "rv1.ub.90", "rv1.ub.95",
        "rv2.lb.95", "rv2.lb.90", "rv2.med", "rv2.ub.90", "rv2.ub.95",
        "rv3.lb.95", "rv3.lb.90", "rv3.med", "rv3.ub.90", "rv3.ub.95",
        "var1.mean",
        "var1.lb.95", "var1.lb.90", "var1.med", "var1.ub.90", "var1.ub.95",
        "var1.mode",
        "var1.sd",
        "var2.mean",
        "var2.lb.95", "var2.lb.90", "var2.med", "var2.ub.90", "var2.ub.95",
        "var2.mode",
        "var2.sd",
        "model.time",
        "tot.time"),"results_inla_simulation_study3_set_prior.txt", 
      ncolumns=100, 
      append=T)

for (i in 1:run){
  #generate data
  set.seed(seeds[i]) #sets the seed to the ith element in the vector seeds
  data = generate_data(fix.effect = fixEffect, 
                       random.intercept = randomIntercept, 
                       random.slope = randomSlope,
                       num.unit = numUnit, 
                       obs.per.unit = obsPerUnit)
  data2 = data.frame(data)
  colnames(data2)[1:6] = c("group", "y", "RV1", "RV2", "RV3", "RV4")
  data2$slopeid = data2$group + numUnit
  prior.fixed <- list(mean.intercept = 0, prec.intercept = 0.0001,
                      mean = 0, prec = 0.0001)
  prec.prior <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01),
                                 initial = 4, fixed = FALSE))
  tic <- proc.time()
  formula = y ~ RV1 + RV2 + RV3 + 
    f(group, model ="iid", hyper = prec.prior) + 
    f(slopeid, RV4, model ="iid", hyper = prec.prior)
  model.inla <- inla(formula,
                     data = data2, family = "poisson",
                     control.fixed = prior.fixed)
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
  prec1 <- model.inla$marginals.hyperpar$`Precision for group`
  marg.var1 <- inla.tmarginal(inv, prec1)
  var1.mean <- inla.emarginal(function(x) x, marg.var1)
  var1.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975), marg.var1)
  var1.mode <- inla.mmarginal(marg.var1)
  mm1 <- inla.emarginal(function(x) x^2, marg.var1)
  var1.sd <- sqrt(mm1 - var1.mean^2)
  
  prec2 <- model.inla$marginals.hyperpar$`Precision for slopeid`
  marg.var2 <- inla.tmarginal(inv, prec2)
  var2.mean <- inla.emarginal(function(x) x, marg.var2)
  var2.quant <- inla.qmarginal(c(0.025, 0.05, 0.5, 0.95, 0.975), marg.var2)
  var2.mode <- inla.mmarginal(marg.var2)
  mm2 <- inla.emarginal(function(x) x^2, marg.var2)
  var2.sd <- sqrt(mm2 - var2.mean^2)
  tot.time <- proc.time() - tic
  tot.time <- tot.time[3]
  
  #write to file
  write(c(i, model.inla$summary.fixed[,1], 
          model.inla$summary.fixed[,2], 
          int.quant,
          rv1.quant,
          rv2.quant,
          rv3.quant,
          var1.mean,
          var1.quant,
          var1.mode,
          var1.sd,
          var2.mean,
          var2.quant,
          var2.mode,
          var2.sd,
          model.time,
          tot.time),"results_inla_simulation_study3_set_prior.txt", 
        ncolumns=100, 
        append=T)
  print(i)
}

dat1<-read.table("results_inla_simulation_study3_set_prior.txt",header=T) #load results


beta = dat1[, c(2:5)]
beta = cbind(beta, dat1$var1.mean, dat1$var2.mean)

model_time = quantile(dat1$model.time, 0.5)
model_time
# 1.627
tot_time = quantile(dat1$tot.time, 0.5)
tot_time
#1.773

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

### Coverage of CI 90 ETI
lower_bound_90.et = cbind(dat1$int.lb.90,
                          dat1$rv1.lb.90,
                          dat1$rv2.lb.90,
                          dat1$rv3.lb.90,
                          dat1$var1.lb.90,
                          dat1$var2.lb.90)

upper_bound_90.et = cbind(dat1$int.ub.90,
                          dat1$rv1.ub.90,
                          dat1$rv2.ub.90,
                          dat1$rv3.ub.90,
                          dat1$var1.ub.90,
                          dat1$var2.ub.90)

ci_90.et = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  ci_90.et[j] <- mean(ifelse((lower_bound_90.et[,j] < beta_true[j]) & 
                               (upper_bound_90.et[,j] > beta_true[j]),1,0 ))
}
# ci_90 = colMeans(ci_table_90)
ci_90.et


## Coverage of CI 95 ETI
lower_bound_95.et = cbind(dat1$int.lb.95,
                          dat1$rv1.lb.95,
                          dat1$rv2.lb.95,
                          dat1$rv3.lb.95,
                          dat1$var1.lb.95,
                          dat1$var2.lb.95)

upper_bound_95.et = cbind(dat1$int.ub.95,
                          dat1$rv1.ub.95,
                          dat1$rv2.ub.95,
                          dat1$rv3.ub.95,
                          dat1$var1.ub.95,
                          dat1$var2.ub.95)

ci_95.et = rep(0,length(beta_true))
for (j in 1:length(beta_true)){
  ci_95.et[j] <- mean(ifelse((lower_bound_95.et[,j] < beta_true[j]) & 
                               (upper_bound_95.et[,j] > beta_true[j]),1,0 ))
}
ci_95.et


datafr <- cbind(bias, ese, mse, ci_90.et*100, ci_95.et*100)
row.names(datafr) <- c("intercept", "beta1", "beta2", "beta3", "var1", "var2")
colnames(datafr) <- c("Bias", "ESE", "MSE", "Cov 90", "Cov 95")

datafr2 = data.frame(datafr)
round(datafr2, 3)

# Bias   ESE   MSE Cov.90 Cov.95
# intercept  0.074 0.056 0.092  100.0  100.0
# beta1     -0.006 0.136 0.136   85.8   92.4
# beta2      0.000 0.039 0.039   86.4   92.8
# beta3      0.004 0.083 0.083   89.4   93.6
# var1       0.262 0.114 0.286   96.8   99.8
# var2      -0.005 0.266 0.266   82.2   84.8

# set prior
# Bias   ESE   MSE Cov.90 Cov.95
# intercept  0.064 0.056 0.085  100.0  100.0
# beta1     -0.003 0.135 0.134   86.4   93.2
# beta2      0.001 0.038 0.038   88.6   93.6
# beta3      0.005 0.081 0.082   90.4   94.6
# var1       0.333 0.122 0.354   90.4   97.6
# var2       0.162 0.222 0.275   91.8   95.8
