rm(list=ls())
library(HSAUR3)
library(microbenchmark)
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
source("labpmm_ci.R")
source("labpmm_map.R")
source("labpmm_grid.R")

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
x6 = ifelse(period=="4", 1, 0)

y=seizure.rate
t = cbind(x2,x3,x4,x5,x6)

#LABPMM
tic()
model1 <- labpmm_ci(y=seizure.rate, t = t, 
                    numUnit=numUnit, obsPerUnit = obsPerUnit, 
                 gammaA = 0.01, gammaB = 0.01, 
                 numFixEff = numFixEff, fixEffVar = 10^4)
toc()
tic()
model1 <- labpmm(y=seizure.rate, t = t, 
                    numUnit=numUnit, obsPerUnit = obsPerUnit, 
                    gammaA = 0.01, gammaB = 0.01, 
                    numFixEff = numFixEff, fixEffVar = 10^4)
toc()
tic()
model2 <- labpmm_grid(y=seizure.rate, t = t, 
                      numUnit=numUnit, obsPerUnit = obsPerUnit, 
                    gammaA = 0.01, gammaB = 0.01, 
                    numFixEff = numFixEff, fixEffVar = 10^4)
toc()
#Likelihood
library(Matrix)
library(lme4)
library(tictoc)
library(boot)
tic()
model3 <- glmer(seizure.rate ~ x2+x3+x4+x5+x6 +
                          (1|subject), family="poisson")
toc()
tic()
llh_ci_95 <- confint.merMod(model3, method = "profile")
toc()

#MCMC
library(rjags)
library(coda)
library(readr)
library(runjags)

cat("model
  {
    for (j in 1:J){
      b[j] ~ dnorm(0, tau)
    }
    for (i in 1:N){
      rate[i] ~ dpois(phi[i])
      log(phi[i]) <- beta0 + beta1*base[i] + beta2*treatment[i] + beta3*base_treatment[i]
                  + beta4*age[i] + beta5*fourth_visit[i] + b[subject[i]]
    }
    tau ~ dgamma(0.01, 0.01)
    variance = 1/tau
    beta0 ~ dnorm(0,1.0E-4)
    beta1 ~ dnorm(0,1.0E-4)
    beta2 ~ dnorm(0,1.0E-4)
    beta3 ~ dnorm(0,1.0E-4)
    beta4 ~ dnorm(0,1.0E-4)
    beta5 ~ dnorm(0,1.0E-4)
  }", file="epilepsy.model.txt")

# Prepare data:

N <- length(y)
J <- n_distinct(subject)
t2 <- data.frame(t)
colnames(t2)[1:5] = c("base", "treatment", "base_treatment", "age", "fourth_visit")
my.data <- list(rate=y, base = t2$base, treatment = t2$treatment, 
                base_treatment = t2$base_treatment, age = t2$age, 
                fourth_visit = t2$fourth_visit, N = N, J = J, subject = subject)

# Initial parameters: 3 lists for 3 chains
my.inits <- list(
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 111, beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, beta4=0.1,beta5=0.2,tau=1),
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 112, beta0=0.1,beta1=0.1, beta2=0.1,beta3=0.1, beta4=0.1,beta5=0.1, tau=1),
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 113, beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, beta4=0.1,beta5=0.2,tau=0.1)
)

parameters <- c("beta0","beta1","beta2","beta3","beta4","beta5","tau","variance")

## Running JAGS:
tic()
jags <- jags.model(file="epilepsy.model.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 3)
update(jags,5000)
epi.sim <- coda.samples(model = jags,
                        variable.names = parameters,
                        n.iter=10000, 
                        thin=10)
toc()
summary(epi.sim)
epi.mcmc <- as.mcmc.list(epi.sim)
epi.combined <- combine.mcmc(epi.mcmc)
HPDinterval(epi.combined)
quantile(epi.combined[,7], c(0.025, 0.975))
# for (i in 1:numFixEff+1){
#   quantile(epi.combined[,i], c(0.025, 0.975))
# }

# INLA
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(INLA)

data2 <- data.frame(cbind(t, subject, seizure.rate))
colnames(data2) = c("base1", "treatment1", "base_treatment1", 
                     "age1", "fourth_visit1", "subject1", "seizure")
inv <- function(x) 1/x
prior.fixed <- list(mean.intercept = 0, prec.intercept = 0.0001,
                    mean = 0, prec = 0.0001)
prec.prior <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01),
                               initial = 4, fixed = FALSE))
formula = seizure ~ base1 + treatment1 + 
  base_treatment1 + age1 + fourth_visit1 + 
  f(subject1, model ="iid", hyper = prec.prior) 
# formula = seizure.rate ~ x2+x3+x4+x5+x6 + 
#   f(subject, model ="iid", hyper = prec.prior) 
tic()
model.inla <- inla(formula,
                   data = data2, family = "poisson",
                   control.fixed = prior.fixed)
toc()
summary(model.inla)

prec <- model.inla$marginals.hyperpar$`Precision for subject`
marg.var <- inla.tmarginal(inv, prec)
var.inla.mean <- inla.emarginal(function(x) x, marg.var)
var.inla.quant <- inla.qmarginal(c(0.025, 0.975), marg.var)
var.mode <- inla.mmarginal(marg.var)
mm <- inla.emarginal(function(x) x^2, marg.var)
var.inla.sd <- sqrt(mm - var.inla.mean^2)


#####
lapmm1_fixeff <- model1$fix_effect
lapmm1_sd_fe <- model1$sd_fixEffect
lapmm2_fixeff <- model2$fix_effect
lapmm2_sd_fe <- model2$sd_fixEffect
llh_fixeff <- fixef(model3)
vcov.matrix = vcov(model3)
llh_sd_fe <- sqrt(vcov.matrix[row(vcov.matrix) == col(vcov.matrix)])
jag_model <- summary(epi.combined)
jag_fixeff <- jag_model$statistics[,1][1:6]
jag_sd_fe <- jag_model$statistics[,2][1:6]
datafr <- cbind(lapmm1_fixeff, lapmm1_sd_fe,lapmm2_fixeff, 
                lapmm2_sd_fe, llh_fixeff, llh_sd_fe, jag_fixeff,
                jag_sd_fe)
row.names(datafr) <- c("intercept", "baseline", "treatment", "treatbase", 
                       "age", "V4")
round(datafr,3)

round(lapmm1_fixeff - 1.96* lapmm1_sd_fe, 3)
round(lapmm1_fixeff + 1.96* lapmm1_sd_fe, 3)

round(lapmm2_fixeff - 1.96* lapmm2_sd_fe, 3)
round(lapmm2_fixeff + 1.96* lapmm2_sd_fe, 3)

round(llh_fixeff - 1.96* llh_sd_fe, 3)
round(llh_fixeff + 1.96* llh_sd_fe, 3)

round(quantile(epi.combined[,6], c(0.025, 0.975)),3)

var.inla.quant
var.inla.mean
var.inla.sd


#### Plot
par(mfrow = c(1,1))
xaxis = seq(-3,3, by=0.001)
xlabel = c(expression(beta[0]), expression(beta[1]), expression(beta[2]), 
           expression(beta[3]), expression(beta[4]), expression(beta[5]))
for (i in 1:length(xlabel)){
  plot(epi.combined[,i], trace = FALSE, density=TRUE,
       type = "l", 
       xlab = xlabel[i],
       ylab = "posterior",
       main = "")
  lines(xaxis, dnorm(xaxis, model2$fix_effect[i], model2$sd_fixEffect[i]), 
        col = "blue",
        type = "l")
  lines(xaxis, dnorm(xaxis, model1$fix_effect[i], model1$sd_fixEffect[i]), 
        col = "red", type = "l", lty = 2)
  lines(model.inla$marginals.fixed[[i]], col = 'green', type = 'l', lty = 2)
  abline(v = fixef(model3)[i], lty = 2)
  legend("topright", legend=c("LABPMM-grid", 
                              "LABPMM-MAP",
                              "INLA",
                              "MCMC"),
         col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)
}

plot(epi.combined[,6], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0,8),
     xlab = expression(beta[5]),
     ylab = "posterior",
     main = "")
lines(xaxis, dnorm(xaxis, model2$fix_effect[6], model2$sd_fixEffect[6]), 
      col = "blue",
      type = "l")
lines(xaxis, dnorm(xaxis, model1$fix_effect[6], model1$sd_fixEffect[6]), 
      col = "red", type = "l", lty = 2)
lines(model.inla$marginals.fixed[[6]], col = 'green', type = 'l', lty = 2)
abline(v = fixef(model3)[6], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)

inv <- function(x) 1/x
prec <- model.inla$marginals.hyperpar$`Precision for subject`
marg.var <- inla.tmarginal(inv, prec)

plot(epi.combined[,8], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0,6),
     xlab = expression(sigma^2),
     ylab = "posterior",
     main = "")
set.seed(12)
y2 <- 1/exp(rnorm(xaxis, model2$mean_v, model2$sd_v))
lines(density(y2), 
      col = "blue",
      type = "l")
lines(marg.var, col = 'green', type = 'l', lty = 2)
abline(v = as.data.frame(VarCorr(model3))$vcov, lty = 2)
legend("topright", legend=c("LABPMM",
                            "INLA",
                            "MCMC"),
       col=c("blue",  "green", "black"), lty=c(1,2,2,1), cex=0.8)


## Microbenchmark
model1_time <- function() labpmm_ci(y=seizure.rate, t = t, 
                                    numUnit=numUnit, obsPerUnit = obsPerUnit, 
                                    gammaA = 0.01, gammaB = 0.01, 
                                    numFixEff = numFixEff, fixEffVar = 10^4)
model2_time <- function() labpmm_grid(y=seizure.rate, t = t, 
                                      numUnit=numUnit, obsPerUnit = obsPerUnit, 
                                      gammaA = 0.01, gammaB = 0.01, 
                                      numFixEff = numFixEff, fixEffVar = 10^4)
model1a_time <- function() labpmm(y=seizure.rate, t = t, 
                                  numUnit=numUnit, obsPerUnit = obsPerUnit, 
                                  gammaA = 0.01, gammaB = 0.01, 
                                  numFixEff = numFixEff, fixEffVar = 10^4)
model3_time <- function() glmer(seizure.rate ~ x2+x3+x4+x5+x6 +
                                  (1|subject), family="poisson")

inla_time <- function() inla(formula,
                             data = data2, family = "poisson",
                             control.fixed = prior.fixed)

microbenchmark(model1a_time(), model1_time(), model2_time(),
               model3_time(), inla_time())
######## Model with one random intercept and one random slope

rm(list=ls())
library(HSAUR3)
data = epilepsy
attach(epilepsy)
source("labpmm2_ci.R")
source("labpmm2_map.R")
source("labpmm2_grid.R")

n = nrow(data)
numUnit = 59
obsPerUnit = 4
numFixEff = 6

x1 = subject
x2 = log(base/4)
x2 = x2-mean(x2)
x3 = ifelse(treatment=="placebo", 0, 1)
x4 = log(base/4)*x3
x4 = x4-mean(x4)
x5 = log(age)
x5 = (x5-mean(x5))
x6 = ifelse(period=="4", 1, 0)
x6 = (as.numeric(period)-2.5)*2/10

y=seizure.rate
t = cbind(x2,x3,x4,x5,x6)

#LABPMM
tic()
model4 <- labpmm2_ci(y=seizure.rate, t = t, r = x6,
                     numUnit=numUnit, obsPerUnit = obsPerUnit, 
                     gammaA = 0.01, gammaB = 0.01, 
                     fixEffVar = 10^4)
toc()
tic()
model4a <- labpmm2(y=seizure.rate, t = t, r = x6,
                   numUnit=numUnit, obsPerUnit = obsPerUnit, 
                   gammaA = 0.01, gammaB = 0.01, 
                   fixEffVar = 10^4)
toc()
tic()
model5 <- labpmm2_grid(y=seizure.rate, t = t, r = x6,
                       numUnit=numUnit, obsPerUnit = obsPerUnit, 
                       gammaA = 0.01, gammaB = 0.01, 
                       fixEffVar = 10^4)
toc()
#Likelihood
library(Matrix)
library(lme4)
library(tictoc)
library(boot)
tic()
model6 <- glmer(seizure.rate ~ x2+x3+x4+x5+x6 +(1+x6|subject), 
                family="poisson")
toc()
tic()
llh.ci.95.model2 <- confint.merMod(model6, method = "profile", level = 0.95)
toc()
# bootstrap CI and sd
# lpmm <- function(formula, data, indices) {
#   d <- data[indices,] # allows boot to select sample
#   fit <- glmer(formula, data=d, family="poisson")
#   return(as.data.frame(VarCorr(fit))$vcov)
# }
# data3 = cbind(seizure.rate, t, subject)
# data3 <- data.frame(data3)
# 
# set.seed(1)
# # bootstrapping with 1000 replications
# 
# boot.res <- boot(data=data3,
#                 statistic=lpmm,
#                 R=5000,
#                 formula=seizure.rate ~ x2+x3+x4+x5+x6+(1+x6|subject))
# boot.ci(boot.res, type = "all")
#MCMC
library(rjags)
library(coda)
library(readr)
library(runjags)
cat("model
  {
    for (j in 1:J){
      b0[j] ~ dnorm(0, tau1)
      b1[j] ~ dnorm(0, tau2)
    }
    for (i in 1:N){
      rate[i] ~ dpois(phi[i])
      log(phi[i]) <- beta0 + beta1*base[i] + beta2*treatment[i] + beta3*base_treatment[i]
                  + beta4*age[i] + beta5*fourth_visit[i] + b0[subject[i]] + b1[subject[i]]*fourth_visit[i]
    }
    tau1 ~ dgamma(0.01, 0.01)
    tau2 ~ dgamma(0.01, 0.01)
    var1 = 1/tau1
    var2 = 1/tau2
    beta0 ~ dnorm(0,1.0E-4)
    beta1 ~ dnorm(0,1.0E-4)
    beta2 ~ dnorm(0,1.0E-4)
    beta3 ~ dnorm(0,1.0E-4)
    beta4 ~ dnorm(0,1.0E-4)
    beta5 ~ dnorm(0,1.0E-4)
  }", file="epilepsy_random_slope_gamma.model.txt")

# Prepare data:

N <- length(y)
J <- n_distinct(subject)
t <- data.frame(t)
colnames(t)[1:5] = c("base", "treatment", "base_treatment", "age", "fourth_visit")
my.data <- list(rate=y, base = t$base, treatment = t$treatment, 
                base_treatment = t$base_treatment, age = t$age, 
                fourth_visit = t$fourth_visit, N = N, J = J, 
                subject = subject)

my.inits <- list(
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 111, beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, beta4=0.1,beta5=0.2,tau1=1, tau2=2),
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 112, beta0=0.1,beta1=0.1, beta2=0.1,beta3=0.1, beta4=0.1,beta5=0.1, tau1=1, tau2=2),
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 113, beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, beta4=0.1,beta5=0.2,tau1=0.1, tau2=2)
)
parameters <- c("beta0","beta1","beta2","beta3","beta4","beta5","tau1", "tau2", "var1", "var2")
tic()
jags2 <- jags.model(file="epilepsy_random_slope_gamma.model.txt",
                    data = my.data,
                    inits = my.inits,
                    n.chains = 3)
update(jags2,10000)
epi.sim2 <- coda.samples(model = jags2,
                         variable.names = parameters,
                         n.iter=20000, 
                         thin=10)
toc()
summary(epi.sim2)
epi.mcmc2 <- as.mcmc.list(epi.sim2)
gelman.diag(epi.mcmc2)
epi.combined2 <- combine.mcmc(epi.mcmc2)
HPDinterval(epi.combined2)
quantile(epi.combined2[,7], c(0.025, 0.975))

###INLA
library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(INLA)
data2 <- data.frame(cbind(t, subject))
colnames(data2) = c("base1", "treatment1", "base_treatment1", 
                    "age1", "fourth_visit1", "subject1")

subject.id <- as.numeric(data2$subject)
data2$slopeid = subject.id + numUnit
prior.fixed <- list(mean.intercept = 0, prec.intercept = 0.0001,
                    mean = 0, prec = 0.0001)
prec.prior <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01),
                               initial = 4, fixed = FALSE))
formula = seizure.rate ~ base1 + treatment1 + 
  base_treatment1 + age1 + fourth_visit1 + 
  f(subject1, model ="iid", hyper = prec.prior) + 
  f(slopeid, fourth_visit1, model ="iid", hyper = prec.prior)
tic()
model.inla <- inla(formula,
                   data = data2, family = "poisson",
                   control.fixed = prior.fixed)
toc()
summary(model.inla)

inv <- function(x) 1/x
prec1 <- model.inla$marginals.hyperpar$`Precision for subject1`
marg.var1 <- inla.tmarginal(inv, prec1)
var1.inla.mean <- inla.emarginal(function(x) x, marg.var1)
var1.inla.quant <- inla.qmarginal(c(0.025, 0.975), marg.var1)
var1.inla.mode <- inla.mmarginal(marg.var1)
mm1 <- inla.emarginal(function(x) x^2, marg.var1)
var1.inla.sd <- sqrt(mm1 - var1.inla.mean^2)

prec2 <- model.inla$marginals.hyperpar$`Precision for slopeid`
marg.var2 <- inla.tmarginal(inv, prec2)
var2.inla.mean <- inla.emarginal(function(x) x, marg.var2)
var2.inla.quant <- inla.qmarginal(c(0.025, 0.975), marg.var2)
var2.inla.mode <- inla.mmarginal(marg.var2)
mm2 <- inla.emarginal(function(x) x^2, marg.var2)
var2.inla.sd <- sqrt(mm2 - var2.inla.mean^2)

#####
lapmm5_fixeff <- model5$fix_effect
lapmm5_sd_fe <- model5$sd_fixEffect
lapmm4_fixeff <- model4$fix_effect
lapmm4_sd_fe <- model4$sd_fixEffect
llh_fixeff <- fixef(model6)
vcov.matrix = vcov(model6)
llh_sd_fe <- sqrt(vcov.matrix[row(vcov.matrix) == col(vcov.matrix)])
jag_model <- summary(epi.combined2)
jag_fixeff <- jag_model$statistics[,1][1:6]
jag_sd_fe <- jag_model$statistics[,2][1:6]
datafr <- cbind(lapmm5_fixeff, lapmm5_sd_fe,lapmm4_fixeff, 
                lapmm4_sd_fe, llh_fixeff, llh_sd_fe, jag_fixeff,
                jag_sd_fe)
row.names(datafr) <- c("intercept", "baseline", "treatment", "treatbase", 
                       "age", "V4")
round(datafr,3)

round(lapmm5_fixeff - 1.96* lapmm5_sd_fe, 3)
round(lapmm5_fixeff + 1.96* lapmm5_sd_fe, 3)

round(lapmm4_fixeff - 1.96* lapmm4_sd_fe, 3)
round(lapmm4_fixeff + 1.96* lapmm4_sd_fe, 3)

round(llh_fixeff - 1.96* llh_sd_fe, 3)
round(llh_fixeff + 1.96* llh_sd_fe, 3)

round(quantile(epi.combined2[,10], c(0.025, 0.975)),3)

var1.inla.quant
var1.inla.mean
var1.inla.sd
var1.inla.quant
var1.inla.mean
var1.inla.sd

####Model II

par(mfrow = c(1,1))
xaxis = seq(-2,3, by=0.001)
xlabel = c(expression(beta[0]), expression(beta[1]), expression(beta[2]), 
           expression(beta[3]), expression(beta[4]), expression(beta[5]))
for (i in 1:length(xlabel)){
  plot(epi.combined2[,i], trace = FALSE, density=TRUE,
       type = "l", 
       xlab = xlabel[i],
       ylab = "posterior",
       main = "")
  lines(xaxis, dnorm(xaxis, model5$fix_effect[i], model5$sd_fixEffect[i]), 
        col = "blue",
        type = "l")
  lines(xaxis, dnorm(xaxis, model4$fix_effect[i], model4$sd_fixEffect[i]), 
        col = "red", type = "l", lty = 2)
  lines(model.inla$marginals.fixed[[i]], col = 'green', type = 'l', lty = 2)
  abline(v = fixef(model6)[i], lty = 2)
  legend("topright", legend=c("LABPMM-grid", 
                              "LABPMM-MAP",
                              "INLA",
                              "MCMC"),
         col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)
}

plot(epi.combined2[,2], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0,3),
     xlab = expression(beta[1]),
     ylab = "posterior",
     main = "")
lines(xaxis, dnorm(xaxis, model5$fix_effect[2], model5$sd_fixEffect[2]), 
      col = "blue",
      type = "l")
lines(xaxis, dnorm(xaxis, model4$fix_effect[2], model4$sd_fixEffect[2]), 
      col = "red", type = "l", lty = 2)
lines(model.inla$marginals.fixed[[2]], col = 'green', type = 'l', lty = 2)
abline(v = fixef(model6)[2], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)

xaxis = seq(-1, 3, by=0.001)
plot(epi.combined2[,9], trace = FALSE, density=TRUE,
     type = "l", 
     xlab = expression(sigma[0]^2),
     ylab = "posterior",
     main = "")
set.seed(12)
y5 <- 1/exp(rnorm(xaxis, model5$mean_v[1], model5$sd_v[1]))
lines(density(y5), 
      col = "blue",
      type = "l")
set.seed(12)
y4 <- 1/exp(rnorm(xaxis, model4$mean_v[1], model4$sd_v[1]))
lines(density(y4),
      col = "red", type = "l", lty = 2)
lines(marg.var1, col = 'green', type = 'l', lty = 2)
abline(v = as.data.frame(VarCorr(model6))$vcov[1], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)

plot(epi.combined2[,10], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0,2),
     xlim = c(0,2),
     xlab = expression(sigma[1]^2),
     ylab = "posterior",
     main = "")
set.seed(12)
y5 <- 1/exp(rnorm(xaxis, model5$mean_v[2], model5$sd_v[2]))
lines(density(y5), 
      col = "blue",
      type = "l")
set.seed(12)
y4 <- 1/exp(rnorm(xaxis, model4$mean_v[2], model4$sd_v[2]))
lines(density(y4),
      col = "red", type = "l", lty = 2)
lines(marg.var2, col = 'green', type = 'l', lty = 2)
abline(v = as.data.frame(VarCorr(model6))$vcov[2], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)

model4_time <- function() labpmm2_ci(y=seizure.rate, t = t, r = x6,
                                     numUnit=numUnit, obsPerUnit = obsPerUnit, 
                                     gammaA = 0.01, gammaB = 0.01, 
                                     fixEffVar = 10^4)
model5_time <- function() labpmm2_grid(y=seizure.rate, t = t, r = x6,
                                      numUnit=numUnit, obsPerUnit = obsPerUnit, 
                                      gammaA = 0.01, gammaB = 0.01, 
                                      fixEffVar = 10^4)
model4a_time <- function() labpmm2(y=seizure.rate, t = t, r = x6,
                                  numUnit=numUnit, obsPerUnit = obsPerUnit, 
                                  gammaA = 0.01, gammaB = 0.01, 
                                  fixEffVar = 10^4)
model6_time <- function() glmer(seizure.rate ~ x2+x3+x4+x5+x6 +(1+x6|subject), 
                                family="poisson")

inla_time <- function() inla(formula,
                             data = data2, family = "poisson",
                             control.fixed = prior.fixed)

microbenchmark(model4a_time(), model4_time(), model5_time(),
               model6_time(), inla_time())

### EDA
rm(list = ls())
library(HSAUR3)

#load data
data = epilepsy
attach(epilepsy)
### Spaghetti plot over time
library(ggplot2)

### 
# spaghetti plot
ggplot(data = data, aes(x = period, y = seizure.rate, color = treatment)) +       
  geom_line(aes(group = subject)) + geom_point() + labs(x = "Period", y = "Seizure rate",
                                                        color = "Treatment")
# mean plot
library(dplyr)
data %>% group_by(treatment, period) %>% summarise(ave_rate = mean(seizure.rate)) %>% 
  ggplot(aes(x = period, y = ave_rate, color = treatment)) +       
  geom_line(aes(group = treatment)) + geom_point() + labs(x = "Period", y = "Average seizure rate",
                                                          color = "Treatment")



ggplot(data = data, aes(x = period, y = seizure.rate)) +       
  geom_line(aes(group = subject)) + geom_point() + labs(x = "Period", 
                                                        y = "Seizure rate")



