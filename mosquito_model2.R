rm(list=ls())

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

#Fit model with one random intercept and one random slope
data = read.csv("mosquito_data.csv",  header = TRUE, sep = ";")
data$Time = (data$Month - mean(data$Month))/6
data$SettlementFac <- as.factor(data$Settlement)
data$TimeSetInt <- data$Time * data$Settlement
data$Settlement = data$Settlement - mean(data$Settlement)
data$TimeSetInt = data$TimeSetInt - mean(data$TimeSetInt)

source("labpmm2_ci.R")
source("labpmm2_map.R")
source("labpmm2_grid.R")

t = cbind(data$Settlement, data$Time, data$TimeSetInt)
r = data$Time
numUnit = 40
obsPerUnit = 6
numFixEff = 4
tic()                                                      
model1 <- labpmm2_ci(y=data$Total, t = t, r=r, numUnit=numUnit, obsPerUnit = obsPerUnit, 
                    gammaA = 0.01, gammaB = 0.01)
toc()
tic()                                                      
model1a <- labpmm2(y=data$Total, t = t, r=r, numUnit=numUnit, obsPerUnit = obsPerUnit, 
                  gammaA = 0.01, gammaB = 0.01)
toc()
tic()                                                      
model2 <- labpmm2_grid(y=data$Total, t = t, r=r, numUnit=numUnit, obsPerUnit = obsPerUnit, 
                      gammaA = 0.01, gammaB = 0.01)
toc()


library(Matrix)
library(lme4)
tic()
fpmm <- glmer(Total ~ Settlement +Time + TimeSetInt +(1+Time|Housecode), 
              data= data, family="poisson")
toc()
tic()
llh.ci.95 <- confint.merMod(fpmm, level = 0.95, method = "profile")
toc()

summary(fpmm)

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
      log(phi[i]) <- beta0 + beta1*Settlement[i] + beta3*TimeSetInt[i] + 
                  beta2*Time[i] + b0[Housecode[i]] + b1[Housecode[i]]*Time[i]
    }
    tau1 ~ dgamma(0.01, 0.01)
    tau2 ~ dgamma(0.01, 0.01)
    var1 = 1/tau1
    var2 = 1/tau2
    beta0 ~ dnorm(0,1.0E-4)
    beta1 ~ dnorm(0,1.0E-4)
    beta2 ~ dnorm(0,1.0E-4)
    beta3 ~ dnorm(0,1.0E-4)
  }", file="mosquito.model3.txt")

# Prepare data:

tic()
N <- length(Total)
J <- n_distinct(Housecode)
my.data <- list(rate=data$Total, Settlement = data$Settlement, 
                TimeSetInt = data$TimeSetInt, Time = data$Time,
                N = N, J = J, Housecode = as.factor(data$Housecode))

# Initial parameters: 3 lists for 3 chains

my.inits <- list(
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 137, beta0=0.1,beta1=0.2, beta2=0.1, beta3 = 0.2, tau1=1, tau2 = 0.5),
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 253, beta0=0.1,beta1=0.1, beta2=0.1, beta3 = 0.2, tau1=1, tau2 = 0.5),
  list(".RNG.name" = "base::Wichmann-Hill",
       ".RNG.seed" = 121, beta0=0.1,beta1=0.2, beta2=0.1, beta3 = 0.2, tau1=0.1, tau2 = 0.5)
)

parameters <- c("beta0","beta1","beta2", "beta3", "tau1","tau2", 
                "var1", "var2")

## Running JAGS:

jags <- jags.model(file="mosquito.model3.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 3)
update(jags,5000)
epi.sim <- coda.samples(model = jags,
                        variable.names = parameters,
                        n.iter=10000, 
                        thin=10)
toc()


epi.mcmc <- as.mcmc.list(epi.sim)
epi.combined <- combine.mcmc(epi.mcmc)
gelman.diag(epi.mcmc)
#1.
round(quantile(epi.combined[,2], c(0.025, 0.975)),3)


lapmm1_fixeff <- model1$fix_effect
lapmm1_sd_fe <- model1$sd_fixEffect
lapmm2_fixeff <- model2$fix_effect
lapmm2_sd_fe <- model2$sd_fixEffect
llh_fixeff <- fixef(fpmm)
vcov.matrix = vcov(fpmm)
llh_sd_fe <- sqrt(vcov.matrix[row(vcov.matrix) == col(vcov.matrix)])
jag_model <- summary(epi.combined)
jag_fixeff <- jag_model$statistics[,1][1:4]
jag_sd_fe <- jag_model$statistics[,2][1:4]
datafr <- cbind(lapmm1_fixeff, lapmm1_sd_fe, lapmm2_fixeff, lapmm2_sd_fe, 
                llh_fixeff, llh_sd_fe, jag_fixeff,
                jag_sd_fe)
row.names(datafr) <- c("intercept", "settlement", "timeSetInt", "Time")
round(datafr,3)

# lapmm_fixeff lapmm_sd_fe  llh_fixeff llh_sd_fe  jag_fixeff jag_sd_fe
# intercept   -0.02327572   0.2260674 -0.08096902 0.2132056 -0.08885304 0.2202063
# settlement   0.65087182   0.3135941  0.67928213 0.2926455  0.66653513 0.3069054
# timeSetInt  -1.64213742   0.2010463 -1.64213920 0.1988593 -1.64976929 0.1996114

lb_lapmm1 = lapmm1_fixeff - 1.96*lapmm1_sd_fe
ub_lapmm1 = lapmm1_fixeff + 1.96*lapmm1_sd_fe

lb_lapmm2 = lapmm2_fixeff - 1.96*lapmm1_sd_fe
ub_lapmm2 = lapmm2_fixeff + 1.96*lapmm1_sd_fe

lb_rjags = jag_model$quantiles[,1][1:4]
ub_rjags = jag_model$quantiles[,5][1:4]

lb_glmer = llh.ci.95[4:7,1]
ub_glmer = llh.ci.95[4:7,2]

datafr2 <- rbind(cbind(lapmm1_fixeff, lb_lapmm1, ub_lapmm1, lapmm1_sd_fe), 
                 cbind(lapmm2_fixeff, lb_lapmm2, ub_lapmm2, lapmm2_sd_fe), 
                 cbind(llh_fixeff, lb_glmer, ub_glmer,  llh_sd_fe), 
                 cbind(jag_fixeff, lb_rjags, ub_rjags, jag_sd_fe))
datafr2 = data.frame(datafr2)
datafr2 %>% mutate(across(where(is.numeric), round, digits=3))


library(Matrix)
library(foreach)
library(parallel)
library(sp)
library(INLA)

data$Houseid <- as.factor(data$Housecode)
# dat.inla <- cbind(data$Total, data$Settlement, data$TimeSetInt, data$Houseid)
# dat.inla <- data.frame(dat.inla)
data$slopeid = as.numeric(data$Houseid) + numUnit
prior.fixed <- list(mean.intercept = 0, prec.intercept = 0.0001,
                    mean = 0, prec = 0.0001)
prec.prior <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01),
                               initial = 4, fixed = FALSE))
tic()
model.inla <- inla(Total ~ Settlement + Time + TimeSetInt + 
                     + f(Houseid, model ="iid", hyper = prec.prior) + 
                     f(slopeid, Time, model ="iid", hyper = prec.prior),
                   data = data, family = "poisson",
                   control.fixed = prior.fixed) 
toc()

summary(model.inla)

# fixed effect
# model.inla$summary.fixed[,1]
# # sd fixed effect
# model.inla$summary.fixed[,2]

inv <- function(x) 1/x
prec1 <- model.inla$marginals.hyperpar$`Precision for Houseid`
marg.var1 <- inla.tmarginal(inv, prec1)
var1.inla.mean <- inla.emarginal(function(x) x, marg.var1)
var1.inla.quant <- inla.qmarginal(c(0.025, 0.975), marg.var1)
var1.mode <- inla.mmarginal(marg.var1)
mm1 <- inla.emarginal(function(x) x^2, marg.var1)
var1.inla.sd <- sqrt(mm1 - var1.inla.mean^2)
round(var1.inla.quant,3)
var1.inla.mean
var1.inla.sd

prec2 <- model.inla$marginals.hyperpar$`Precision for slopeid`
marg.var2 <- inla.tmarginal(inv, prec2)
var2.inla.mean <- inla.emarginal(function(x) x, marg.var2)
var2.inla.quant <- inla.qmarginal(c(0.025, 0.975), marg.var2)
var2.mode <- inla.mmarginal(marg.var2)
mm2 <- inla.emarginal(function(x) x^2, marg.var2)
var2.inla.sd <- sqrt(mm2 - var2.inla.mean^2)
round(var2.inla.quant,3)
var2.inla.mean
var2.inla.sd

### Plot posterior
par(mfrow = c(1,1))
xaxis = seq(-3,3, by=0.001)
xlabel = c(expression(beta[0]), expression(beta[1]), expression(beta[2]),
           expression(beta[3]))
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
  abline(v = fixef(fpmm)[i], lty = 2)
  legend("topright", legend=c("LABPMM-grid", 
                              "LABPMM-MAP",
                              "INLA",
                              "MCMC"),
         col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)
}

# plot(epi.combined[,1], trace = FALSE, density=TRUE,
#      type = "l", 
#      ylim = c(0,3),
#      xlab = expression(beta[0]),
#      ylab = "posterior",
#      main = "")
# lines(xaxis, dnorm(xaxis, model2$fix_effect[1], model2$sd_fixEffect[1]), 
#       col = "blue",
#       type = "l")
# lines(xaxis, dnorm(xaxis, model1$fix_effect[1], model1$sd_fixEffect[1]), 
#       col = "red", type = "l", lty = 2)
# lines(model.inla$marginals.fixed[[1]], col = 'green', type = 'l', lty = 2)
# abline(v = fixef(fpmm)[1], lty = 2)
# legend("topright", legend=c("LABPMM-grid", 
#                             "LABPMM-MAP",
#                             "INLA",
#                             "MCMC"),
#        col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)
# 
# plot(epi.combined[,4], trace = FALSE, density=TRUE,
#      type = "l", 
#      ylim = c(0,1.2),
#      xlab = expression(beta[3]),
#      ylab = "posterior",
#      main = "")
# lines(xaxis, dnorm(xaxis, model2$fix_effect[4], model2$sd_fixEffect[4]), 
#       col = "blue",
#       type = "l")
# lines(xaxis, dnorm(xaxis, model1$fix_effect[4], model1$sd_fixEffect[4]), 
#       col = "red", type = "l", lty = 2)
# lines(model.inla$marginals.fixed[[4]], col = 'green', type = 'l', lty = 2)
# abline(v = fixef(fpmm)[4], lty = 2)
# legend("topright", legend=c("LABPMM-grid", 
#                             "LABPMM-MAP",
#                             "INLA",
#                             "MCMC"),
#        col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)
# 
# plot(epi.combined[,3], trace = FALSE, density=TRUE,
#      type = "l", 
#      ylim = c(0,1.4),
#      xlab = expression(beta[2]),
#      ylab = "posterior",
#      main = "")
# lines(xaxis, dnorm(xaxis, model2$fix_effect[3], model2$sd_fixEffect[3]), 
#       col = "blue",
#       type = "l")
# lines(xaxis, dnorm(xaxis, model1$fix_effect[3], model1$sd_fixEffect[3]), 
#       col = "red", type = "l", lty = 2)
# lines(model.inla$marginals.fixed[[3]], col = 'green', type = 'l', lty = 2)
# abline(v = fixef(fpmm)[3], lty = 2)
# legend("topright", legend=c("LABPMM-grid", 
#                             "LABPMM-MAP",
#                             "INLA",
#                             "MCMC"),
#        col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)

plot(epi.combined[,7], trace = FALSE, density=TRUE,
     type = "l", 
     ylim = c(0, 2),
     xlab = expression(sigma[0]^2),
     ylab = "posterior",
     main = "")
set.seed(12)
y2 <- 1/exp(rnorm(xaxis, model2$mean_v[1], model2$sd_v[1]))
lines(density(y2), 
      col = "blue",
      type = "l")
set.seed(12)
y1 <- 1/exp(rnorm(xaxis, model1$mean_v[1], model1$sd_v[1]))
lines(density(y1),
      col = "red", type = "l", lty = 2)
lines(marg.var1, col = 'green', type = 'l', lty = 2)
abline(v = as.data.frame(VarCorr(fpmm))$vcov[1], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)


plot(epi.combined[,8], trace = FALSE, density=TRUE,
     type = "l", 
     xlim = c(0,0.8),
     ylim = c(0, 14),
     xlab = expression(sigma[1]^2),
     ylab = "posterior",
     main = "")
set.seed(12)
y2 <- 1/exp(rnorm(xaxis, model2$mean_v[2], model2$sd_v[2]))
lines(density(y2), 
      col = "blue",
      type = "l")
set.seed(12)
y1 <- 1/exp(rnorm(xaxis, model1$mean_v[2], model1$sd_v[2]))
lines(density(y1),
      col = "red", type = "l", lty = 2)
lines(marg.var2, col = 'green', type = 'l', lty = 2)
abline(v = as.data.frame(VarCorr(fpmm))$vcov[2], lty = 2)
legend("topright", legend=c("LABPMM-grid", 
                            "LABPMM-MAP",
                            "INLA",
                            "MCMC"),
       col=c("blue", "red", "green", "black"), lty=c(1,2,2,1), cex=0.8)

library(microbenchmark)
model1_time <- function() labpmm2_ci(y=Total, t = t, r=r, numUnit=numUnit, obsPerUnit = obsPerUnit, 
                                     gammaA = 0.01, gammaB = 0.01)
model1a_time <- function() labpmm2(y=Total, t = t, r=r, numUnit=numUnit, obsPerUnit = obsPerUnit, 
                                   gammaA = 0.01, gammaB = 0.01)
model2_time <- function() labpmm2_grid(y=Total, t = t, r=r, numUnit=numUnit, obsPerUnit = obsPerUnit, 
                                       gammaA = 0.01, gammaB = 0.01)
model3_time <- function() glmer(Total ~ Settlement +Time + TimeSetInt +(1+Time|Housecode), 
                                data= data, family="poisson")
model4_time <- function() inla(Total ~ Settlement + Time + TimeSetInt + 
                                 + f(Houseid, model ="iid", hyper = prec.prior) + 
                                 f(slopeid, Time, model ="iid", hyper = prec.prior),
                               data = data, family = "poisson",
                               control.fixed = prior.fixed)
microbenchmark(model1_time(), model1a_time(), model2_time(), model3_time(),
               model4_time())


