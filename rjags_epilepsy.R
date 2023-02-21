library("rjags")
library("coda")
library("readr")
library("dplyr")

library(HSAUR3)

#load data
data = epilepsy
attach(epilepsy)

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

y=seizure.rate
t = data.frame(cbind(x1, x2,x3,x4,x5,x6))
colnames(t)[1:6] = c("subject", "base", "treatment", "base_treatment", 
                     "age", "fourth_visit")

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
    tau ~ dgamma(0.5, 0.0164)
    sigma = 1/tau
    beta0 ~ dnorm(0,1.0E-4)
    beta1 ~ dnorm(0,1.0E-4)
    beta2 ~ dnorm(0,1.0E-4)
    beta3 ~ dnorm(0,1.0E-4)
    beta4 ~ dnorm(0,1.0E-4)
    beta5 ~ dnorm(0,1.0E-4)
  }", file="epilepsy.model.txt")

# Prepare data:

N <- length(y)
J <- n_distinct(t$subject)

my.data <- list(rate=y, base = t$base, treatment = t$treatment, 
                base_treatment = t$base_treatment, age = t$age, 
                fourth_visit = t$fourth_visit, N = N, J = J, subject = t$subject)

# Initial parameters: 3 lists for 3 chains

my.inits <- list(
  list(beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, beta4=0.1,beta5=0.2, tau=1),
  list(beta0=0.1,beta1=0.1, beta2=0.1,beta3=0.1, beta4=0.1,beta5=0.1, tau=1),
  list(beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, beta4=0.1,beta5=0.2, tau=0.1)
)

#my.inits <- function()(list(beta0=rnorm(1,0.4,1), beta1=rnorm(1), tau=runif(1,0.5,1)))

# Specify parameters to monitor

parameters <- c("beta0","beta1","beta2","beta3","beta4","beta5","tau","sigma")

## Running JAGS:

jags <- jags.model(file="epilepsy.model.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 3)
update(jags,100)
epi.sim <- coda.samples(model = jags,
                          variable.names = parameters,
                          n.iter=5000, 
                          thin=10)
# Produce general summary of obtained MCMC sampling

#print(epi.sim,digits=3)
#plot(epi.sim)
summary(epi.sim)

# Convert osteo.sim into mcmc.list for processing with CODA package

epi.mcmc <- as.mcmc.list(epi.sim)

# Produce general summary of obtained MCMC sampling

plot(epi.mcmc)
summary(epi.mcmc)


# Specific output obtained from CODA functions

par(mfrow=c(2,2))
traceplot(epi.mcmc)

cumuplot(epi.mcmc,ask=FALSE)
acfplot(epi.mcmc)
autocorr.plot(epi.mcmc)

par(mfrow=c(1,1))
crosscorr.plot(epi.mcmc)


par(mfrow=c(2,2))
densplot(epi.mcmc)
effectiveSize(epi.mcmc)
HPDinterval(epi.mcmc)

# Convergence tests

gelman.diag(epi.mcmc)
gelman.plot(epi.mcmc,ask=FALSE)

geweke.diag(epi.mcmc)
geweke.plot(epi.mcmc,ask=FALSE)


library(runjags)
epi.combined <- combine.mcmc(epi.mcmc)
HPDinterval(epi.combined)


# obtain DIC
dic <- dic.samples(model = jags,
                   n.iter=1500, 
                   thin=1)

### Simulations


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
  
  data = cbind(group, y, x)
  return(data)
}

fixEffect = c(1,2,3,4,5)
pRecision = 10
numUnit = 50
obsPerUnit = 6
numFixEff = length(fixEffect)
set.seed(2)
randomIntercept = rnorm(numUnit, 0, sd=sqrt(1/pRecision))

data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, num.unit = numUnit, obs.per.unit = obsPerUnit)
t = data.frame(data[,-c(1:3)])
y = data[, 2]
group = data[,1]

colnames(t)[1:4] = c("RV1", "RV2", "RV3", "RV4")

cat("model
  {
    for (j in 1:J){
      b[j] ~ dnorm(0, tau)
    }
    for (i in 1:N){
      rate[i] ~ dpois(phi[i])
      log(phi[i]) <- beta0 + beta1*RV1[i] + beta2*RV2[i] + beta3*RV3[i]
                  + beta4*RV4[i] + b[group[i]]
    }
    tau ~ dgamma(0.5,0.0164)
    sigma = 1/tau
    beta0 ~ dnorm(0,1.0E-4)
    beta1 ~ dnorm(0,1.0E-4)
    beta2 ~ dnorm(0,1.0E-4)
    beta3 ~ dnorm(0,1.0E-4)
    beta4 ~ dnorm(0,1.0E-4)
  }", file="simulation.model.txt")

# Prepare data:

N <- length(y)
J <- n_distinct(group)

my.data <- list(rate=y, RV1 = t$RV1, RV2 = t$RV2, 
                RV3 = t$RV3, RV4 = t$RV4, 
                N = N, J = J, group = group)

# Initial parameters: 3 lists for 3 chains

my.inits <- list(
  list(beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, beta4=0.1, tau=1),
  list(beta0=0.1,beta1=0.1, beta2=0.1,beta3=0.1, beta4=0.1, tau=1),
  list(beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, beta4=0.1, tau=0.1)
)

#my.inits <- function()(list(beta0=rnorm(1,0.4,1), beta1=rnorm(1), tau=runif(1,0.5,1)))

# Specify parameters to monitor

parameters <- c("beta0","beta1","beta2","beta3","beta4","tau","sigma")

## Running JAGS:

jags <- jags.model(file="simulation.model.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 3)
update(jags,100)
sim.sim <- coda.samples(model = jags,
                        variable.names = parameters,
                        n.iter=1500, 
                        thin=10)
# Produce general summary of obtained MCMC sampling

#print(epi.sim,digits=3)
#plot(epi.sim)
model_obj = summary(sim.sim)
model_obj
model_obj$statistics[,1]

# Convert osteo.sim into mcmc.list for processing with CODA package

sim.mcmc <- as.mcmc.list(sim.sim)

# Produce general summary of obtained MCMC sampling

plot(sim.mcmc)
summary(sim.mcmc)


# Specific output obtained from CODA functions

par(mfrow=c(2,2))
traceplot(sim.mcmc)

cumuplot(sim.mcmc,ask=FALSE)
acfplot(sim.mcmc)
autocorr.plot(sim.mcmc)

par(mfrow=c(1,1))
crosscorr.plot(sim.mcmc)


par(mfrow=c(2,2))
densplot(sim.mcmc)
effectiveSize(sim.mcmc)
HPDinterval(sim.mcmc)

# Convergence tests

gelman.diag(sim.mcmc)
gelman.plot(sim.mcmc,ask=FALSE)

geweke.diag(sim.mcmc)
geweke.plot(sim.mcmc,ask=FALSE)


library(runjags)
sim.combined <- combine.mcmc(sim.mcmc)
HPDinterval(sim.combined)


# obtain DIC
dic <- dic.samples(model = jags,
                   n.iter=1500, 
                   thin=1)



