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
      b[j] ~ dnorm(0, tau1)
    }
    for (i in 1:N){
      rate[i] ~ dpois(phi[i])
      e[i] ~ dnorm(0, tau2)
      log(phi[i]) <- beta0 + beta1*base[i] + beta2*treatment[i] + beta3*base_treatment[i]
                  + beta4*age[i] + beta5*fourth_visit[i] + b[subject[i]] + e[i]
    }
    tau1 ~ dgamma(0.001, 0.001)
    tau2 ~ dgamma(0.001, 0.001)
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
  list(beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, beta4=0.1,beta5=0.2, tau1=1, tau2=1),
  list(beta0=0.1,beta1=0.1, beta2=0.1,beta3=0.1, beta4=0.1,beta5=0.1, tau1=1, tau2=1),
  list(beta0=0.1,beta1=0.2, beta2=0.1,beta3=0.2, beta4=0.1,beta5=0.2, tau1=0.1, tau2=0.1)
)

#my.inits <- function()(list(beta0=rnorm(1,0.4,1), beta1=rnorm(1), tau1=runif(1,0.5,1), tau2=runif(1,0.5,1)))

# Specify parameters to monitor

parameters <- c("beta0","beta1","beta2","beta3","beta4","beta5","tau1", "tau2")

## Running JAGS:

jags <- jags.model(file="epilepsy.model.txt",
                   data = my.data,
                   inits = my.inits,
                   n.chains = 3)
update(jags,1000)
epi.sim <- coda.samples(model = jags,
                        variable.names = parameters,
                        n.iter=1500, 
                        thin=10)
# Produce general summary of obtained MCMC sampling

#print(epi.sim,digits=3)
plot(epi.sim)
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


