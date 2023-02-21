rm(list = ls())
library(HSAUR3)

#load data
data = epilepsy
attach(epilepsy)

n = nrow(data)
numUnit = 59
obsPerUnit = 4
numFixEff = 6

x1 = rep(1,n)
x2 = log(base/4)
x2 = x2-mean(x2)
x3 = ifelse(treatment=="placebo", 0, 1)
x4 = log(base/4)*x3
x4 = x4-mean(x4)
x5 = log(age)
x5 = (x5-mean(x5))
x6 = ifelse(period=="4", 1, 0)

y=seizure.rate
t = cbind(x1,x2,x3,x4,x5,x6)
colnames(t)[1:6] = c("intercept", "base", "treatment", "base_treatment", "age", "fourth_visit")

data2 <- cbind(y,t[,1:6], subject)
data2 <- data.frame(data2)
### Spaghetti plot over time
library(ggplot2)

### 
ggplot(data = data, aes(x = factor(period), y = seizure.rate)) +       
  geom_line(aes(group = subject)) + geom_point()

### Frequentist Poisson Mixed Model 
library(Matrix)
library(lme4)
library(tictoc)

fpmm <- glmer(seizure.rate ~ x2+x3+x4+x5+x6 +(1|subject), family="poisson")
summary(fpmm)

## Get the random intercept
getME(fpmm, "b")

## Model 1: Simulation study
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

data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, num.unit = numUnit, obs.per.unit = obsPerUnit)
t = data[,-c(1:3)]
y = data[, 2]
data2 = data.frame(data)

library(Matrix)
library(lme4)
library(tictoc)

tic("Running time for likelihood method")
fpmm <- glmer(y ~ V4 + V5+ V6+V7+(1|group), data=data2, family="poisson")
toc()
summary(fpmm)

## Get the random intercept
getME(fpmm, "b")