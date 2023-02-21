run = 100
fixEffect = c(1,2,3,4,5)
pRecision = 10
numUnit = 50
obsPerUnit = 6
numFixEff = length(fixEffect)
set.seed(2)
randomIntercept = rnorm(numUnit, 0, sd=sqrt(1/pRecision))
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

sample = function(run){
  beta = matrix(0,100,5)
  for (i in 1:run){
    data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, num.unit = numUnit, obs.per.unit = obsPerUnit)
    t = data[,-c(1:3)]
    y = data[, 2]
    model <- lapmm(y=y, t=t, numUnit=numUnit, obsPerUnit=obsPerUnit, numFixEff=numFixEff,
                    gammaA = 0.5, gammaB = 0.0164, fixEffVar = 10^4)
    beta[i, ] = model$fix_effect
  }
  return(beta)
}
beta = sample(run = run)
beta_true = c(1,2,3,4,5)
beta_mean = colMeans(beta_mean)
bias = beta_mean - beta_true

res = sweep(beta, 2, beta_true)
mse = sqrt(colMeans(res^2))

eres = sweep(beta, 2, beta_mean)
ese = sqrt(colSums(eres^2)/99)

#frequentist
sample_glmer = function(run){
  beta = matrix(0,100,5)
  for (i in 1:run){
    data = generate_data(fix.effect = fixEffect, random.intercept = randomIntercept, num.unit = numUnit, obs.per.unit = obsPerUnit)
    t = data[,-c(1:3)]
    y = data[, 2]
    data2 = data.frame(data)
    fpmm <- glmer(y ~ V4 + V5+ V6+V7+(1|group), data=data2, family="poisson")
    beta[i, ] = fixef(fpmm)
  }
  return(beta)
}
beta_freq = sample_glmer(run = run)
beta_true = c(1,2,3,4,5)
beta_mean = colMeans(beta_freq)
bias = beta_mean - beta_true

res = sweep(beta_freq, 2, beta_true)
mse = sqrt(colMeans(res^2))

eres = sweep(beta_freq, 2, beta_mean)
ese = sqrt(colSums(eres^2)/99)


