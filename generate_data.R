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