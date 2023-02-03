get_hessian <- function(phi, inv.cov.mat){
  # return the hessian matrix of log likelihood of model given eta and data
  # phi: vector of fixed effects and random effects
  # inv.cov.mat: inverse of covariance matrix of phi which is known given eta
  t(t)%*%(t*c(exp(t%*%phi))) + inv.cov.mat
}
